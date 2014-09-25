# Simulate electron cusp losses in a Polywell. 
# (c) John Hedditch, 2014
module CuspConfinement
  using MaxRel
  using Polywell
  using Roots
  importall Constants

  # for determinism, seed the RNG with the same value each time
  function seedRNG()
    srand(1010283)
  end

  function CylindricalBeamAlongXAxis(R,nparticles,v,radius)
    Particles::Vector{MaxRel.Particle} = []
    for i=1:nparticles
      # generate a particle just inside the trap (at x = -R)
      while(true)
        y = rand()
        z = rand()
        if (y^2 + z^2) <= 1.0
          Position = [ -R, 2.0*radius*y - radius, 2.0*radius*z - radius ]        
          # 1% scatter in velocities
          Scatter = v * 0.01 * (rand(3) - 0.5)
          Velocity = [ v, 0.0, 0.0 ] + Scatter
          append!(Particles, [MaxRel.Electron(Position, Velocity, zero_vec)])
          break;
        end
      end
    end
  end

  function IntersectingCylindricalBeamsAlongXAxis(R,nparticles,v,radius)

    xscatter = 0.01

    Particles::Vector{MaxRel.Particle} = []
    for i=1:nparticles
      # generate a particle just inside the trap (at x = -R)
      while(true)
        y = rand()
        z = rand()
        if (y^2 + z^2) <= 1.0
          Position = [ -R + xscatter*(rand()-0.5), 2.0*radius*y - radius, 2.0*radius*z - radius ]        
          # 1% scatter in velocities
          Scatter = v * 0.01 * (rand(3) - 0.5)
          Velocity = [ v, 0.0, 0.0 ] + Scatter
          append!(Particles, [MaxRel.Electron(Position, Velocity, zero_vec)])
          break;
        end
      end
    end

    for i=1:nparticles
      # generate a particle just inside the trap (at x = R)
      while(true)
        y = rand()
        z = rand()
        if (y^2 + z^2) <= 1.0
          Position = [ R + xscatter*(rand()-0.5), 2.0*radius*y - radius, 2.0*radius*z - radius ]        
          # 1% scatter in velocities
          Scatter = v * 0.01 * (rand(3) - 0.5)
          Velocity = [ -v, 0.0, 0.0 ] + Scatter
          append!(Particles, [MaxRel.Electron(Position, Velocity, zero_vec)])
          break;
        end
      end
    end

    return Particles
  end

  function ContrivedTestCase(R, v)
    Particles::Vector{MaxRel.Particle} = []
    append!(Particles, [MaxRel.Electron( [ -R, 0.01, 0.0 ], [v, 0.0, 0.0], zero_vec ) ])
    append!(Particles, [MaxRel.Electron( [ R, 0.01, 0.0 ], [-v, 0.0, 0.0], zero_vec ) ])
    return Particles
  end

  function LogParticles(Particles, f)
    write(f, "x, v, a, m, q\n")
    for i=1:length(Particles)
      p = Particles[i]
      loc = p.loc[1]
      x = string(loc.x)
      v = string(loc.v)
      a = string(loc.a)
      m = string(p.m)
      q = string(p.q)
      write(f, "$x,$v,$a,$m,$q\n")
    end
  end

  function LogLocations(Particles, f, t)
    nparticles = length(Particles)
    write(f, "$nparticles, $t\n")
    for i=1:nparticles
      x = string(Particles[i].loc[1].x)
      #v = string(Particles[i].loc[1].v)
      #a = string(Particles[i].loc[1].a)
      write(f, "$x\n")
    end
  end

  function RemoveEscapedParticles(Particles, escape_radius)
    NewParticles::Vector{MaxRel.Particle} = []
    n_lost = 0
    for i=1:length(Particles)
      p = Particles[i]
      if ( dot(p.loc[1].x, p.loc[1].x) > escape_radius^2 )
        n_lost += 1
      else
        append!(NewParticles, [Particles[i]])
      end
    end
   
    return (n_lost, NewParticles)
  end

  function ComputeStepSize(Particles, upper_bound)
    # Make sure time step * velocity <  1/2 smallest interparticle distance
    min_dsq = 1e34
    max_v = 0.0
    for i=1:length(Particles)
      for j=i+1:length(Particles)
        R = Particles[j].loc[1].x - Particles[i].loc[1].x
        min_dsq = min(min_dsq, dot(R,R))
      end
      max_v = max(max_v, sqrt(dot(Particles[i].loc[1].v, Particles[i].loc[1].v)))
    end
    step = min(upper_bound, 0.5 * sqrt(min_dsq) / max_v)
    return step
  end

  function ev_to_metres_per_second(electron_volts, charge, mass)
    function f(v)
      return v^2*(1.0 / sqrt(1.0 - (v^2/c^2))) - (2.0*abs(charge)*electron_volts/mass)
    end
    return Roots.fzero(f,0.0, Constants.c)
  end

  # Compute intersecting beams with given field strength and electron energy etc.
  function DemoScenario(Bmax, nparticles, R, a, electron_volts)
    I = 2*R*Bmax / mu_0
    v = ev_to_metres_per_second(electron_volts, q_e, m_e)
    MeasureCuspLoss(R, a, I, 0.0, v, nparticles, 1.0e-9, 1e-7, "results/confinement")
  end

  function MeasureCuspLoss(R,a,I,Q,v,nparticles,max_step,runtime,output_prefix,log_interval=1e-10)

    seedRNG()
    # set up particles in a 5mm-radius beam
    Particles = IntersectingCylindricalBeamsAlongXAxis(R, nparticles, v, 5e-3)
    #Particles = ContrivedTestCase(R, v)

    # Open files, record simulation parameters
    paramsfile = open("$output_prefix.params", "w")
    write(paramsfile, "--- Parameter Summary ---\n")
    write(paramsfile, "R: $R\n")
    write(paramsfile, "a: $a\n")
    write(paramsfile, "I: $I\n")
    write(paramsfile, "Q: $Q\n")
    write(paramsfile, "nparticles: $nparticles\n")
    write(paramsfile, "max_step: $max_step\n")
    write(paramsfile, "runtime: $runtime\n")
    close(paramsfile)

    particle_file = open("$output_prefix.particles", "w")
    LogParticles(Particles, particle_file)
    close(particle_file)

    trajectory_file = open("$output_prefix.trajectories", "w")
    loss_file = open("$output_prefix.loss", "w")
    write(loss_file, "t,             n_lost\n")

    function external_field(pos)
      (x,y,z) = pos
      return Polywell.polyBE(I, Q, R, a, x, y, z)
    end

    # Begin! 
    t = 0.0
    last_log = 0.0
    while(t < runtime)
      if ( (t -last_log) >= log_interval)
        last_log = t
        LogLocations(Particles, trajectory_file, t)
      end
      step = ComputeStepSize(Particles, max_step)
      MaxRel.green("t = $t")
      # Record positions to disk every nsteps_per_log steps
       #(measure number of particles 'lost'), remove from simulation
       #(we take 'loss' as escaping to a distance of 2R from the origin)
      (n_lost, Particles) = RemoveEscapedParticles(Particles, 2.0*R)
      if (n_lost > 0)
        println("Lost $n_lost particles")
        write(loss_file, "$t,  $n_lost\n")
      end
      # do update ( 10ns max history length )
      NewParticles = MaxRel.nbody(Particles, t, step, external_field, 1.0e-8)
      Particles = deepcopy(NewParticles)
      t += step
    end
    close(loss_file)
    close(trajectory_file)

    # Summarise
  end

end 
