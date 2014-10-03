# A relativistic n-body Maxwell code
# (c) John Hedditch, 2014
module MaxRel
using NLopt
importall Constants
using Base.Test

        custom_handler(r::Test.Success) = green("Success on $(r.expr)")
        custom_handler(r::Test.Failure) = error("test failed: $(r.expr)")
        custom_handler(r::Test.Error) = rethrow(r)


        # A relativistic N-Body code
        # Uses the retarded potential formulation to compute exact
        # evolution of particles under external and self-fields.

        # This is a dimensional code, using SI units.
	
        export Coordinate, Particle, selftest, green, red


	type Coordinate
	  x::Vector{Float64}  # at full times
          v::Vector{Float64}  # at half-times
	  a::Vector{Float64}  # at full times
          t::Float64
	end

	type Particle
	  m::Float64 # mass
	  q::Float64 # charge 
	  loc::Vector{Coordinate} 
	end 

        function Electron(pos, vel, acc)
          return Particle(m_e, q_e, [ Coordinate(pos, vel, acc, 0.0)])
        end

        function TestParticle(pos, vel, acc)
          return Particle(1, 1, [ Coordinate(pos, vel, acc, 0.0)])
        end

        function red(str) 
          Base.print_with_color(:red, str, "\n") 
        end

        function redn(str) 
          Base.print_with_color(:red, str) 
        end
        
        function green(str)
          Base.print_with_color(:green, str, "\n")
        end

        function selftest()
          tic()
          # test that our integration step works as expected.
          print("Testing electrostatics: ")
          test_electrostatic_particle_evolution(1000)
          green("PASS")
 
          print("Testing Electric field due to particle: ")
          test_particle_electric_field() 
          green("PASS")

          print("Testing Magnetic field due to particle: ")
          test_particle_magnetic_field() 
          green("PASS")

          print("Testing magnetostatics: ")
          test_magnetostatic_particle_evolution(1e7, 2000)
          green("PASS")
          toc();
        end

        function test_particle_electric_field()
          # Check Coulomb repulsion
          p1 = TestParticle(zero_vec, zero_vec, zero_vec) 
          p2 = TestParticle([1.0, 0.0, 0.0], zero_vec, zero_vec)
          p3 = TestParticle([2.0, 0.0, 0.0], zero_vec, zero_vec)

          # E as measured by p1 should be [ - 1 / (4*pi*epsilon_0), 0, 0]
          expected_E = [-1.0 / (4*pi*epsilon_0), 0.0, 0.0]
          err = E(p1, p2, 0.0, 1.0) - expected_E
          @Test.test dot(err,err) < 1.0e-14

          # E as measured by p2 should be [ 1 / (4*pi*epsilon_0), 0, 0]
          expected_E = [1.0 / (4*pi*epsilon_0), 0.0, 0.0]
          err = E(p2, p1, 0.0, 1.0) - expected_E
          @Test.test dot(err,err) < 1.0e-14

          # E due to p1 as measured by p3 should be [ 1 / (16*pi*epsilon_0), 0, 0]
          expected_E = [1.0 / (16*pi*epsilon_0), 0.0, 0.0]
          err = E(p3, p1, 0.0, 1.0) - expected_E
          @Test.test dot(err,err) < 1.0e-14


           # E as measured by p4 due to p5 should be [  q_e / (4*pi*epsilon_0), 0, 0]
           p4 = Electron(zero_vec, zero_vec, zero_vec)
           p5 = Electron([-1.0, 0.0, 0.0], zero_vec, zero_vec)
           expected_E = [ q_e / (4*pi*epsilon_0), 0.0, 0.0]
           err = E(p4, p5, 0.0, 1.0) - expected_E
           @Test.test dot(err,err) < 1.0e-14

          # Should see no magnetic field.
          @Test.test dot(B(p1,p2,0.0, 0.0), B(p1,p2,0.0,0.0)) < 1.0e-6
          @Test.test dot(B(p2,p1,0.0, 0.0), B(p1,p2,0.0,0.0)) < 1.0e-6
        end

        function test_particle_magnetic_field()
          p1 = TestParticle(zero_vec, [0.0, 1.0, 0.0], zero_vec) 
          p2 = TestParticle([1.0, 0.0, 0.0], zero_vec, zero_vec)

          # E as measured by p1 should be [ - 1 / (4*pi*epsilon_0), 0, 0]
          expected_E = [-1.0 / (4*pi*epsilon_0), 0.0, 0.0]
          err = E(p1, p2, 0.0, 1.0) - expected_E
          @Test.test dot(err,err) < 1.0e-13

          # B as measured by p2 should be (mu_0 / 4*pi) * [ 0, 0, -1 ]
          expected_B = (mu_0 / (4.0*pi)) * [0.0, 0.0, -1.0]
          found_B = B(p2, p1, 0.0, 1.0)
          err = found_B - expected_B
          @Test.test dot(err, err) < 1.0e-13

        end

        function test_electrostatic_particle_evolution(nsteps)

          # Set up a particle of unit charge-to-mass at rest at the origin
          step = 1e-3 / nsteps
          # 1 - (1/2)*a*t^2 where a = 1 and t^2 = (half * step)^2
          vinit = 1.0 - 0.125 * step * step

          p = TestParticle(zero_vec, zero_vec, [1.0, 0.0, 0.0])
          # Set up x-axis aligned electric field of 1V/m 
          E = [1.0, 0.0, 0.0]
          # No magnetic field
          B = [0.0, 0.0, 0.0]

          # Compute velocity after 1ms, with 0.1 microsecond timesteps
          # acceleration should be arbitrarily close to 1m/s^2, so end velocity
          # should be 1e-3 m/s  (relativistic corrections less than machine precision)
          # position should be 1/2 * t^2 = 0.5e-6
          t = 0.0
          for i=1:nsteps
            p1prime = do_update(p, E, B, step)
            p1prime.t = i*step
            # store new result
            prepend!(p.loc, [p1prime])
          end 
  
          v = p.loc[1].v[1]
          x = p.loc[1].x[1]
          @Test.test abs(p.loc[1].v[1] - 1e-3) < 1e-7
          @Test.test abs(p.loc[1].x[1] - 0.5e-6) < 1e-7 

        end

        function test_magnetostatic_particle_evolution(vy, nsteps)

          # Check that we produce the correct gyroradius and gyrofrequency for a test particle
          # of unit charge to mass with given y-velocity in unit z-directed magnetic field.

          # No electric field
          zero_E = [0.0, 0.0, 0.0]
          # Unit magnetic field
          Bz = [0.0, 0.0, 1.0]


          # relativistic case
          r_g = gamma(vy) * vy

          omega_g = 1 / gamma(vy) 

          step = (2*pi*gamma(vy)) / nsteps;

          vinit = vy * [sin(-0.5*step * omega_g ), cos( 0.5 * step * omega_g ), 0 ]
          gam = gamma(vinit)
          p = TestParticle([-r_g, 0.0, 0.0], vinit, zero_vec)
          t = 0.0
          
          for i=1:nsteps
            p1prime = do_update(p, zero_E, Bz, step)
            p1prime.t = i*step
            # store new result
            prepend!(p.loc, [p1prime])
          end 

          x = p.loc[1].x
          err = abs(x + r_g)
           
          @Test.test abs(x[1] + r_g) < 1.0/nsteps
          
        end

        function get_retarded_position(p, t)

          # find the last index where the time was
          # greater than the desired time t
          le_index = 0
          for i=1:length(p.loc) # going backwards in time
            if (p.loc[i].t >= t) then
              le_index = le_index + 1 
            else
              break
            end
          end

          # if no time recorded is greater than time t, then use the
          # last available time
          if( le_index == 0 ) then
            return deepcopy(p.loc[1])
          end
 
          # Otherwise
          if (le_index < length(p.loc)) then     # interpolate between le_index+1 and le_index
            p1 = p.loc[le_index+1]
            p2 = p.loc[le_index]
            if ((t - p1.t) > 0) then
              alpha = (t - p1.t) / (p2.t - p1.t)
            else # avoid trouble in the case that p2.t is v v close to p1.t
              alpha = 0.0
            end
            result = deepcopy(p1)
            result.x = alpha * p2.x + (1-alpha)*p1.x
            result.v = alpha * p2.v + (1-alpha)*p1.v
            result.a = alpha * p2.a + (1-alpha)*p1.a
            result.t = t
          else # we found no such index! Extrapolate - it is the best we can do.
            p_last = p.loc[le_index]
            result = deepcopy(p_last)
            result.x = p_last.x - ((p_last.t - t)*p_last.v)
            result.t = t
          end

          return result
        end


        function gamma(v)
          beta = v * invc;
          return 1.0 / sqrt( 1.0 - dot(beta,beta)) 
        end

	function E(p1,p2,t, step) # electric field due to p2 as seen by p1 at time t
          d = p1.loc[1].x - p2.loc[1].x
          delta_t = sqrt(dot(d,d)) / c  # how far back in time do we need to look?
          retarded_position = get_retarded_position(p2, t - delta_t)
          retarded_velocity = get_retarded_position(p2, t + 0.5*step - delta_t)
          R = p1.loc[1].x - retarded_position.x
          v = retarded_velocity.v
          a = retarded_position.a
          rsq = dot(R,R)
          r = sqrt(rsq)
          kappa = 1 - invc*(dot(v,R) / r)
          vsq = dot(v, v)

          Rhat = R / r

          return (p2.q / (4.0 * pi * epsilon_0)) * (
            Rhat / (kappa*rsq) 
            - (v / (kappa^2 * rsq * c)) 
            - (R / (kappa^3 * r^3 * c^2))*( (r*c + 2*kappa*c)*dot(v,Rhat) + r*(vsq -dot(a,R)) )
            - a / (kappa^2 * r * c)
            - (v / (kappa^3 * r^2 * c^2))*( dot(v, Rhat) - invc*(vsq - dot(a,R)) )
          )

	end

	function B(p1,p2,t, step) # magnetic field due to p2 as seen by p1 at time t
          d = p1.loc[1].x - p2.loc[1].x
          delta_t = sqrt(dot(d,d)) / c  # how far back in time do we need to look?
          retarded_position = get_retarded_position(p2, t - delta_t)
          retarded_velocity = get_retarded_position(p2, t + 0.5*step - delta_t)
          R = p1.loc[1].x - retarded_position.x
          v = retarded_velocity.v
          a = retarded_position.a
          rsq = dot(R,R)
          r = sqrt(rsq)
          kappa = 1 - invc*(dot(v,R) / r)
          vsq = dot(v, v)

          Rhat = R / r
          #moved = p2.loc[1].x - retarded_position.x
          #println("moved = $moved")
          #println("d = $d")
          #println("R = $R")
          #println("kappa = $kappa")
          #println("Rhat = $Rhat")
          #println("rsq = $rsq")

          return (p2.q * mu_0  / (4.0 * pi)) * (
            cross(v, Rhat) / (kappa*rsq)
            + (1.0 / (kappa*c))* (
                (cross(a, Rhat)/(kappa*r))
                + (cross(v, Rhat)/(kappa^2*rsq))* (      
                  dot(v, Rhat) - invc*(vsq - dot(a,R))
                ) 
              )
          )
	end

        # Perform one n-body iteration.
        # Particles : array of Particle objects
        #         t : time at which to perform the calculation
        #      step : time step 
        #      tmax : length of history to keep for each particle
        function nbody(Particles, t, step, external_field, tmax = 1.0e-8)
          NewParticles::Vector{Particle} = []
          n = length(Particles)
          for i=1:n 
            # compute E and B fields acting on particle i
            p1 = deepcopy(Particles[i])
            (extB,extE) = external_field(p1.loc[1].x)
            self_E = zero_vec
            self_B = zero_vec

            for j=1:n
              if (j!=i)
                self_E += E(p1, Particles[j], t, step)         
                self_B += B(p1, Particles[j], t, step)         
              end
            end

            # perform coordinate update 
            p1prime = do_update(p1, extE + self_E, extB + self_B, step)
            p1prime.t = t + step
            # store new result
            oldest = last(p1.loc)
            if ( t - oldest.t ) > tmax
              p1.loc = p1.loc[1:length(p1.loc)-1]
            end
            prepend!(p1.loc, [p1prime])
            append!(NewParticles, [p1])
          end

          return NewParticles
        end

        function vsqminuscsq(x::Vector, grad::Vector)
          return dot(x,x) - c^2
        end

        function do_update(p1, E, B, step) 
          # Perform an implicit update of the particle location.

          count::Int = 0
          charge_to_mass = (p1.q / p1.m)
          # compute acceleration
          function a(v::Vector)
            return (charge_to_mass / gamma(v))*(E + cross(v, B) - invcsquared*dot(E,v)*v)
          end

          vm = p1.loc[1].v  # v^{n-1/2}
          am = a(vm)
         
          # the function we wish to minimise, as a function of x == v^{n+1/2}
          # We do not calculate all the derivatives here, so we are restricted
          # to using derivative-free solvers
          function f(x::Vector, grad::Vector)
            count += 1
            if (dot(x,x)) > c^2
              return 1e6
            end
            err = (x - vm) - 0.5*(am + a(x))*step
            return dot(err,err)
          end 

          # use NLopt to find v^{n+1/2} as close as possible to the zero of f.
          opt = Opt(:LN_COBYLA, 3)
          min_objective!(opt, f)
          #ftol_rel!(opt, 1e-12)
          #xtol_rel!(opt, 1e-12)
          lower_bounds!(opt, [-c, -c, -c])
          upper_bounds!(opt, [c, c, c])
          inequality_constraint!(opt,vsqminuscsq, 1e-8)
          #println("step = $step")

          # compute the (constant-gamma) Boris estimate as a starting point
          invgvm = 1.0 / gamma(vm)
          E_eff = E - invcsquared*dot(E,vm)*vm
          vmb = vm + (E_eff * 0.5 * step * charge_to_mass)*invgvm
          t = (charge_to_mass * 0.5 * B * step) *invgvm
          s = (2.0 *  t) / (1.0 + dot(t,t))
          vprb = vmb + cross(vmb, t)
          vpb = vmb + cross(vprb, s)
          E_eff_final = E - invcsquared*dot(E,vpb)*vpb
          v_guess = vpb + (E_eff_final * 0.5 * step * charge_to_mass)*invgvm

          if ( f(v_guess, [0.0, 0.0, 0.0]) > 1e-6 ) then
            (minf, minx, ret) = optimize(opt, v_guess)
          else
            minx = v_guess
          end

          x = p1.loc[1].x
          xn = x + minx*step
          retval = Coordinate(xn, minx, (minx - vm)/step, 0)

          return retval

        end
end
