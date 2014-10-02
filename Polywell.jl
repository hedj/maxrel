# Computes Polywell fields in cylindrical polar or Cartesian coordinates.
# Provision for departures from ideality.
# (c) John Hedditch, 2014

module Polywell
  using Elliptic
  importall Constants

  # A value sufficient close to zero.
  machine_precision = 2.3e-15


# Electric field due to a charged ring of given radius and charge in cylindrical coords.
# Arguments:
#   Q is the total charge on the ring (Coulombs)
#   R is the radius of the coil (metres)
#   z is the distance along the axial direction from the center of the ring (metres)
#   rho is the radial distance from the axis (metres)
# Returns:
#  [E_z, E_radial]    (units of volts / metre)
function current_loop_E_cyl(Q, R, z, rho)
  n = -4.0*R*rho / ((rho -R)^2 + z^2)
  #println("n = $n")
  rplus = (rho + R)^2 + z^2
  rminus = (rho - R)^2 + z^2
  tsq = rminus + 2.0*rho*R
  R_on_l = 1.0 / (2.0*rho)

  # If rho=0, use the limiting forms for accuracy
  if (rho < machine_precision) # ~ machine epsilon for double precision
    Ez = Q*z / (4.0 * pi * epsilon_0 * (z^2 + R^2)^1.5 )
    Erho = 0.0
  else 

    # Approximate K, E using power-series, as Elliptic Library doesn't support negative arguments
    K = 0.5 * pi*(1.0 + (x/4.0) + (9.0*x^2/64.0) + (25.0*x^3/256.0) + (1225.0*x^4/16384.0) )
    E = 0.5 * pi*(1.0 - (x/4.0) - (3.0*x^2/64.0) - (5.0*x^3/256.0) -(175.0*x^4/16384.0) )

    Ez = Q*z*sqrt(1.0 - n)*E / ( 2.0*pi^2*epsilon_0 * rplus^1.5 )
    Erho = ( Q / (2.0*pi^2*epsilon_0*sqrt(rminus)) ) * (  (R_on_l * K) +  (rho / rplus)*E - (R_on_l * tsq / rplus)* E  )
  end

  return [Ez, Erho]
end

# Magnetic field due to a current loop of given radius and current in cylindrical coords.
# Arguments:
#   I is the total current through the ring (Coulombs)
#   R is the radius of the coil (metres)
#   z is the distance along the axial direction from the center of the ring (metres)
#     positive z corresponds to direction of field at centre due to a positive current.
#   rho is the radial distance from the axis (metres)
# Returns:
#  [B_z, B_radial]    (units of Tesla)
function current_loop_B_cyl(I, R, z, rho)
  m = 4*R*rho / (z^2 + (R + rho)^2)

  # If rho=0, use the limiting forms.
  if (rho < machine_precision)
    Brho = 0.0
    Bz = 0.5 * mu_0 * I * R^2 / (( z^2 + R^2 )^1.5)
  else
    K,E = Elliptic.ellipke(m)
    scale_factor = mu_0 * I * sqrt(m / (R*rho^3)) / (4.0*pi) 
    Brho = scale_factor * z * ( (2.0 - m)*E/(2.0 - 2.0*m) - K )
    Bz = scale_factor * ( rho*K + E*(R*m - rho*(2.0 -m))/(2.0 - 2.0*m) )
  end

  return [Bz, Brho]
end

# Helper function - convert from cartesians to 
# local cylindrical coords (z, rho) with respect to a given point and direction
# Arguments:
#   x,y,z         : coordinates to convert
#   x_c, y_c, z_c : origin of cylindrical coordinate system
#   x_n, y_n, z_n : direction of cylindrical coordinate z-axis
function cart_to_cyl( x,y,z,
                      x_c,y_c,z_c,
                      x_n,y_n,z_n )

  r_vec = [x-x_c, y-y_c, z-z_c]
  zcyl = dot(r_vec, [x_n, y_n, z_n])
  rho = sqrt( dot(r_vec, r_vec) - zcyl^2 )

  rho_vec = r_vec - zcyl*[x_n, y_n, z_n]
  if ( dot(rho_vec,rho_vec) > machine_precision )
    rho_hat = rho_vec / sqrt(dot(rho_vec, rho_vec))
  else
    # on-axis; no meaningful radial normal direction exists.
    rho_hat = [0.0, 0.0, 0.0]
  end

  return (zcyl, rho, rho_hat)
end 

# Compute the fields at a given position in space for a simple polywell
# centred at the origin, with coils of radius R, at plus and minus R + a
# charge Q and current I through each coil.
# Arguments:
# Returns: B, E
#   B = [Bx, By, Bz],    (units of Tesla)
#   E = [Ex, Ey, Ez]     (units of V/m) 

function polyBE(I, Q, R, a, x, y, z)
  
  locations = ( [ -(R+a),    0.0,    0.0 ],
                [  (R+a),    0.0,    0.0 ],
                [    0.0, -(R+a),    0.0 ],
                [    0.0,  (R+a),    0.0 ],
                [    0.0,    0.0, -(R+a) ],
                [    0.0,    0.0,  (R+a) ]  )


  normals =   ( [  1.0,  0.0,  0.0 ], 
                [ -1.0,  0.0,  0.0 ],
                [  0.0,  1.0,  0.0 ],
                [  0.0, -1.0,  0.0 ],
                [  0.0,  0.0,  1.0 ],
                [  0.0,  0.0, -1.0 ] )


  # Just for paranoia, make sure we have normal direction for each centre.
  @assert length(locations) == length(normals)

  B = zero_vec
  E = zero_vec

  for i=1:length(locations) 
    (x_c,y_c,z_c) = locations[i]
    (x_n,y_n,z_n) = normals[i]
    zcyl, rho, rho_hat = cart_to_cyl(x,y,z, x_c,y_c,z_c, x_n,y_n,z_n)
    (Baxial, Brho) = current_loop_B_cyl(I, R, zcyl, rho)
    (Eaxial, Erho) = current_loop_E_cyl(Q, R, z, rho)

    B += Baxial * [x_n, y_n, z_n] + Brho.*rho_hat
    E += Eaxial * [x_n, y_n, z_n] + Erho.*rho_hat

  end

  return (B,E)
end

end # module Polywell
