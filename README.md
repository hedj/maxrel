maxrel
======

A relativistic N-body code.

# Bad
 - It is very young and almost certainly has lurking bugs.
 - It also has a lot of room for optimization.

# Good
 - It comes with code for computing E and B fields of coils and polywells
 - It comes with tests
 - It computes properly relativistic and time-delayed fields.

# Using maxrel
 - Install [Julia](https://julialang.org)
 - Install the deps
 -- [Elliptic](https://github.com/nolta/Elliptic.jl)
 -- [NLopt](https://github.com/JuliaOpt/NLopt.jl)
 -- [Roots](https://github.com/JuliaLang/Roots.jl)
 -- (The above can be done using 'Pkg.update(); Pkg.add("Elliptic"); Pkg.add("NLopt"); Pkg.add("Roots")
 - See e.g CuspConfinment.jl for a simple usage example.
# Contact
j.hedditch@gmail.com (John Hedditch)
