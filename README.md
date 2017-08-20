# Hop.jl
[![Build Status](https://travis-ci.org/mistguy/Hop.jl.svg?branch=master)](https://travis-ci.org/mistguy/Hop.jl)
[![codecov](https://codecov.io/gh/mistguy/Hop.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mistguy/Hop.jl)

Hop.jl is a simple julia tight binding library.

# Examples
## Hofstadter Butterfly

```julia
"""
We are using square lattice to exemplify the Hofstadter Butterfly.
"""
using Hop
using Plots

# Size of the lattice is 15x15
sz = 15
# Magnetic flux is varied between 0 flux quantum and 1 flux quantum by 51 divisions.
nmags = 51
# Lattice vector. Notice that lattice vectors are stored by column.
lat = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
# For a square lattice, there are only one atom in one unit cell.
positions = Array{Float64}(3, 1)
positions[:, 1] = [0.0, 0.0, 0.0]
# Create a tight binding model.
t = TightBindingModel(lat, positions)
# Add hoppings between neighbourhood unit cells
sethopping!(t, 1, 1, [1, 0, 0], 1.0) # ⟨0|H|(1, 0, 0)⟩=1
sethopping!(t, 1, 1, [0, 1, 0], 1.0) # ⟨0|H|(0, 1, 0)⟩=1
# Create a supercell and then create a cluster from the supercell
c = makecluster(makesupercell(t, [sz, sz, 1]))
# List of magnetic field.
Bs = linspace(0, 1, nmags)
# Eigenvalues.
egvals = []
for i in 1:nmags
    hofstadter = deepcopy(c)
    # Add magnetic field. Laudau gauge is used internally.
    addmagneticfield!(hofstadter, Bs[i])
    # Calculate eigenvalues.
    push!(egvals, caleig(hofstadter, [0.0, 0.0, 0.0]))
end
# Plot the Butterfly
gr()
p = plot(size=(2000,2000))
for i in 1:nmags
    B = zeros(size(egvals[1]))
    fill!(B, B_list[i])
    plot!(B, egvals[i], seriestype=:scatter, markersize=1, markercolor=:black,
          legend=false, markeralpha=0.05, axis=nothing)
end
savefig("butterfly.png")
```

output figure (Notice the following butterfly is produced with a 100x100 lattice):

![butterfly](http://i.imgur.com/IBzSsXV.png)
