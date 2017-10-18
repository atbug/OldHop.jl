# Examples

## Graphene Bands

```julia
using Hop
using Pyplot
# Lattice vector. Notice that lattice vectors are stored by column.
lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
# Reduced positions of orbits. Notice that positions are stored by column.
positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]

graphene = TightBindingModel(lat, positions)

sethopping!(graphene, 1, 2, [0, 0, 0], -1.0) # ⟨1|H|(0, 0, 0)2⟩=-1
sethopping!(graphene, 2, 1, [1, 0, 0], -1.0) # ⟨2|H|(1, 0, 0)1⟩=-1
sethopping!(graphene, 2, 1, [0, 1, 0], -1.0) # ⟨2|H|(0, 1, 0)1⟩=-1

kdist, egvals = calband(graphene, [1 0; 0 1; 0 0], 3)

plot(kdist, egvals[1, :])
plot(kdist, egvals[2, :])
savfig("bands.png")
```
output figure:

![graphene](https://i.imgur.com/NowwGtr.png)


## Hofstadter Butterfly

```julia
"""
We are using square lattice to exemplify the Hofstadter Butterfly.
"""
using Hop
using Plots

# Size of the lattice is 15x15
sz = 15
# Magnetic flux is varied between 0 flux quantum and 1 flux quantum
# by 51 divisions.
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
c = cutedge(cutedge(cutedge(
    makesupercell(t, [sz 0 0; 0 sz 0; 0 0 1]), 1), 2), 3)
# List of magnetic field.
Bs = linspace(0, 1, nmags)
# Eigenvalues.
egvals = []
for i in 1:nmags
    # Add magnetic field. Laudau gauge is used internally.
    hofstadter = addmagneticfield(c, Bs[i])
    # Calculate eigenvalues.
    push!(egvals, caleig(hofstadter, [0.0, 0.0, 0.0]))
end
# Plot the Butterfly
gr()
p = plot(size=(2000,2000))
for i in 1:nmags
    B = zeros(size(egvals[1]))
    fill!(B, Bs[i])
    plot!(B, egvals[i], seriestype=:scatter, markersize=1, markercolor=:black,
          legend=false, markeralpha=0.05, axis=nothing)
end
savefig("butterfly.png")
```

output figure:
(Notice the following butterfly is produced with a 100x100 lattice)

![butterfly](http://i.imgur.com/IBzSsXV.png)
