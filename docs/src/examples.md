# Examples

## Graphene Bands

```julia
using Hop
using Plots
# Lattice vector. Notice that lattice vectors are stored by column.
lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
# Reduced positions of orbits. Notice that positions are stored by column.
positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]

graphene = TightBindingModel(lat, positions)

sethopping!(graphene, 1, 2, [0, 0, 0], -1.0) # ⟨1|H|(0, 0, 0)2⟩=-1
sethopping!(graphene, 2, 1, [1, 0, 0], -1.0) # ⟨2|H|(1, 0, 0)1⟩=-1
sethopping!(graphene, 2, 1, [0, 1, 0], -1.0) # ⟨2|H|(0, 1, 0)1⟩=-1

kdist, egvals = calband(graphene, [1 0; 0 1; 0 0], 100)

plot(kdist, egvals', label=["",""], lw=2, color=:black, dpi=400,
    xaxis=(raw"$k$", font(15, "sans-serif")),
    yaxis=("E(eV)", font(15, "sans-serif")),
    xticks=[],
    size=(400, 300)
    )
savefig("bands.png")
```
output figure:

![graphene](https://i.imgur.com/E6Yg6Nx.png)
