# Examples

## Graphene Bands

```julia
using Hop
using PyPlot
rc("font", family="Times New Roman")
rc("mathtext", fontset="cm")
rc("font", size=10)
# Lattice vector. Notice that lattice vectors are stored by column.
lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
# Reduced positions of orbits. Notice that positions are stored by column.
positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]

graphene = TightBindingModel(lat, positions)

sethopping!(graphene, 1, 2, [0, 0, 0], -1.0) # ⟨1|H|(0, 0, 0)2⟩=-1
sethopping!(graphene, 2, 1, [1, 0, 0], -1.0) # ⟨2|H|(1, 0, 0)1⟩=-1
sethopping!(graphene, 2, 1, [0, 1, 0], -1.0) # ⟨2|H|(0, 1, 0)1⟩=-1

kdist, egvals = getband(graphene, [1 0; 0 1; 0 0], 100)
figure(figsize=(3, 2.5))
plot(kdist, egvals[1, :], "k")
plot(kdist, egvals[2, :], "k")
xlabel("")
xticks([])
ylabel("E(eV)")
tight_layout()
savefig("band.png", dpi=500)
```
output figure:

<img src="https://i.imgur.com/PuCVnmC.png" style="width:400px;"/>
