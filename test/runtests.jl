using Hop, Base.Test

lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]

graphene = TightBindingModel(lat, positions)

@test graphene.norbits == 2

@test graphene.rlat ≈ [
    2π 0.0 0.0;
    -2π*tan(π/6) 2π/cos(π/6) 0.0;
    0.0 0.0 2π;
]

sethopping!(graphene, 1, 2, [0, 0, 0], -1.0)
sethopping!(graphene, 2, 1, [1, 0, 0], -1.0)
sethopping!(graphene, 2, 1, [0, 1, 0], -1.0)

@test caleig(graphene, [0.0, 0.0, 0.0]) ≈ [-3.0, 3.0]
@test calhamiltonian(graphene, [0.0, 0.0, 0.0]) * caleig(graphene, [0.0, 0.0, 0.0], true)[2][:, 1] ≈ -3.0*caleig(graphene, [0.0, 0.0, 0.0], true)[2][:, 1]

graphene2x2 = makesupercell(graphene, [2, 2, 1])

@test graphene2x2.norbits == 8
@test caleig(graphene2x2, [0.0, 0.0, 0.0]) ≈ sort!(
    [
        caleig(graphene, [0.0, 0.0, 0.0]);
        caleig(graphene, [0.5, 0.0, 0.0]);
        caleig(graphene, [0.0, 0.5, 0.0]);
        caleig(graphene, [0.0, 0.5, 0.0])
    ]
)

grapheneisolated = makecluster(graphene)

@test haskey(grapheneisolated.hoppings, [1, 2, 0, 0, 0])
@test !haskey(grapheneisolated.hoppings, [2, 1, 1, 0, 0])
