using Hop, Test

let
    lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
    positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]
    graphene = TightBindingModel(lat, positions)
    sethopping!(graphene, 1, 2, [0, 0, 0], -1.0)
    sethopping!(graphene, 2, 1, [1, 0, 0], -1.0)
    sethopping!(graphene, 2, 1, [0, 1, 0], -1.0)

    @test graphene.norbits == 2

    @test graphene.rlat ≈ [
        2π 0.0 0.0;
        -2π*tan(π/6) 2π/cos(π/6) 0.0;
        0.0 0.0 2π;
    ]

    @test caleigvals(graphene, [0.0, 0.0, 0.0]) ≈ [-3.0, 3.0]
    @test calhamiltonian(graphene, [0.0, 0.0, 0.0])*caleig(graphene, [0.0, 0.0, 0.0])[2][:, 1] ≈
        -3.0*caleig(graphene, [0.0, 0.0, 0.0])[2][:, 1]
    kdist, egvals = calband(graphene, [1 0; 0 1; 0 0], 3)
    @test kdist ≈ [0, 2π, 4π]
    @test egvals ≈ [-3.0 -1.0 -3.0; 3.0 1.0 3.0]
end
