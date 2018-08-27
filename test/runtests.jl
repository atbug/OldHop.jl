using Hop, Test, LinearAlgebra

function getgraphene()
    lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
    positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]
    graphene = TightBindingModel(lat, positions)
    sethopping!(graphene, 1, 2, [0, 0, 0], -1.0)
    sethopping!(graphene, 2, 1, [1, 0, 0], -1.0)
    sethopping!(graphene, 2, 1, [0, 1, 0], -1.0)
    return graphene
end

let
    graphene = getgraphene()
    @test graphene.norbits == 2
    @test graphene.rlat ≈ [
        2π 0.0 0.0;
        -2π*tan(π/6) 2π/cos(π/6) 0.0;
        0.0 0.0 2π;
    ]
    @test geteigvals(graphene, [0.0, 0.0, 0.0]) ≈ [-3.0, 3.0]
    @test gethamiltonian(graphene, [0.0, 0.0, 0.0])*geteig(graphene, [0.0, 0.0, 0.0])[2][:, 1] ≈
        -3.0*geteig(graphene, [0.0, 0.0, 0.0])[2][:, 1]
    kp = KPath(graphene.rlat, [1 0 1 0; 0 1 0 1; 0 0 0 0], 3)
    bands = getband(graphene, kp)
    @test kp.distances ≈ [0, 2π, 4π, 4π, 6π, 8π]
    @test bands[:, 1:3] ≈ [-3.0 -1.0 -3.0; 3.0 1.0 3.0]
end

let    # Kane-Mele
    lat = [1 0.5 0; 0 (√3)/2 0; 0 0 1]
    positions = [1/3 2/3; 1/3 2/3; 0 0]
    tm = TightBindingModel(lat, positions, spinful=true)
    σ0 = [1 0; 0 1]
    σ1 = [0 1; 1 0]
    σ2 = [0 -im; im 0]
    σ3 = [1 0; 0 -1]
    onsite = 1.0
    t = 1.0
    so = 0.6*t*0.5
    sethopping!(tm, 1, 1, [0, 0, 0], σ0*onsite)
    sethopping!(tm, 2, 2, [0, 0, 0], -σ0*onsite)
    sethopping!(tm, 1, 2, [0, 0, 0], σ0*t)
    sethopping!(tm, 1, 2, [0, -1, 0], σ0*t)
    sethopping!(tm, 1, 2, [-1, 0, 0], σ0*t)
    sethopping!(tm, 1, 1, [0, 1, 0], -im*so*σ3)
    sethopping!(tm, 1, 1, [1, 0, 0], im*so*σ3)
    sethopping!(tm, 1, 1, [1, -1, 0], -im*so*σ3)
    sethopping!(tm, 2, 2, [0, 1, 0], im*so*σ3)
    sethopping!(tm, 2, 2, [1, 0, 0], -im*so*σ3)
    sethopping!(tm, 2, 2, [1, -1, 0], im*so*σ3)
    @test geteigvals(tm, [2/3, 1/3, 0.0]) ≈ [-2.55885, -0.558846, 0.558846, 2.55885] atol=1.0e-5
end

include("floquet.jl")
