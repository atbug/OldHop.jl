using Hop, Base.Test


function graphene_test()
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
    @test calhamiltonian(graphene, [0.0, 0.0, 0.0])*caleig(graphene, [0.0, 0.0, 0.0], true)[2][:, 1] ≈
        -3.0*caleig(graphene, [0.0, 0.0, 0.0], true)[2][:, 1]
    kdist, egvals = calband(graphene, [1 0; 0 1; 0 0], 3)
    @test isapprox(kdist, [0.0, 1.0, 2.0], atol=1.0e-5)
    @test isapprox(egvals, [-3.0 -1.0 -3.0; 3.0 1.0 3.0], atol=1.0e-5)

    graphenesc = makesupercell(graphene, [2 0 0; 0 2 0; 0 0 1])

    @test graphenesc.norbits == 8
    @test caleig(graphenesc, [0.0, 0.0, 0.0]) ≈ sort(
        [
            caleig(graphene, [0.0, 0.0, 0.0]);
            caleig(graphene, [0.5, 0.0, 0.0]);
            caleig(graphene, [0.0, 0.5, 0.0]);
            caleig(graphene, [0.0, 0.5, 0.0])
        ]
    )

    graphenesc = makesupercell(graphene, [2 -1 0; 0 1 0; 0 0 1])

    @test graphenesc.norbits == 4
    @test caleig(graphenesc, [0.0, 0.0, 0.0]) ≈ sort(
        [
            caleig(graphene, [0.0, 0.0, 0.0]);
            caleig(graphene, [0.5, 0.0, 0.0])
        ]
    )

    graphenerb = cutedge(graphene, 1)
    @test caleig(graphenerb, [0.0, 0.0, 0.5]) == caleig(graphenerb, [0.0, 0.0, 0.0])
end


function line_test()
    line = TightBindingModel(
        [1.0 0 0; 0 1 0; 0 0 1],
        [0 1/3 2/3; 0 0 0; 0 0 0]
    )
    sethopping!(line, 1, 2, [0, 0, 0], 1)
    sethopping!(line, 2, 3, [0, 0, 0], 1)
    sethopping!(line, 3, 1, [1, 0, 0], 1)
    segment = cutedge(line, 1, true)
    @test caleig(segment, [0.0, 0, 0]) ≈ [-1, -1, 2]
end


function cluster_test()
    lat = [2 0 0.0; 0 2 0; 0 0 1]
    positions = [0.5 0.5 0.0; 0 0.5 0.0; 0 0 0.0]
    cluster = TightBindingModel(lat, positions)
    sethopping!(cluster, 1, 2, [0, 0, 0], 1.0)
    sethopping!(cluster, 1, 3, [0, 0, 0], 1.0)
    sethopping!(cluster, 2, 3, [0, 0, 0], 1.0)
    clusterm = addmagneticfield(cluster, 0.5)
    @test all(abs.(caleig(clusterm, [0.0, 0.0, 0.0]) - [-1.73205, 0, 1.73205]) .< 1.0e-6)
end


function test_spin()
    # graphene
    lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
    positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]
    graphene = TightBindingModel(lat, positions, spinful=true)
    hopping = [-1.0 0; 0 -1.0]
    sethopping!(graphene, 1, 2, [0, 0, 0], hopping)
    sethopping!(graphene, 2, 1, [1, 0, 0], hopping)
    sethopping!(graphene, 2, 1, [0, 1, 0], hopping)
    @test caleig(graphene, [0.0, 0.0, 0.0]) ≈ [-3.0, -3.0, 3.0, 3.0]

    # graphene
    lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
    positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]
    graphene = TightBindingModel(lat, positions, spinful=true)
    hopping = [-1.0 0; 0 -1.0]
    sethopping!(graphene, 1, 2, [0, 0, 0], hopping/2)
    sethopping!(graphene, 1, 2, [0, 0, 0], hopping/2)
    sethopping!(graphene, 2, 1, [1, 0, 0], hopping)
    sethopping!(graphene, 2, 1, [0, 1, 0], hopping)
    @test caleig(graphene, [0.0, 0.0, 0.0]) ≈ [-3.0, -3.0, 3.0, 3.0]

    # Kane-Mele
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
    @test isapprox(caleig(tm, [2/3, 1/3, 0.0]),
        [-2.55885, -0.558846, 0.558846, 2.55885], atol=1.0e-5)

    lat = eye(3)
    positions = zeros(3, 1)
    atom = TightBindingModel(lat, positions, spinful=true)
    sethopping!(atom, 1, 1, [0, 0, 0], [0 -im; im 0])
    @test caleig(atom, [0.0, 0.0, 0.0]) ≈ [-1.0, 1.0]
end

graphene_test()
line_test()
cluster_test()
test_spin()
