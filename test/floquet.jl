let
    f = Hop.Floquet.FloquetHamiltonian([ones(2, 2), 2*ones(2, 2)*im], 1, 2)
    @test f.Hf[1:2, 1:2] ≈ ones(2, 2)+2*I
    @test f.Hf[5:6, 3:4] ≈ 2*ones(2, 2)*im
    @test f.Hf[3:4, 5:6] ≈ -2*ones(2, 2)*im
end

let
    graphene = getgraphene()
    A = [1.0, 0.0, 0.0]
    Ω = 1.0
    Γ = 1.0
    μ = 0.0
    f = Hop.Floquet.get_illuminated_hamiltonian(graphene, [0.2, 0.0, 0.0], A=A, Ω=Ω, harmonics_cutoff=5)
    @test Hop.Floquet.get_floquet_occupation(f; μ=μ, Γ=Γ) ≈ [0.23, 0.77] atol=1.0e-2
end

let
    graphene = getgraphene()
    kp = KPath(graphene.rlat, [1 0; 0 1; 0 0], 10)
    bands = Hop.Floquet.get_illuminated_band(graphene, A=[1.0, 0.0, 0.0],
        Ω=1.0, kp=kp, harmonics_cutoff=5)
    @test bands[21, 4] ≈ 5.589253 atol=1.0e-6
end
