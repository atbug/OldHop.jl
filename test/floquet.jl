let
    f = Hop.Floquet.FloquetHamiltonian([ones(2, 2), 2*ones(2, 2)*im], 1, 2)
    @test f.Hf[1:2, 1:2] ≈ ones(2, 2)+2*I
    @test f.Hf[5:6, 3:4] ≈ 2*ones(2, 2)*im
    @test f.Hf[3:4, 5:6] ≈ -2*ones(2, 2)*im
end

let
    lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
    positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]
    graphene = TightBindingModel(lat, positions)
    sethopping!(graphene, 1, 2, [0, 0, 0], -1.0) # ⟨1|H|(0, 0, 0)2⟩=-1
    sethopping!(graphene, 2, 1, [1, 0, 0], -1.0) # ⟨2|H|(1, 0, 0)1⟩=-1
    sethopping!(graphene, 2, 1, [0, 1, 0], -1.0) # ⟨2|H|(0, 1, 0)1⟩=-1

    A = [1.0, 0.0, 0.0]
    Ω = 1.0
    Γ = 1.0
    μ = 0.0
    f = Hop.Floquet.get_illuminated_hamiltonian(graphene, [0.2, 0.0, 0.0], A=A, Ω=Ω, harmonics_cutoff=5)
    @test Hop.Floquet.get_floquet_occupation(f; μ=μ, Γ=Γ) ≈ [0.23, 0.77] atol=1.0e-2
end

let
    lat = [1.0 0.5 0.0; 0.0 (√3)/2 0.0; 0.0 0.0 1.0]
    positions = [1/3 2/3; 1/3 2/3; 0.0 0.0]
    graphene = TightBindingModel(lat, positions)
    sethopping!(graphene, 1, 2, [0, 0, 0], -1.0) # ⟨1|H|(0, 0, 0)2⟩=-1
    sethopping!(graphene, 2, 1, [1, 0, 0], -1.0) # ⟨2|H|(1, 0, 0)1⟩=-1
    sethopping!(graphene, 2, 1, [0, 1, 0], -1.0) # ⟨2|H|(0, 1, 0)1⟩=-1

    kdist, egvals = Hop.Floquet.get_illuminated_band(graphene, A=[1.0, 0.0, 0.0],
        Ω=1.0, kpath=[1 0; 0 1; 0 0], harmonics_cutoff=5, ndiv=100)

    @test kdist[42] ≈ 5.204254 atol=1.0e-6
    @test egvals[21, 42] ≈ 5.626495 atol=1.0e-6
end
