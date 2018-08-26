module Floquet
using ..Hop, LinearAlgebra, QuadGK

export FloquetHamiltonian, get_floquet_occupation, get_illuminated_hamiltonian, get_illuminated_band


struct FloquetHamiltonian
    Hs::Vector{Matrix{ComplexF64}}
    Hf::Matrix{ComplexF64}
    Ω::Float64
    nstates::Int64
    nfstates::Int64
    harmonics_cutoff::Int64
end

function Base.show(io::IO, f::FloquetHamiltonian)
    print(io, "FloquetHamiltonian: $(f.nstates) states with $(2*f.harmonics_cutoff+1) harmonics.")
end


"""
FloquetHamiltonian(Hs::Vector{Matrix{T}}, Ω::Real, harmonics_cutoff::Integer) where T<:Number
```
Floquet Hamiltonian.

Generally, Floquet theory deals with time periodic
Hamiltonians ``H(t+Ω)=H(t)``.

# Constructor Arguments

 - `Hs` is a list of Hamiltonian looking like ``[H_0, H_1, H_2, ...]``.
   ``H_n`` is defined as: ``H_n=\\frac{1}{T}\\int_0^T e^{inΩt} H(t) dt``.
 - `Ω` is the frequency (provided in eV) of the Floquet system.
 - `harmonics_cutoff` determines how many harmonics are used in Floquet Hamiltonian

# Fields
 - `Hs::Vector{Matrix{ComplexF64}}`, see above explanations.
 - `Hf::Matrix{ComplexF64}` is Floquet Hamiltonian.
   With basis functions defined as ``u_n=e^{-inΩt}``, ``H^F_{mn}=H_{m-n}-mΩδ_{mn}``.
   The basis of returned Floquet Hamiltonian is
   ``(e^{iNΩt}, e^{i(N-1)Ωt}, ..., e^{-iNΩt})``
   where internal states are omitted.
 - `Ω::Float64`, see above explanations.
 - `nstates::Int64` is number of internal states.
 - `nfstates::Int64` is number of Floquet states.
 - `harmonics_cutoff::Int64`, see above explanations.
"""
function FloquetHamiltonian(Hs::Vector{Matrix{T}}, Ω::Real, harmonics_cutoff::Integer) where T<:Number
    nstates = size(Hs[1], 1)
    nfstates = nstates*(2*harmonics_cutoff+1)
    nhams = length(Hs)
    Hf = zeros(ComplexF64, (nfstates, nfstates))
    for j = 0:(nhams-1)
        for i = 1:(2harmonics_cutoff+1-j)
            if j == 0
                Hf[nstates*(i-1)+1:nstates*i, nstates*(i-1)+1:nstates*i] = Hs[1]-(i-harmonics_cutoff-1)*Ω*I
            else
                Hf[nstates*(i-1)+1:nstates*i, nstates*(i-1+j)+1:nstates*(i+j)] = Hs[j+1]'
                Hf[nstates*(i-1+j)+1:nstates*(i+j), nstates*(i-1)+1:nstates*i] = Hs[j+1]
            end
        end
    end
    @assert ishermitian(Hf)
    return FloquetHamiltonian(Hs, Hf, Ω, nstates, nfstates, harmonics_cutoff)
end


function get_floquet_GR(f::FloquetHamiltonian, ω::Real; Γ::Real=1.0)
    return inv((ω+im*Γ/2)*I-f.Hf)
end


function get_floquet_Σless(f::FloquetHamiltonian, ω::Real; Γ::Real=1.0, μ::Real=1.0)
    Σ = zeros(ComplexF64, size(f.Hf))
    for i=1:(2*f.harmonics_cutoff+1)
        Σ[(i-1)*f.nstates+1:i*f.nstates, (i-1)*f.nstates+1:i*f.nstates] =
            im*Γ*fermi(ω+(i-f.harmonics_cutoff-1)*f.Ω, μ)*Matrix(1.0I, f.nstates, f.nstates)
    end
    return Σ
end


function get_floquet_Gless(f::FloquetHamiltonian, ω::Real; Γ::Real=1.0, μ::Real=1.0)
    GR = get_floquet_GR(f, ω, Γ=Γ)
    Σless = get_floquet_Σless(f, ω, Γ=Γ, μ=μ)
    return GR*Σless*(GR')
end


"""
```julia
get_floquet_occupation(f::FloquetHamiltonian; Γ::Real=1.0, μ::Real=0.0,
    atol::Real=1.0e-3) --> Vector{Float64}
```

Calculate occupation number of `f` assuming coupling to a heat bath with coupling
constant `Γ`. Chemical potential of the heat bath is `μ`. `atol` is the integration
absolute tolerance.
"""
function get_floquet_occupation(f::FloquetHamiltonian; Γ::Real=1.0, μ::Real=0.0,
    atol::Real=1.0e-3)
    _, egvecs = eigen(f.Hf)
    occs = zeros(f.nstates)
    for i=(f.harmonics_cutoff*f.nstates+1):(f.harmonics_cutoff+1)*f.nstates
        occ = -im*quadgk(
        ω->egvecs[:, i]'*get_floquet_Gless(f, ω, Γ=Γ, μ=μ)*egvecs[:, i],
        -Inf, Inf, atol=atol)[1]/(2π)
        @assert imag(occ) < 1.0e-6
        occs[i-f.harmonics_cutoff*f.nstates] = real(occ)
    end
    return occs
end


"""
```julia
get_illuminated_hamiltonian(t::TightBindingModel, k::Vector{<:Real};
    A::Vector{<:Number}, Ω::Real, harmonics_cutoff::Integer=2)
    --> FloquetHamiltonian
```

Calculate Floquet Hamiltonian of `t` at `k` point under the illumination of light with frequency Ω.
Peierls substitution is performed to the lowest order. Only electric field is taken
into account. The parameter `A` looks like ``[A_x, A_y, A_z]`` denoting vector potential
``A(t)=[A_x, A_y, A_z]e^{-iΩt}+c.c.``. Floquet Hamiltonian is truncated
up to `harmonics_cutoff` harmonics.
"""
function  get_illuminated_hamiltonian(t::TightBindingModel, k::Vector{<:Real};
    A::Vector{<:Number}, Ω::Real, harmonics_cutoff::Integer=2)
    @assert length(A) == 3 "Length of A is not correct."

    H0 = calhamiltonian(t, k)

    # Unfortunately we cannot calculate Hamiltonian by calhamiltonian since
    # H1 is not necessarily Hermitian.
    cpositions = t.lat*t.positions # Cartesian positions
    H1 = zeros(ComplexF64, size(H0))
    for (R, hopping) = t.hoppings
        cR = t.lat*R
        Peierlshopping = zeros(ComplexF64, size(H0))
        for n = 1:t.norbits
            for m = 1:t.norbits
                Peierlshopping[m, n] = im*hopping[m, n]*((cpositions[:, n]-cpositions[:, m]+cR)⋅A)
            end
        end
        H1 += exp(2π*im*(k⋅R))*Peierlshopping
    end

    return FloquetHamiltonian([H0, H1], Ω, harmonics_cutoff)
end


"""
```julia
function get_illuminated_band(t::TightBindingModel; A::Vector{<:Number}, Ω::Real,
    kpath::Matrix{<:Real}, harmonics_cutoff::Integer=2, ndiv::Int64=100)
    --> (Vector{Float64}, Matrix{Float64})
```

Calculate bands of `t` under the illumination of light with frequency Ω.
`kpath` is a (3, x) size matrix where x is an even number and
should be provided in reduced coordinates. Peierls substitution is performed to
the lowest order. Only electric field is taken
into account. The parameter `A` looks like ``[A_x, A_y, A_z]`` denoting vector potential
``A(t)=[A_x, A_y, A_z]e^{-iΩt}+c.c.``. Floquet Hamiltonian is truncated
up to `harmonics_cutoff` harmonics.

This function returns (`kdist`, `egvals`). `kdist` is the distance of k points and
`egvals` is the Floquet energies of band stored in column at each k.
"""
function get_illuminated_band(t::TightBindingModel; A::Vector{<:Number}, Ω::Real,
    kpath::Matrix{<:Real}, harmonics_cutoff::Integer=2, ndiv::Int64=100)
    @assert iseven(size(kpath, 2))
    npaths = size(kpath, 2)÷2
    nkpts = ndiv*npaths
    kdist = zeros(nkpts)
    egvals = zeros(t.norbits*(2harmonics_cutoff+1), nkpts)
    for ipath = 1:npaths
        dk = (kpath[:, 2*ipath]-kpath[:, 2*ipath-1])/(ndiv-1)
        dkn = norm(t.rlat*dk) # real norm of dk
        if ipath == 1
            kdist0 = 0
        else
            kdist0 = kdist[(ipath-1)*ndiv]
        end
        for ikpt = 1:ndiv
            kdist[(ipath-1)*ndiv+ikpt] = dkn*(ikpt-1) + kdist0
            k = kpath[:, 2*ipath-1]+dk*(ikpt-1)
            egvals[:, (ipath-1)*ndiv+ikpt] = eigvals(
                get_illuminated_hamiltonian(t, k, A=A, Ω=Ω, harmonics_cutoff=harmonics_cutoff).Hf
            )
        end
    end
    return (kdist, egvals)
end

fermi(E, μ; T=0) = E<μ ? 1.0 : 0.0

end
