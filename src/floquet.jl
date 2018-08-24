module Floquet
using ..Hop, LinearAlgebra

export cal_floquet_hamiltonian, cal_illuminated_hamiltonian, cal_illuminated_band

"""
```julia
cal_floquet_hamiltonian(Hs::Vector{Matrix{Number}}, Ω::Real, N::Integer)
    --> Matrix{ComplexF64}
```

Calculate Floquet Hamiltonian. Generally, Floquet theory deals with time periodic
Hamiltonians ``H(t+Ω)=H(t)``.

`Hs` is a list of Hamiltonian looking like ``[H_0, H_1, H_2, ...]``.
``H_i`` is coefficient of Fourier expansion of ``H(t)``:
``H(t)=H_0+(H_1 e^{iΩt}+H_2 e^{2iΩt}+...+c.c.)``.

The basis of returned Floquet Hamiltonian is
``(e^{-iNΩt}, e^{-i(N-1)Ωt}, ..., e^{iNΩt})``
where internal degree of freedom is omitted. Therefore, ``H^F_{mn}=H_{n-m}-mΩδ_{mn}``
"""
function cal_floquet_hamiltonian(Hs::Vector{Matrix{T}}, Ω::Real, N::Integer) where T<:Number
    dim = size(Hs[1], 1)
    dimF = dim*(2N+1)
    nhams = length(Hs)
    HF = zeros(ComplexF64, (dimF, dimF))
    for j = 0:(nhams-1)
        for i = 1:(2N+1-j)
            if j == 0
                HF[dim*(i-1)+1:dim*i, dim*(i-1)+1:dim*i] = Hs[1] + (i-N-1)*Ω*I
            else
                HF[dim*(i-1)+1:dim*i, dim*(i-1+j)+1:dim*(i+j)] = Hs[j+1]
                HF[dim*(i-1+j)+1:dim*(i+j), dim*(i-1)+1:dim*i] = Hs[j+1]'
            end
        end
    end
    return Hermitian(HF)
end


"""
```julia
cal_illuminated_hamiltonian(t::TightBindingModel, k::Vector{<:Real}; A::Vector{<:Number},
    Ω::Real, N::Integer=2) --> Matrix{ComplexF64}
```

Calculate Floquet Hamiltonian of `t` at `k` point under the illumination of light with frequency Ω.
Peierls substitution is performed to the lowest order. Only electric field is taken
into account. The parameter `A` looks like ``[A_x, A_y, A_z]`` denoting vector potential
``A(t)=[A_x, A_y, A_z]e^{iΩt}+c.c.``. Floquet Hamiltonian is truncated
up to `N` harmonics.
"""
function  cal_illuminated_hamiltonian(t::TightBindingModel, k::Vector{<:Real}; A::Vector{<:Number},
    Ω::Real, N::Integer=2)
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
    return cal_floquet_hamiltonian([convert(Array, H0), H1], Ω, N)
end


"""
```julia
cal_illuminated_band(t::TightBindingModel; A::Vector{<:Number}, Ω::Real,
    kpath::Matrix{<:Real}, N::Integer=2, ndiv::Int64=100) --> (Vector{Float64}, Matrix{Float64})
```

Calculate bands of `t` under the illumination of light with frequency Ω.
`kpath` is a (3, x) size matrix where x is an even number and
should be provided in reduced coordinates. Peierls substitution is performed to
the lowest order. Only electric field is taken
into account. The parameter `A` looks like ``[A_x, A_y, A_z]`` denoting vector potential
``A(t)=[A_x, A_y, A_z]e^{iΩt}+c.c.``. Floquet Hamiltonian is truncated
up to `N` harmonics.

This function returns (`kdist`, `egvals`). `kdist` is the distance of k points and
`egvals` is the Floquet energies of band stored in column at each k.
"""
function cal_illuminated_band(t::TightBindingModel; A::Vector{<:Number}, Ω::Real,
    kpath::Matrix{<:Real}, N::Integer=2, ndiv::Int64=100)
    @assert iseven(size(kpath, 2))
    npaths = size(kpath, 2)÷2
    nkpts = ndiv*npaths
    kdist = zeros(nkpts)
    egvals = zeros(t.norbits*(2N+1), nkpts)
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
            egvals[:, (ipath-1)*ndiv+ikpt] = eigvals(cal_illuminated_hamiltonian(t, k, A=A, Ω=Ω, N=N))
        end
    end
    return (kdist, egvals)
end

end
