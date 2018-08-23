module Floquet
using ..Hop, LinearAlgebra

export cal_floquet_hamiltonian, cal_illuminated_hamiltonian

"""
```julia
cal_floquet_hamiltonian(Hs::Vector{Matrix{Number}}, Ω::Real, N::Integer)
    --> Matrix{ComplexF64}
```

Calculate Floquet Hamiltonian. Generally, Floquet theory deals with time periodic
Hamiltonians ``H(t+Ω)=H(t)``.

`Hs` is a list of Hamiltonian looking like ``[H_0, H_1, H_2, ...]``.
``H_i`` is coefficient of Fourier expansion of ``H(t)``:
``H(t)=H_0+(H_1 e^(iΩt)+H_2 e^(2iΩt)+...+c.c.)``.

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
cal_illuminated_hamiltonian(t::TightBindingModel, A::Vector{<:Number}, Ω::Real, k::Vector{<:Real}, N::Integer)
    --> Matrix{ComplexF64}
```

Calculate Floquet Hamiltonian of `t` at `k` point under the illumination of light with frequency Ω.
Peierls substitution is performed to the lowest order. Only electric field is taken
into account. The parameter `A` looks like ``[A_x, A_y, A_z]`` denoting vector potential
``A(t)=[A_x, A_y, A_z]e^(iΩt)+c.c.``. Floquet Hamiltonian is truncated
up to `N` harmonics.
"""
function  cal_illuminated_hamiltonian(t::TightBindingModel, A::Vector{<:Number}, Ω::Real, k::Vector{<:Real}, N::Integer)
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
        println(R)
        println(Peierlshopping)
    end
    println(H1)
    return calfloquethamiltonian([convert(Array, H0), H1], Ω, N)
end

end
