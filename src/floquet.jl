module Floquet
using ..Hop, LinearAlgebra

export calfloquethamiltonian

"""
```julia
calfloquethamiltonian(Hs::Vector{Matrix{Number}}, Ω::Real, N::Integer)
    --> Matrix{ComplexF64}
```

Calculate Floquet Hamiltonian. Generally, Floquet theory deals with time periodic
Hamiltonians H(t+Ω)=H(t).

`Hs` is a list of Hamiltonian looking like [H0, H1, H2, ...]. Hi is coefficient
of Fourier expansion of H(t): H(t)=H0+H1*exp(iΩt)+H2*exp(iΩt)+...

The basis of returned Floquet Hamiltonian is (e^{-iNΩt}, e^{-i(N-1)Ωt}, ..., e^{iNΩt})
where internal degree of freedom is omitted. Therefore, HF_mn=H_{n-m}-mΩδ_{mn}
"""
function calfloquethamiltonian(Hs::Vector{Matrix{T}}, Ω::Real, N::Integer) where T<:Number
    dim = size(Hs[1], 1)
    dimF = dim*(2N+1)
    nhams = length(Hs)
    HF = zeros(ComplexF64, (dimF, dimF))
    for j=0:(nhams-1)
        for i=1:(2N+1-j)
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

end
