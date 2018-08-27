module Hop
using StaticArrays, LinearAlgebra

export TightBindingModel, sethopping!,
    gethamiltonian, geteig, geteigvals,
    KPath, getband


struct TightBindingModel
    "number of orbits"
    norbits::Int
    "lattice vector"
    lat::Matrix{Float64}
    "reciprocal lattice vector"
    rlat::Matrix{Float64}
    "atom positions"
    positions::Matrix{Float64}
    "hoppings"
    hoppings::Dict{SVector{3,Int64},Matrix{ComplexF64}}
end


"""
```julia
TightBindingModel(lat::Matrix{Float64}, positions::Matrix{Float64}; spinful::Bool=false)
```
A tight binding model.

Generally, every tight binding model is three dimensional.
Therefore, every lattice vector should have three components.
Lower dimensional models should be simulated by vacuum layer.

# Constructor Arguments
- `lat`: lattice vector.
  Lattice vectors should be provided in columns.
- `positions`: atom positions in reduced coordinate.
  Atom positions should be provided in columns.
- `spinful`: false for spinless systems and true for spinful systems. If `spinful`
  is true, `norbits` will be twice the number of `size(positions, 1)`. Orbits are
  ordered as (|1↑⟩, |1↓⟩, |2↑⟩, |2↓⟩, ...).

# Fields
- `norbits::Int`: number of orbits.
- `lat::Matrix{Float64}`: lattice vectors stored in columns.
- `rlat::Matrix{Float64}`: reciprocal lattice vectors stored in columns.
- `positions::Matrix{Float64}`: position of orbits in reduced coordinate stored in columns.
- `hoppings::Dict{SVector{3,Int64},Matrix{ComplexF64}}`: hoppings stored as
   R->⟨0m|H|Rn⟩.
"""
function TightBindingModel(lat::Matrix{Float64}, positions::Matrix{Float64}; spinful::Bool=false)
    @assert size(lat) == (3, 3) "Size of lat is not correct."
    @assert size(positions, 1) == 3 "Size of positions is not correct."
    rlat = 2*π*inv(lat)'
    if !spinful
        return TightBindingModel(size(positions, 2), lat, rlat, positions, Dict())
    else
        spinpositions = zeros(3, size(positions, 2)*2)
        for i in 1:size(positions, 2)
            spinpositions[:, 2*i-1] = positions[:, i]
            spinpositions[:, 2*i] = positions[:, i]
        end
        return TightBindingModel(size(positions, 2)*2, lat, rlat, spinpositions, Dict())
    end
end


function Base.show(io::IO, t::TightBindingModel)
    print(io, "TightBindingModel: $(t.norbits) orbitals.")
end


"""
```julia
sethopping!(t::TightBindingModel, m::Int64, n::Int64, R, hopping;
    mode::Symbol=:a)
```

Set ⟨0m|H|Rn⟩ to `hopping`. `hopping::Number` for spinless models and
`hopping::Matrix{<:Number}` for spinful models. For spinful models,
`size(hopping)` should be (2, 2) and the basis for `hopping` is (|↑⟩, |↓⟩).
`mode` has two possible values: :a for add mode and :s for set or reset mode.
"""
function sethopping!(t::TightBindingModel, m::Int64, n::Int64, R,
    hopping::Number; mode::Symbol=:a)
    @assert (m in 1:t.norbits) && (n in 1:t.norbits) "No such orbit."
    @assert length(R) == 3 "Size of R is not correct."
    if !(R in keys(t.hoppings))
        @assert !(-R in keys(t.hoppings))
        t.hoppings[R] = zeros(ComplexF64, (t.norbits, t.norbits))
        t.hoppings[-R] = zeros(ComplexF64, (t.norbits, t.norbits))
    end
    if mode == :s
        t.hoppings[R][m, n] = hopping
        t.hoppings[-R][n, m] = conj(hopping)
    elseif mode == :a
        t.hoppings[R][m, n] += hopping
        if (m, n, R) != (n, m, -R) # not onsite energy
            t.hoppings[-R][n, m] += conj(hopping)
        end
    else
        error("Unknown mode.")
    end
    return nothing
end


# For spinful model
function sethopping!(t::TightBindingModel, m::Int64, n::Int64, R,
    hopping::Matrix{<:Number}; mode::Symbol=:a)
    @assert iseven(t.norbits) "Not a spinful model."
    @assert (m in 1:(t.norbits÷2)) && (n in 1:(t.norbits÷2)) "No such orbit."
    @assert size(hopping) == (2, 2) "Size of hopping is not correct."
    sethopping!(t, 2m-1, 2n-1, R, hopping[1, 1], mode=mode)
    sethopping!(t, 2m, 2n-1, R, hopping[2, 1], mode=mode)
    if (m, n, R) != (n, m, -R) # not onsite energy
        sethopping!(t, 2m-1, 2n, R, hopping[1, 2], mode=mode)
    end
    sethopping!(t, 2m, 2n, R, hopping[2, 2], mode=mode)
    return nothing
end


"""
```julia
gethamiltonian(t::TightBindingModel, k::Vector{<:Real}) --> Matrix{ComplexF64}
```

Calculate Hamiltonian of a TightBindingModel t for a specific k point. k should
be provided in reduced coordinate.
"""
function gethamiltonian(t::TightBindingModel, k::Vector{<:Real})
    @assert size(k) == (3,) "Size of k is not correct."
    h = zeros(ComplexF64, (t.norbits, t.norbits))
    for (R, hopping) = t.hoppings
        h += exp(2π*im*(k⋅R))*hopping
    end
    @assert ishermitian(h)
    return h
end


"""
```julia
geteig(t::TightBindingModel, k::Vector{<:Real})
```

Calculate eigenvalues and eigenvectors of t. k should be provided in reduced coordinate.

# Return
`(egvals::Vector{Float64}, egvecs::Matrix{ComplexF64})`,
Eigenvectors are stored in columns and eigenvalues are sorted from small to large.
"""
function geteig(t::TightBindingModel, k::Vector{<:Real})
    @assert size(k) == (3,) "Size of k is not correct."
    hamiltonian = gethamiltonian(t, k)
    egvals, egvecs = eigen(Hermitian(hamiltonian))
    perm = sortperm(egvals)
    return (egvals[perm], egvecs[:, perm])
end


"""
```julia
geteigvals(t::TightBindingModel, k::Vector{<:Real})
  --> Vector{Float64}
```

Calculate eigenvalues and eigenvectors of t. k should be provided in reduced coordinate.

Eigenvalues are sorted from small to large.
"""
function geteigvals(t::TightBindingModel, k::Vector{<:Real})
    @assert size(k) == (3,) "Size of k is not correct."
    hamiltonian = gethamiltonian(t, k)
    return sort(eigvals(Hermitian(hamiltonian)))
end


struct KPath
    rlat::Matrix{Float64}
    nodes::Matrix{Float64}
    distances::Vector{Float64}
    kpts::Matrix{Float64}
    ndiv::Int64
    nkpts::Int64
end


"""
```julia
KPath(rlat::Matrix{<:Real}, nodes::Matrix{<:Real}, ndiv::Int64)
```

A path in reciprocal space.

# Constructor Arguments
 - `rlat`: reciprocal lattice vector stored in columns.
 - `nodes`: start and end points of the path stored in columns in reduced coordinates.
 - `ndiv`: number of divisions on each section of the path.

# Fields
 - `rlat::Matrix{Float64}`: see explanations above.
 - `nodes::Matrix{Float64}`: see explanations above.
 - `distances::Vector{Float64}`: distance along the path for each k point.
 - `kpts::Matrix{Float64}`: k points stored in columns  in reduced coordinates.
 - `ndiv::Int64`: see explanations above.
 - `nkpts::Int64`: number of k points on the path.
"""
function KPath(rlat::Matrix{<:Real}, nodes::Matrix{<:Real}, ndiv::Int64)
    @assert iseven(size(nodes, 2))
    npaths = size(nodes, 2)÷2
    nkpts = ndiv*npaths
    distances = zeros(nkpts)
    kpts = zeros(3, nkpts)
    for ipath = 1:npaths
        dk = (nodes[:, 2*ipath]-nodes[:, 2*ipath-1])/(ndiv-1)
        dkn = norm(rlat*dk) # real norm of dk
        if ipath == 1
            distance0 = 0
        else
            distance0 = distances[(ipath-1)*ndiv]
        end
        for ikpt = 1:ndiv
            distances[(ipath-1)*ndiv+ikpt] = dkn*(ikpt-1) + distance0
            kpts[:, (ipath-1)*ndiv+ikpt] = nodes[:, 2*ipath-1]+dk*(ikpt-1)
        end
    end
    return KPath(rlat, nodes, distances, kpts, ndiv, nkpts)
end


function Base.show(io::IO, kp::KPath)
    print(io, "KPath: $(kp.nkpts÷kp.ndiv) line(s).")
end


"""
```julia
getband(t::TightBindingModel, kp::KPath)
    --> Matrix{Float64}
```

Calculate bands of `t` on `kp`.

Return energies at each k point stored in columns.
"""
function getband(t::TightBindingModel, kp::KPath)
    bands = zeros(t.norbits, kp.nkpts)
    for ik=1:kp.nkpts
        bands[:, ik] = geteigvals(t, kp.kpts[:, ik])
    end
    return bands
end

include("floquet.jl")

end
