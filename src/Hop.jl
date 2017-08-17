module Hop

export TightBindingModel, sethopping!, calhamiltonian, caleig

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
    hoppings::Dict{Vector{Int64}, Complex128}
end


"""
A tight binding model.

Generally, every tight binding model is three dimensional.
Thus every lattice vector should have three components.
Lower dimensional models should be simulated by vacuum layer.

# Arguments
- `lat::Matrix{Float64}`: Lattice vector.
  Lattice vectors should be provided as lat[:, i]
- `positions::Matrix{Float64}`: Atom postions.
  Atom positions should be provided as positions[:, i]
"""
function TightBindingModel(lat::Matrix{Float64}, positions::Matrix{Float64})
    @assert size(lat) == (3, 3) "Shape of lat is not correct."
    @assert size(positions, 1) == 3 "Shape of positions is not correct."
    a1 = lat[:, 1]
    a2 = lat[:, 2]
    a3 = lat[:, 3]
    rlat = zeros((3, 3))
    rlat[:, 1] = 2π*((a2×a3)/(a1⋅(a2×a3)))
    rlat[:, 2] = 2π*((a3×a1)/(a2⋅(a3×a1)))
    rlat[:, 3] = 2π*((a1×a2)/(a3⋅(a1×a2)))
    TightBindingModel(size(positions, 2), lat, rlat, positions, Dict())
end

"""
Set hoppings for a TightBindingModel t. Hoppings are expressed as ⟨0n|H|Rm⟩,
where R is a 3-element Vector{Int} representing lattice vector.
"""
function sethopping!(t::TightBindingModel, n::Int, m::Int, R::Vector{Int}, hopping::Union{Float64,Complex128})
    @assert n in 1:t.norbits && m in 1:t.norbits "No such orbit."
    t.hoppings[[[n, m]; R]] = hopping
    t.hoppings[[[m, n]; -R]] = hopping
    return
end

"""
Calculate Hamiltonian of a TightBindingModel t for a specific k point.
k is a 3-element Vector{Float} representing k point in relative coordinate.
"""
function calhamiltonian(t::TightBindingModel, k::Vector{Float64})
    @assert size(k) == (3,) "k point is not in correct shape"
    h = zeros((t.norbits, t.norbits))
    for (label, hopping) in t.hoppings
        h[label[1], label[2]] += exp(2π*im*(k⋅label[3:5]))*hopping
    end
    return h
end

"""
Calculate Hamiltonian of a TightBindingModel t for a specific k point.
k is a 3-element Vector{Float} representing k point in relative coordinate.
"""
function caleig(t::TightBindingModel, k::Vector{Float64}, calegvecs::Bool=false)
    hamiltonian = calhamiltonian(t, k)
    if calegvecs
        (egvals, egvecs) = eig(hamiltonian)
        egvals = real(egvals)
        perm = sortperm(egvals)
        sortedegvecs = zeros(size(egvecs))
        for i in 1:t.norbits
            sortedegvecs[:, i] = egvecs[:, perm[i]]
        end
        return (egvals[perm], sortedegvecs)
    else
        return sort(real(eigvals(hamiltonian)))
    end
end

end
