module Hop

export TightBindingModel, sethopping!, addmagneticfield!, makesupercell, makecluster, calhamiltonian, caleig


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
- `lat::Matrix{Float64}`: lattice vector.
  Lattice vectors should be provided as lat[:, i]
- `positions::Matrix{Float64}`: atom postions.
  Atom positions should be provided as positions[:, i]
"""
function TightBindingModel(lat::Matrix{Float64}, positions::Matrix{Float64})
    @assert size(lat) == (3, 3) "Size of lat is not correct."
    @assert size(positions, 1) == 3 "Size of positions is not correct."
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
Add constant magnetic field in z direction for a TightBindingModel.

# Arguments
- `t::TightBindingModel`: a TightBindingModel. t has to be a cluster.
- `B::Float64`: magnetic field in z direction. B is actually Be/h, thus
  its unit is 1/[length]^2. e here is fundamental charge.
  Since electron charge is -e, positive B means -z direction for electron system.
"""
function addmagneticfield!(t::TightBindingModel, B::Float64)
    for (label, hopping) in t.hoppings
        # landau gauge
        absolute_position_n = t.lat*t.positions[:,label[1]]
        absolute_position_m = t.lat*t.positions[:,label[2]]
        t.hoppings[label] = hopping*exp(
            im*2π*B*(absolute_position_n[2]-absolute_position_m[2])*(absolute_position_n[1]+absolute_position_m[1])/2
        )
    end
    return
end


"""
Create a supercell TightBindingModel from original TightBindingModel.

# Arguments
- `t::TightBindingModel`: the original TightBindingModel.
- `ncells::Vector{Int}`: a 3-element vector representing number of unit cells
  in three spatial directions.
"""
function makesupercell(t::TightBindingModel, ncells::Vector{Int})
    @assert size(ncells) == (3,) "Size of ncells is not correct."
    lat = [t.lat[:, 1]*ncells[1] t.lat[:, 2]*ncells[2] t.lat[:, 3]*ncells[3]]
    positions = Vector{Float64}()
    for k in 1:ncells[3]
        for j in 1:ncells[2]
            for i in 1:ncells[1]
                for n in 1:t.norbits
                    positions = [positions; t.positions[:, n] + [i-1, j-1, k-1]]
                end
            end
        end
    end
    positions = reshape(positions, (3, prod(ncells)*t.norbits))
    positions[1, :] /= ncells[1]
    positions[2, :] /= ncells[2]
    positions[3, :] /= ncells[3]
    sc = TightBindingModel(lat, positions)

    function getsupercellpositionindex(i, j, k, n)
        return ((i-1) + (j-1)*ncells[1] + (k-1)*ncells[1]*ncells[2])*t.norbits + n
    end

    for k in 1:ncells[3]
        for j in 1:ncells[2]
            for i in 1:ncells[1]
                for (label, hopping) in t.hoppings
                    sethopping!(
                        sc,
                        getsupercellpositionindex(i, j, k, label[1]),
                        getsupercellpositionindex(mod((i+label[3]-1),ncells[1])+1, mod((j+label[4]-1),ncells[2])+1, mod((k+label[5]-1),ncells[3])+1, label[2]),
                        [fld(i+label[3]-1, ncells[1]), fld(j+label[4]-1, ncells[2]), fld(k+label[5]-1, ncells[3])],
                        hopping
                    )
                end
            end
        end
    end
    return sc
end


"""
Create cluster by cutting off all hoppings between cells.

# Arguments
- `t::TightBindingModel`: a TightBindingModel.
"""
function makecluster(t::TightBindingModel)
    c = deepcopy(t)
    for (label, hopping) in c.hoppings
        if label[3] != 0 || label[4] != 0 || label[5] != 0
            pop!(c.hoppings, label)
        end
    end
    return c
end


"""
Calculate Hamiltonian of a TightBindingModel t for a specific k point.
k is a 3-element Vector{Float} representing k point in relative coordinate.
"""
function calhamiltonian(t::TightBindingModel, k::Vector{Float64})
    @assert size(k) == (3,) "Size of k is not correct."
    h = zeros(Complex128, (t.norbits, t.norbits))
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
    @assert size(k) == (3,) "Size of k is not correct."
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
