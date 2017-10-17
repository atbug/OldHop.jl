module Hop
using StaticArrays

export TightBindingModel, sethopping!, calhamiltonian, caleig, makesupercell, cutedge, addmagneticfield


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
    hoppings::Dict{Tuple{Int64,Int64,SVector{3,Int64}},Complex128}
end


"""
A tight binding model.

Generally, every tight binding model is three dimensional.
Thus every lattice vector should have three components.
Lower dimensional models should be simulated by vacuum layer.

# Constructor Arguments
- `lat::Matrix{Float64}`: lattice vector.
  Lattice vectors should be provided in columns.
- `positions::Matrix{Float64}`: atom positions in reduced coordinate.
  Atom positions should be provided in columns.

# Fields
- `norbits::Int`: number of orbits.
- `lat::Matrix{Float64}`: lattice vectors stored in columns.
- `rlat::Matrix{Float64}`: reciprocal lattice vectors stored in columns.
- `positions::Matrix{Float64}`: position of orbits in reduced coordinate stored in columns.
- `hoppings::Dict{Tuple{Int64,Int64,Array{Int64,1}}, Complex128}`: hoppings.
   Hopping example: (1, 1, [1, 0, 0]) => 1.0+0.0im indicates hopping from orbit 1
   in unit cell labeled by [1, 0, 0] to orbit 1 in home unit cell is 1.0.
"""
function TightBindingModel(lat::Matrix{Float64}, positions::Matrix{Float64})
    @assert size(lat) == (3, 3) "Size of lat is not correct."
    @assert size(positions, 1) == 3 "Size of positions is not correct."
    rlat = 2*π*inv(lat)'
    TightBindingModel(size(positions, 2), lat, rlat, positions, Dict())
end


"""
    sethopping!(t::TightBindingModel, n::Int64, m::Int64, R::Vector{Int64}, hopping::Number)

Set hoppings for a TightBindingModel t. Hoppings are labeled as ⟨0n|H|Rm⟩,
where R is a 3-element Vector{Int} representing lattice vector. Hamiltonian is
guaranteed to be Hermitian.
"""
function sethopping!(t::TightBindingModel, n::Int64, m::Int64, R::Vector{Int64}, hopping::Number)
    @assert (n in 1:t.norbits) && (m in 1:t.norbits) "No such orbit."
    t.hoppings[(n, m, R)] = hopping
    t.hoppings[(m, n, -R)] = hopping
    return
end


"""
    calhamiltonian(t::TightBindingModel, k::Vector{T}) where T<:Real --> Matrix{Complex128}

Calculate Hamiltonian of a TightBindingModel t for a specific k point. k should
be provided in reduced coordinate.
"""
function calhamiltonian(t::TightBindingModel, k::Vector{T}) where T<:Real
    @assert size(k) == (3,) "Size of k is not correct."
    h = zeros(Complex128, (t.norbits, t.norbits))
    for ((n, m, R), hopping) in t.hoppings
        h[n, m] += exp(2π*im*(k⋅R))*hopping
    end
    return h
end


"""
    caleig(t::TightBindingModel, k::Vector{<:Real}, calegvecs::Bool=false)

Calculate eigenvalues and eigenvectors of t. k should be provided in reduced coordinate.

# Return
If calegvecs is true, `(egvals::Vector{Float64}, egvecs::Matrix{Complex128})`,
otherwise just `egvals::Vector{Float64}`. Eigenvectors are stored in columns
and eigenvalues are sorted from small to large.
"""
function caleig(t::TightBindingModel, k::Vector{<:Real}, calegvecs::Bool=false)
    @assert size(k) == (3,) "Size of k is not correct."
    hamiltonian = calhamiltonian(t, k)
    if calegvecs
        (egvals, egvecs) = eig(hamiltonian)
        egvals = real(egvals)
        perm = sortperm(egvals)
        sortedegvecs = zeros(Complex128, size(egvecs))
        for i in 1:t.norbits
            sortedegvecs[:, i] = egvecs[:, perm[i]]
        end
        return (egvals[perm], sortedegvecs)
    else
        return sort(real(eigvals(hamiltonian)))
    end
end


"""
    makesupercell(t::TightBindingModel, scrdlat::Matrix{Int64}) --> TightBindingModel

Create a supercell out of a TightBindingModel t. scrdlat is a 3x3 matrix representing
supercell reduced lattice vector in columns.
"""
function makesupercell(t::TightBindingModel, scrdlat::Matrix{Int64})
    @assert size(scrdlat) == (3, 3) "Size of scrdlat is not correct."
    @assert det(scrdlat) > 0 "scrdlat is not right handed."
    # find all the unit cells in the supercell
    ucs = []
    # convert scrdlat to Rational{BigInt} to perform accurate linear
    # algebra calculations
    scrdlatinv = inv(convert(Matrix{Rational{BigInt}}, scrdlat))
    @assert typeof(scrdlatinv) == Matrix{Rational{BigInt}}
    for k in minimum(scrdlat[3, :]):maximum(scrdlat[3, :])
        for j in minimum(scrdlat[2, :]):maximum(scrdlat[2, :])
            for i in minimum(scrdlat[1, :]):maximum(scrdlat[1, :])
                sccoord = scrdlatinv*[i, j, k]
                @assert typeof(sccoord) == Vector{Rational{BigInt}}
                if all(sccoord.>=0) && all(sccoord.<1)
                    push!(ucs, [i, j, k])
                end
            end
        end
    end
    nucs = length(ucs)
    @assert nucs == det(scrdlat) "Number of unit cells found is not correct."
    # construct the supercell
    scpositions = zeros((3, nucs*t.norbits))
    for iuc in nucs
        for iorbit in 1:t.norbits
            scpositions[:, (iuc-1)*t.norbits+iorbit] = t.positions[:, iorbit]+ucs[iuc]
        end
    end
    scpositions = inv(scrdlat)*scpositions
    sc = TightBindingModel(t.lat*scrdlat, scpositions)
    # set hoppings
    for i in 1:nucs
        for ((n, m, R), hopping) in t.hoppings
            scR = Int64.(fld.(scrdlatinv*(ucs[i]+R), 1))
            sethopping!(
                sc,
                (i-1)*t.norbits+n,
                (indexin([ucs[i]+R-scrdlat*scR], ucs)[1]-1)*t.norbits+m,
                scR,
                hopping
            )
        end
    end
    return sc
end


"""
    cutedge(t::TightBindingModel, dir::Int64, glueedges::Bool=false) --> TightBindingModel

Create a D-1 dimensional TightBindingModel from a D dimensional one `t`. The finite
direction is represented by `dir` following the convention of 1:x, 2:y, 3:z.
If `glueedges` is true, the returned TightBindingModel will be made periodic in the
`dir` direction.
"""
function cutedge(t::TightBindingModel, dir::Int64, glueedges::Bool=false)
    r = deepcopy(t)
    for ((n, m, R), hopping) in r.hoppings
        if abs(R[dir]) > 0
            @assert abs(R[dir]) < 2  "Cutting an edge with glueedges=true is undefined for
                                     TightBindingModel with next nearest unit cell hopping
                                     is undefined in that direction."
            if glueedges
                pop!(r.hoppings, (n, m, R))
                sethopping!(
                r,
                n,
                m,
                [0, 0, 0],
                hopping
                )

            else
                pop!(r.hoppings, (n, m, R))
            end
        end
    end
    return r
end


"""
    addmagneticfield(t::TightBindingModel, B::Real) --> TightBindingModel

Add constant magnetic field in z direction for a TightBindingModel.
# Arguments
- `t::TightBindingModel`: a TightBindingModel.
- `B::Float64`: magnetic field in z direction. B is actually Be/h, thus
  its unit is 1/[length]^2. e here is fundamental charge.
  Since electron charge is -e, positive B means -z direction for electron system.
"""
function addmagneticfield(t::TightBindingModel, B::Real)
    tm = deepcopy(t)
    for ((n, m, R), hopping) in tm.hoppings
        # landau gauge
        absolute_position_n = tm.lat*tm.positions[:,n]
        absolute_position_m = tm.lat*tm.positions[:,m]
        tm.hoppings[(n, m, R)] = hopping*exp(
            im*2π*B*(absolute_position_n[2]-absolute_position_m[2])*(absolute_position_n[1]+absolute_position_m[1])/2
        )
    end
    return tm
end

end
