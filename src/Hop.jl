module Hop
using StaticArrays

export TightBindingModel, sethopping!, calhamiltonian, caleig, calband, makesupercell, cutedge, addmagneticfield


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
- `nspins::Int64=1`: 1 for spinless systems and 2 for spinful systems. If `nspin`
  is 2, `norbits` will be twice the number of `size(positions, 1)`. Orbits are
  ordered as (|1↑⟩, |1↓⟩, |2↑⟩, |2↓⟩, ...).

# Fields
- `norbits::Int`: number of orbits.
- `lat::Matrix{Float64}`: lattice vectors stored in columns.
- `rlat::Matrix{Float64}`: reciprocal lattice vectors stored in columns.
- `positions::Matrix{Float64}`: position of orbits in reduced coordinate stored in columns.
- `hoppings::Dict{Tuple{Int64,Int64,SVector{3,Int64}},Complex128}`: hoppings.
   Hopping example: (1, 1, [1, 0, 0]) => 1.0+0.0im indicates hopping from orbit 1
   in unit cell labeled by [1, 0, 0] to orbit 1 in home unit cell is 1.0.
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
    print(io, "TightBindingModel: $(t.norbits) orbitals")
end


"""
    sethopping!(t::TightBindingModel, n::Int64, m::Int64, R::Vector{Int64}, hopping)

Set ⟨0n|H|Rm⟩ to `hopping`. `hopping::Number` for spinless systems and
`hopping::Matrix{<:Number}` for spinful systems. For spinful systems,
`size(hopping)` should be (2, 2) and the basis for `hopping` is (|↑⟩, |↓⟩).
"""
function sethopping!(t::TightBindingModel, n::Int64, m::Int64, R::Union{Vector{Int64},SVector{3,Int64}}, hopping::Number)
    @assert (n in 1:t.norbits) && (m in 1:t.norbits) "No such orbit."
    t.hoppings[(n, m, R)] = hopping
    t.hoppings[(m, n, -R)] = conj(hopping)
    return
end


function sethopping!(t::TightBindingModel, n::Int64, m::Int64, R::Union{Vector{Int64},SVector{3,Int64}}, hopping::Matrix{<:Number})
    @assert (n in 1:Int64(t.norbits/2)) && (m in 1:Int64(t.norbits)) "No such orbit."
    @assert size(hopping) == (2, 2) "Size of hopping is not correct."
    sethopping!(t, 2n-1, 2m-1, R, hopping[1, 1])
    sethopping!(t, 2n, 2m-1, R, hopping[2, 1])
    sethopping!(t, 2n-1, 2m, R, hopping[1, 2])
    sethopping!(t, 2n, 2m, R, hopping[2, 2])
    return
end


"""
    calhamiltonian(t::TightBindingModel, k::Vector{<:Real}) --> Matrix{Complex128}

Calculate Hamiltonian of a TightBindingModel t for a specific k point. k should
be provided in reduced coordinate.
"""
function calhamiltonian(t::TightBindingModel, k::Vector{<:Real})
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
    calband(t::TightBindingModel, kpath::Matrix{<:Real}, ndiv::Int64) --> (Vector{Float64}, Matrix{Float64})

Calculate bands. `kpath` is a (3, x) size matrix where x is an even number and
should be provided in reduced coordinates.
This function returns (`kdist`, `egvals`). `kdist` is the distance of k points and
`egvals` is the energies of band stored in column at each k.
"""
function calband(t::TightBindingModel, kpath::Matrix{<:Real}, ndiv::Int64)
    @assert mod(size(kpath, 2), 2) == 0
    npaths = Int64(size(kpath, 2)/2)
    nkpts = Int64(size(kpath, 2)/2)*ndiv
    kdist = zeros(nkpts)
    egvals = zeros(t.norbits, nkpts)
    for ipath in 1:npaths
        dk = (kpath[:, 2*ipath]-kpath[:, 2*ipath-1])/(ndiv-1)
        dkn = norm(t.rlat*dk)/(2π) # real norm of dk
        if ipath == 1
            kdist0 = 0
        else
            kdist0 = kdist[(ipath-1)*ndiv]
        end
        for ikpt in 1:ndiv
            kdist[(ipath-1)*ndiv+ikpt] = dkn*(ikpt-1) + kdist0
            k = kpath[:, 2*ipath-1]+dk*(ikpt-1)
            egvals[:, (ipath-1)*ndiv+ikpt] = caleig(t, k)
        end
    end
    return (kdist, egvals)
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
    ucs = Vector{SVector{3, Int64}}()
    # convert scrdlat to Rational{BigInt} to perform accurate linear
    # algebra calculations
    scrdlatinv = inv(convert(SMatrix{3, 3, Rational{BigInt}}, scrdlat))
    ranges = [
        minimum(scrdlat[1, :]):maximum(scrdlat[1, :]),
        minimum(scrdlat[2, :]):maximum(scrdlat[2, :]),
        minimum(scrdlat[3, :]):maximum(scrdlat[3, :])
    ]
    @assert typeof(scrdlatinv) == SArray{Tuple{3,3},Rational{BigInt},2,9}
    for k in ranges[3]
        for j in ranges[2]
            for i in ranges[1]
                uc = @SVector [i, j, k]
                sccoord = scrdlatinv*uc
                @assert typeof(sccoord) == SVector{3, Rational{BigInt}}
                if all(sccoord.>=0) && all(sccoord.<1)
                    push!(ucs, uc)
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
    function findind(a, c)
        for i in 1:length(c)
            if a == c[i]
                return i
            end
        end
    end
    for i in 1:nucs
        for ((n, m, R), hopping) in t.hoppings
            scR = convert(SVector{3, Int64}, fld.(scrdlatinv*(ucs[i]+R), 1))
            sethopping!(
                sc,
                (i-1)*t.norbits+n,
                (findind(ucs[i]+R-scrdlat*scR, ucs)-1)*t.norbits+m,
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
