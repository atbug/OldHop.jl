module Hop
using StaticArrays

export TightBindingModel, sethopping!, calhamiltonian, caleig, calband,
    makesupercell, cutedge, addmagneticfield, calproj, calwf, calwilson


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
    hoppings::Dict{SVector{3,Int64},Matrix{Complex128}}
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
- `spinful::Bool=false`: false for spinless systems and true for spinful systems. If `spinful`
  is true, `norbits` will be twice the number of `size(positions, 1)`. Orbits are
  ordered as (|1↑⟩, |1↓⟩, |2↑⟩, |2↓⟩, ...).

# Fields
- `norbits::Int`: number of orbits.
- `lat::Matrix{Float64}`: lattice vectors stored in columns.
- `rlat::Matrix{Float64}`: reciprocal lattice vectors stored in columns.
- `positions::Matrix{Float64}`: position of orbits in reduced coordinate stored in columns.
- `hoppings::Dict{SVector{3,Int64},Matrix{Complex128}}`: hoppings stored as
   R->⟨0n|H|Rm⟩.
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
```julia
sethopping!(t::TightBindingModel, n::Int64, m::Int64, R::Vector{Int64}, hopping;
    mode::Char='a')
```

Set ⟨0n|H|Rm⟩ to `hopping`. `hopping::Number` for spinless systems and
`hopping::Matrix{<:Number}` for spinful systems. For spinful systems,
`size(hopping)` should be (2, 2) and the basis for `hopping` is (|↑⟩, |↓⟩).
`mode` has two possible values: 'a' for add mode and 's' for set or reset mode.
"""
function sethopping!(t::TightBindingModel, n::Int64, m::Int64,
    R::Union{Vector{Int64},SVector{3,Int64}}, hopping::Number; mode::Char='a')
    @assert (n in 1:t.norbits) && (m in 1:t.norbits) "No such orbit."
    if !(R in keys(t.hoppings))
        @assert !(-R in keys(t.hoppings))
        t.hoppings[R] = zeros(Complex128, (t.norbits, t.norbits))
        t.hoppings[-R] = zeros(Complex128, (t.norbits, t.norbits))
    end
    if mode == 's'
        t.hoppings[R][n, m] = hopping
        t.hoppings[-R][m, n] = conj(hopping)
    elseif mode == 'a'
        t.hoppings[R][n, m] += hopping
        if (n, m, R) != (m, n, -R) # not onsite energy
            t.hoppings[-R][m, n] += conj(hopping)
        end
    end
    return
end


# For spinful system
function sethopping!(t::TightBindingModel, n::Int64, m::Int64,
    R::Union{Vector{Int64},SVector{3,Int64}}, hopping::Matrix{<:Number}; mode::Char='a')
    @assert (n in 1:Int64(t.norbits/2)) && (m in 1:Int64(t.norbits)) "No such orbit."
    @assert size(hopping) == (2, 2) "Size of hopping is not correct."
    sethopping!(t, 2n-1, 2m-1, R, hopping[1, 1], mode=mode)
    sethopping!(t, 2n, 2m-1, R, hopping[2, 1], mode=mode)
    if (n, m, R) != (m, n, -R) # not onsite energy
        sethopping!(t, 2n-1, 2m, R, hopping[1, 2], mode=mode)
    end
    sethopping!(t, 2n, 2m, R, hopping[2, 2], mode=mode)
    return
end


"""
```julia
calhamiltonian(t::TightBindingModel, k::Vector{<:Real}) --> Matrix{Complex128}
```

Calculate Hamiltonian of a TightBindingModel t for a specific k point. k should
be provided in reduced coordinate.
"""
function calhamiltonian(t::TightBindingModel, k::Vector{<:Real})
    @assert size(k) == (3,) "Size of k is not correct."
    h = zeros(Complex128, (t.norbits, t.norbits))
    for (R, hopping) in t.hoppings
        h += exp(2π*im*(k⋅R))*hopping
    end
    return h
end


"""
```julia
caleig(t::TightBindingModel, k::Vector{<:Real}, calegvecs::Bool=false)
```

Calculate eigenvalues and eigenvectors of t. k should be provided in reduced coordinate.

# Return
If calegvecs is true, `(egvals::Vector{Float64}, egvecs::Matrix{Complex128})`,
otherwise just `egvals::Vector{Float64}`. Eigenvectors are stored in columns
and eigenvalues are sorted from small to large.
"""
function caleig(t::TightBindingModel, k::Vector{<:Real}; calegvecs::Bool=false)
    @assert size(k) == (3,) "Size of k is not correct."
    hamiltonian = calhamiltonian(t, k)
    if calegvecs
        (egvals, egvecs) = eig(hamiltonian)
        egvals = real(egvals)
        perm = sortperm(egvals)
        return (egvals[perm], egvecs[:, perm])
    else
        return sort(real(eigvals(hamiltonian)))
    end
end


"""
```julia
calband(t::TightBindingModel, kpath::Matrix{<:Real}, ndiv::Int64)
    --> (Vector{Float64}, Matrix{Float64})
```

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
```julia
makesupercell(t::TightBindingModel, scrdlat::Matrix{Int64}) --> TightBindingModel
```

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
    for iuc in 1:nucs
        for (R, hopping) in t.hoppings
            scR = convert(SVector{3, Int64}, fld.(scrdlatinv*(ucs[iuc]+R), 1))
            scmbase = (findind(ucs[iuc]+R-scrdlat*scR, ucs)-1)*t.norbits
            for n in 1:t.norbits
                for m in 1:t.norbits
                    sethopping!(
                        sc,
                        (iuc-1)*t.norbits+n,
                        scmbase+m,
                        scR,
                        hopping[n, m],
                        mode='s'
                    )
                end
            end
        end
    end
    return sc
end


"""
```julia
cutedge(t::TightBindingModel, dir::Int64; glueedges::Bool=false)
    --> TightBindingModel
```

Create a D-1 dimensional TightBindingModel from a D dimensional one `t`. The finite
direction is represented by `dir` following the convention of 1:x, 2:y, 3:z.
If `glueedges` is true, the returned TightBindingModel will be made periodic in the
`dir` direction.
"""
function cutedge(t::TightBindingModel, dir::Int64; glueedges::Bool=false)
    r = deepcopy(t)
    for (R, hopping) in r.hoppings
        if abs(R[dir]) > 0
            if glueedges
                @assert abs(R[dir]) < 2  "Cutting an edge with glueedges=true is undefined for
                                         TightBindingModel with next nearest unit cell hopping
                                         is undefined in that direction."
                pop!(r.hoppings, R)
                if R[dir] > 0
                    newR = convert(Vector{Int64}, R)
                    newR[dir] = 0
                    for n in 1:t.norbits
                        for m in 1:t.norbits
                            sethopping!(
                            r,
                            n,
                            m,
                            newR,
                            hopping[n, m],
                            mode='a'
                            )
                        end
                    end
                end
            else
                pop!(r.hoppings, R)
            end
        end
    end
    return r
end


"""
```julia
addmagneticfield(t::TightBindingModel, B::Real) --> TightBindingModel
```

Add constant magnetic field in z direction for a TightBindingModel.
# Arguments
- `t::TightBindingModel`: a TightBindingModel.
- `B::Float64`: magnetic field in z direction. B is actually Be/h, thus
  its unit is 1/[length]^2. e here is fundamental charge.
  Since electron charge is -e, positive B means -z direction for electron system.
"""
function addmagneticfield(t::TightBindingModel, B::Real)
    tm = deepcopy(t)
    for (R, hopping) in tm.hoppings
        for n in 1:t.norbits
            for m in 1:t.norbits
                # landau gauge
                absolute_position_n = tm.lat*tm.positions[:,n]
                absolute_position_m = tm.lat*tm.positions[:,m]
                tm.hoppings[R][n, m] = hopping[n, m]*exp(
                    im*2π*B*(absolute_position_n[2]-absolute_position_m[2])*(absolute_position_n[1]+absolute_position_m[1])/2
                )
            end
        end
    end
    return tm
end


"""
```julia
calproj(t::TightBindingModel, lfs::Dict{Vector{Int64}, Matrix{T1}},
    bands::Vector{Int64}, k::Vector{T2}) where T1<:Number where T2<:Real
    --> Matrix{Complex{Float64}}
```

Calculate ⟨u_nk|g_m⟩, where g is a localized function.
`lfs` is stored in format {R: ⟨Rn|gm⟩}.
"""
function calproj(t::TightBindingModel, lfs::Dict{Vector{Int64}, Matrix{T1}},
    bands::Vector{Int64}, k::Vector{T2}) where T1<:Number where T2<:Real
    nlfs = 0
    for (R, overlap) in lfs
        nlfs = size(overlap, 2)
        break
    end
    proj = zeros(Complex128, (length(bands), nlfs))
    egvecs = caleig(t, k; calegvecs=true)[2][:, bands]
    for (R, overlap) in lfs
        proj += exp(-im*2π*(k⋅R))*(egvecs')*overlap
    end
    return proj
end



function _smooth(t, lfs::Dict{Vector{Int64}, Matrix{T1}},
    bands::Vector{Int64}, k::Vector{T2}) where T1<:Number where T2<:Real
    A = calproj(t, lfs, bands, k)
    S = A'*A
    @assert real(det(S)) > 0.1
    egvecs = caleig(t, k, calegvecs=true)[2]
    return egvecs[:, bands]*A*inv(sqrtm(S))
end


"""
```julia
calwf(t::TightBindingModel, twfs::Dict{Vector{Int64}, Matrix{T}},
    bands::Vector{Int64}, nkmesh::Vector{Int64}, nrmesh::Vector{Int64})
    where T<:Number --> Dict{Vector{Int64}, Matrix{Complex128}}
```

Calculate wannier functions of `bands`.
"""
function calwf(t::TightBindingModel, twfs::Dict{Vector{Int64}, Matrix{T}},
    bands::Vector{Int64}, nkmesh::Vector{Int64}, nrmesh::Vector{Int64}) where T<:Number
    @assert size(nkmesh, 1) == 3
    @assert size(nrmesh, 1) == 3
    # generate number of Wannier functions
    nwfs = 0
    for (R, overlap) in twfs
        nwfs = size(overlap, 2)
        break
    end
    @assert nwfs <= length(bands)
    # generate k points
    nkpts = prod(nkmesh)
    kpts = zeros(3, nkpts)
    i = 1
    for ikx in 1:nkmesh[1]
        for iky in 1:nkmesh[2]
            for ikz in 1:nkmesh[3]
                kpts[:, i] = [(ikx-1)/nkmesh[1], (iky-1)/nkmesh[2], (ikz-1)/nkmesh[3]]
                i += 1
            end
        end
    end
    # generate r points
    nrpts = prod(2*nrmesh+1)
    wf = Dict{Vector{Int64}, Matrix{Complex128}}()
    i = 1
    for Rx in (-nrmesh[1]):nrmesh[1]
        for Ry in (-nrmesh[2]):nrmesh[2]
            for Rz in (-nrmesh[3]):nrmesh[3]
                wf[[Rx, Ry, Rz]] = zeros(Complex128, (t.norbits, nwfs))
                i += 1
            end
        end
    end
    # perform integration
    for ik in 1:nkpts
        Utilde = _smooth(t, twfs, bands, kpts[:, ik])
        for R in keys(wf)
            wf[R] += exp(im*2π*(kpts[:, ik]⋅R))*Utilde
        end
    end

    # normalize Wannier function
    N = zeros(nwfs)
    for R in keys(wf)
        N += diag(wf[R]'*wf[R])
    end
    for R in keys(wf)
        wf[R] = wf[R]./reshape(sqrt.(N), (1, nwfs))
    end

    return wf
end


"""
```julia
calwilson(t::TightBindingModel, bands::Vector{Int64}, kpath::Matrix{<:Real},
    ndiv::Int64) --> Vector{Float64}
```

Calculate Wilson loop. k points in `k_path` are stored in column. `kpath` must
be closed.
"""
function calwilson(t::TightBindingModel, bands::Vector{Int64}, kpath::Matrix{<:Real}, ndiv::Int64)
    @assert ndiv > 1
    bands = sort(bands)
    W = eye(Complex128, length(bands))
    npath = size(kpath, 2)-1

    klist = zeros(3, npath*ndiv+1)
    for ipath in 1:npath
        dk = (kpath[:, ipath+1]-kpath[:, ipath])/ndiv
        for ik in 1:ndiv
            klist[:, (ipath-1)*ndiv+ik] = kpath[:, ipath]+dk*(ik-1)
        end
    end
    klist[:, end] = kpath[:, end]

    egvecski = caleig(t, kpath[:, 1], calegvecs=true)[2]
    egvecsk1 = egvecski
    egvecsk2 = egvecski
    for ik in 1:npath*ndiv
        egvecsk2 = caleig(t, klist[:, ik+1], calegvecs=true)[2]
        dk = klist[:, ik+1]-klist[:, ik]
        rdk = diagm(exp.((im*2π*dk'*(t.positions)))')
        W = egvecsk2[:, bands]'*rdk*(egvecsk1[:, bands])*W
        egvecsk1 = egvecsk2
    end
    W = egvecski[:, bands]'*(egvecsk2[:, bands])*W

    return sort(imag.(log.(eigvals(W))))
end

include("plotting.jl")

end
