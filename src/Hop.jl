module Hop

export TightBindingModel, sethopping!, gethamiltonian, geteig

struct TightBindingModel
    norbits::Int
    lat::Matrix{Float64}
    rlat::Matrix{Float64}
    hoppings::Dict{Vector{Int64}, Float64}
end

function TightBindingModel(norbits, lat)
    a1 = lat[:, 1]
    a2 = lat[:, 2]
    a3 = lat[:, 3]
    rlat = zeros((3, 3))
    rlat[:, 1] = 2π*((a2×a3)/(a1⋅(a2×a3)))
    rlat[:, 2] = 2π*((a3×a1)/(a2⋅(a3×a1)))
    rlat[:, 3] = 2π*((a1×a2)/(a3⋅(a1×a2)))
    TightBindingModel(norbits, lat, rlat, Dict())
end

function sethopping!(t, n, m, R, hopping)
    t.hoppings[[[n, m]; R]] = hopping
    t.hoppings[[[m, n]; -R]] = hopping
    return
end

function gethamiltonian(t, k)
    h = zeros((t.norbits, t.norbits))
    for (label, hopping) in t.hoppings
        h[label[1], label[2]] += exp(2π*im*(k⋅label[3:5]))*hopping
    end
    return h
end

function geteig(t, k, cal_egvecs=false)
    hamiltonian = gethamiltonian(t, k)
    if cal_egvecs
        (egvals, egvecs) = eig(hamiltonian)
        egvals = real(eigvals)
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
