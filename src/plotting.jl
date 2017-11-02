module Plotting

using PyPlot
using ..Hop

"""
```julia
plotwf(t::TightBindingModel, wf::Dict{Vector{Int64}, Matrix{T}}, iwf::Int64;
    scale::Float64=100.0) where T<:Number
```

Plot wannier functions indexed by `iwf`.
"""
function plotwf(t::TightBindingModel, wf::Dict{Vector{Int64}, Matrix{T}}, iwf::Int64; scale::Float64=100.0) where T<:Number
    for (R, overlap) in wf
        for iorbit in 1:t.norbits
            pos = (t.lat*(R+t.positions[:, iorbit]))
            scatter(pos[1], pos[2], s=abs(overlap[iorbit, iwf])*scale, c="k")
        end
    end
end

end
