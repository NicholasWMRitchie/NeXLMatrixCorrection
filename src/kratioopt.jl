"""
    KRatioOptimizer

Defines an optimizeks(kro::KRatioOptimizer, krs::Vector{KRatio})::Vector{KRatio} method which takes a vector
of k-ratios which may have redundant data (more than one KRatio per element) and trims it down to a vector
of k-ratios with one KRatio per element.
"""
abstract type KRatioOptimizer end

"""
    SimpleKRatioOptimizer(overvoltage, favor::Vector{CharXRay} = CharXRay[])

Implements a simple optimizer based on shell first, overvoltage next and brightness last.  Once it picks
an optimum set of lines for an element, it will not change.  You can force the selection of specific
lines by including the brightest in a family in `favor`.
"""
struct SimpleKRatioOptimizer <: KRatioOptimizer
    overvoltage::Float64
    scores::Dict{Vector{CharXRay},AbstractFloat}
    favor::Vector{CharXRay}

    SimpleKRatioOptimizer(overvoltage, favor::Vector{CharXRay} = CharXRay[]) = new(overvoltage, Dict{Vector{CharXRay}, AbstractFloat}(), favor)
end

Base.show(io::IO, skro::SimpleKRatioOptimizer) = print(io,"SimpleKRatioOptimizer[U>$(skro.overvoltage), favor=$(skro.favor)")

function optimizeks(skro::SimpleKRatioOptimizer, krs::AbstractVector{T})::Vector{T} where  T <: Union{KRatio, KRatios}
    function score(kr) # Larger is better....
        if brightest(kr.xrays) in skro.favor
            sc = 1.0e100
        else
            sc =  get(skro.scores, kr.xrays, -1.0)
        end
        if sc==-1.0
            br = brightest(kr.xrays)
            ov = min(kr.stdProps[:BeamEnergy], kr.unkProps[:BeamEnergy]) / edgeenergy(br)
            sc = convert(Float64, 5 - n(shell(br))) - # Line K->4, L->3, M->2, N->1
                skro.overvoltage / ov + # Overvoltage (<1 if ov > over)
                0.1*sum(weight.(kr.xrays)) # line weight (favor brighter)
            skro.scores[kr.xrays]=sc
        end
        return sc
    end
    res = Vector{T}()
    for elm in elms(krs)
        elmkrs=filter(k->k.element==elm, krs)
        push!(res, elmkrs[findmax(score.(elmkrs))[2]])
    end
    return res
end
