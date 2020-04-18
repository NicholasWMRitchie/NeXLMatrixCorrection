"""
KRatioOptimizer abstract type

Defines an optimizeks(kro::KRatioOptimizer, krs::Vector{KRatio})::Vector{KRatio} method which takes a vector
of k-ratios which may have redundant data (more than one KRatio per element) and trims it down to a vector
of k-ratios with one KRatio per element.
"""
abstract type KRatioOptimizer end

"""
    SimpleKRatioOptimizer

Implements a simple optimizer based on shell first, overvoltage next and brightness last.  Once it picks
an optimum set of lines for an element, it will not change.
"""
struct SimpleKRatioOptimizer <: KRatioOptimizer
    overvoltage::Float64
    scores::Dict{Vector{CharXRay},AbstractFloat}
    SimpleKRatioOptimizer(overvoltage) = new(overvoltage, Dict{Vector{CharXRay}, AbstractFloat}())
end

function optimizeks(skro::SimpleKRatioOptimizer, krs::Vector{KRatio})::Vector{KRatio}
    function score(kr) # Larger is better....
        sc = get(skro.scores, kr.lines, -1.0)
        if sc==-1.0
            br = brightest(kr.lines)
            ov = min(kr.stdProps[:BeamEnergy], kr.unkProps[:BeamEnergy]) / edgeenergy(br)
            sc = convert(Float64, 5 - n(shell(br))) - # Line K->4, L->3, M->2, N->1
                skro.overvoltage / ov + # Overvoltage (<1 if ov > over)
                0.1*sum(weight.(kr.lines)) # line weight (favor brighter)
            skro.scores[kr.lines]=sc
        end
        return sc
    end
    res = Vector{KRatio}()
    for elm in elms(krs)
        elmkrs=filter(k->k.element==elm, krs)
        push!(res, elmkrs[findmax(score.(elmkrs))[2]])
    end
    return res
end
