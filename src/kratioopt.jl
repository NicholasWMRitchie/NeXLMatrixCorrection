"""
    KRatioOptimizer

Defines an optimizeks(kro::KRatioOptimizer, krs::Vector{KRatio})::Vector{KRatio} method which takes a vector
of k-ratios which may have redundant data (more than one KRatio per element) and trims it down to a vector
of k-ratios with one KRatio per element.
"""
abstract type KRatioOptimizer end

"""
    SimpleKRatioOptimizer(overvoltage, favor::Vector{CharXRay} = CharXRay[], minOvervoltage=1.2)

Implements a simple optimizer based on shell first, overvoltage next and brightness last.  Once it picks
an optimum set of lines for an element, it will not change.  You can force the selection of specific
lines by including the brightest in a family in `favor`.  Lines with overvoltages less than `minOvervoltage`
are not considered as available.
"""
struct SimpleKRatioOptimizer <: KRatioOptimizer
    overvoltage::Float64
    scores::Dict{Vector{CharXRay},Float64}
    favor::Vector{CharXRay}
    minOver::Float64

    function SimpleKRatioOptimizer(overvoltage, favor::Vector{CharXRay} = CharXRay[], minOvervoltage=1.2) 
        elms = map(element, favor)
        @assert length(unique(elms)) == length(elms) "SimpleKRatioOptimizer: Please specify only one characteristic X-ray per element to favor."
        new(overvoltage, Dict{Vector{CharXRay}, Float64}(), favor, minOvervoltage)
    end
end

Base.show(io::IO, skro::SimpleKRatioOptimizer) = print(io,"SimpleKRatioOptimizer[U>$(skro.overvoltage), favor=[$(skro.favor)]]")

function optimizeks(skro::SimpleKRatioOptimizer, krs::AbstractVector{T})::Vector{T} where  T <: Union{KRatio, KRatios}
    function score(kr) # Larger is better....
        br = brightest(kr.xrays)
        sc = any(f->f in kr.xrays, skro.favor) && (kr.unkProps[:BeamEnergy] / energy(inner(br)) > skro.minOver) ? # 
            1.0e100 : get(skro.scores, kr.xrays, -1.0)
        if sc == -1.0
            ov = min(kr.stdProps[:BeamEnergy], kr.unkProps[:BeamEnergy]) / energy(inner(br))
            sc = if ov > skro.minOver
                5.0 - n(shell(br)) - # Line K->8, L->6, M->4, N->2
                    skro.overvoltage / ov + # Overvoltage (<1 if ov > over)
                    0.1*sum(weight.(NormalizeToUnity, kr.xrays)) # line weight (favor brighter)
            else
                0.0 # Avoid it...
            end
            skro.scores[kr.xrays]=sc
        end
        return sc
    end
    return map(collect(elms(krs))) do elm
        elmkrs=filter(k->k.element==elm, krs)
        maxi = findmax(score.(elmkrs))
        maxi[1] == 0.0 && @warn "There is no k-ratio with an overvoltage greater than $(skro.minOver) for $elm."
        elmkrs[maxi[2]]
    end
end

"""
    NullKRatioOptimizer

Does nothing.  Returns the input k-ratios.  For when, the k-ratios are already optimal.
"""
struct NullKRatioOptimizer <: KRatioOptimizer end

function optimizeks(::NullKRatioOptimizer, krs::AbstractVector{T})::Vector{T} where  T <: Union{KRatio, KRatios}
    return krs
end

struct ElementalMap end

"""
    asa(ElementalMap, krs::AbstractVector{KRatios}, scale=Log3Band)

Converts `Vector{KRatios}` into a Dict that maps each Element present to an image.

    scale = [ Log3Band , Log3BandColorblind, Log3BandBright, LogScale, Gray ]

"""
function NeXLUncertainties.asa(::Type{ElementalMap}, krs::AbstractVector{KRatios}, scale=Log3Band, kro=SimpleKRatioOptimizer(2.0))
    opt=optimizeks(kro, krs)
    nks=normalizek(opt)
    Dict(map(nk-> nk.element => scale.(nk.kratios), nks))
end
