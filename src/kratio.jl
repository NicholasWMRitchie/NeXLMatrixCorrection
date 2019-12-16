using NeXLUncertainties
using NeXLCore

"""
    KRatio

The k-ratio is the result of two intensity measurements - one on a standard
with known composition and one on an unknown. Each measurement has properties
like :BeamEnergy (req), :TakeOffAngle (req), :Coating (opt) that characterize
the measurement.

Properties: (These Symbols are intentionally the same used in NeXLSpectrum)

    :BeamEnergy incident beam energy in eV
    :TakeOffAngle in radians
    :Coating A NeXLCore.Film object describing a conductive coating
"""
struct KRatio
    element::Element
    lines::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratio::UncertainValue

    function KRatio(
        lines::Vector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::AbstractFloat
    )
        if length(lines) < 1
            error("Must specify at least one characteristic X-ray.")
        end
        elm = element(lines[1])
        if !all(element(l) == elm for l in lines)
            error("The characteristic X-rays must all be from the same element.")
        end
        if standard[elm] <= 1.0e-4
            error("The standard must contain the element $(elm).  $(standard[elm])")
        end
        return new(elm, lines, unkProps, stdProps, standard, convert(UncertainValue, kratio))
    end
end

NeXLCore.element(kr::KRatio) = kr.element
nonnegk(kr::KRatio) = max(zero(typeof(kr.kratio), kr.kratio))

Base.show(io::IO, kr::KRatio) = print(io, "k[$(name(kr.standard)), $(name(kr.lines))] = $(kr.kratio)")

function NeXLUncertainties.asa(::Type{DataFrame}, krs::AbstractVector{KRatio})::DataFrame
    elms, zs, lines, e0u = Vector{String}(), Vector{Int}(), Vector{Vector{CharXRay}}(), Vector{Float64}()
    e0s, toau, toas, mat = Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{String}()
    celm, krv = Vector{Float64}(), Vector{Float64}()
    for kr in krs
        push!(elms, kr.element.symbol)
        push!(zs, z(kr.element))
        push!(lines, kr.lines)
        push!(e0u, get(kr.unkProps, :BeamEnergy, -1.0))
        push!(e0s, get(kr.stdProps, :BeamEnergy, -1.0))
        push!(toau, get(kr.unkProps, :TakeOffAngle, -1.0))
        push!(toas, get(kr.stdProps, :TakeOffAngle, -1.0))
        push!(mat, name(kr.standard))
        push!(celm, kr.standard[kr.element])
        push!(krv, kr.kratio)
    end
    return DataFrame(
        Element = elms,
        Z = zs,
        Lines = lines,
        E0unk = e0u,
        E0std = e0s,
        θunk = toau,
        θstd = toas,
        Standard = mat,
        Cstd = celm,
        KRatio = krv
    )
end

"""
    elms(krs::Vector{KRatio})::Vector{Element}

Returns a vector containing the elements present in krs (no duplicate elements).
"""
function NeXLCore.elms(krs::Vector{KRatio})::Vector{Element}
    res=Vector{Element}()
    for kr in krs
        if !(kr.element in res)
            push!(res, kr.element)
        end
    end
    return res
end

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
            sc = convert(Float64, 'O'-shell(br)) - # Line K->4, L->3, M->2, N->1
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
