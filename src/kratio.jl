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
    kratio::Float64

    function KRatio(
        lines::Vector{CharXRay},
        unkProps::Dict{Symbol,<:Any},
        stdProps::Dict{Symbol,<:Any},
        standard::Material,
        kratio::Float64,
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
        return new(elm, lines, unkProps, stdProps, standard, kratio)
    end
end

NeXLCore.element(kr::KRatio) = kr.element
nonneg(kr::KRatio) = max(0.0, kr.kratio)

Base.show(io::IO, kr::KRatio) = print(io, "k[$(name(kr.standard)), $(name(kr.lines))] = $(kr.kratio)")

function NeXLCore.tabulate(krs::AbstractVector{KRatio})::DataFrame
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
        Material = mat,
        Celm = celm,
        KRatio = krv
    )
end
