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
        kratio::Float64
    )
        if length(lines)<1
            error("Must specify at least one characteristic X-ray.")
        end
        elm = element(lines[1])
        if !all(element(l)==elm for l in lines)
            error("The characteristic X-rays must all be from the same element.")
        end
        if standard[elm]<=1.0e-4
            error("The standard must contain the element $(elm).  $(standard[elm])")
        end
        return new(elm, lines, unkProps,stdProps, standard, kratio)
    end
end

NeXLCore.element(kr::KRatio) = kr.element
nonneg(kr::KRatio) = max(0.0, kr.kratio)

Base.show(io::IO, kr::KRatio) =
    print(io, "k[$(name(kr.lines))] = $(kr.kratio)")


abstract type UnmeasuredElementRule end

struct NullUnmeasuredRule <: UnmeasuredElementRule end
"""
    compute(::Type{UnmeasuredElementRule}, inp::Dict{Element,Float64})::Dict{Element,Float64}

A null UnmeasuredElementRule.  Just returns the inputs.
"""
compute(::NullUnmeasuredRule, inp::Dict{Element,Float64})::Dict{Element,Float64} =
    inp

abstract type UpdateRule end

struct NaiveUpdateRule <: UpdateRule end

"""
    update(
        ::NaiveUpdateRule,
        prevcomp::Material,
        measured::Vector{KRatio},
        estkrs::Dict{Element, Float64}
    )::Dict{Element,Float64}

Determine the next estimate of the composition that brings the estkrs closer to measured.
"""
function update(::NaiveUpdateRule, prevcomp::Material, measured::Vector{KRatio}, estkrs::Dict{Element, Float64})::Dict{Element,Float64}
    emf = Dict{Element, Float64}()
    for mkr in measured
        emf[mkr.element] = (nonneg(mkr) / estkrs[mkr.element])*prevcomp[mkr.element]
    end
    return emf
end

abstract type ConvergenceTest end

struct RMSBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(rbt::RMSBelowTolerance, meas::Vector{KRatio}, comp::Dict{Element,Float64})::Bool =
    sum( (nonneg(kr)-comp[kr.element])^2 for kr in meas ) < rbt.tolerance^2

struct AllBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(abt::AllBelowTolerance, meas::Vector{KRatio}, comp::Dict{Element,Float64})::Bool =
    all(abs(nonneg(kr)-comp[kr.element])<abt.tolerance for kr in meas)

struct IsApproximate <: ConvergenceTest
    atol::Float64
    rtol::Float64
end

converged(ia::IsApproximate, meas::Vector{KRatio}, comp::Dict{Element,Float64}) =
    all( (abs(1.0 - nonneg(kr)/comp[kr.element])<rtol) || (abs(nonneg(kr)-comp[kr.element]) < atol) for kr in meas)

struct Iteration
    mctype::Type{<:MatrixCorrection}
    fctype::Type{<:FluorescenceCorrection}
    updater::UpdateRule
    converged::ConvergenceTest
    unmeasured::UnmeasuredElementRule
end

"""
    firstEstimate(iter::Iteration)::Material

Make a first estimate at the composition.
"""
function firstEstimate(iter::Iteration, name::String, measured::Vector{KRatio})::Material
    mfs = Dict{Element,Float64}( (kr.element, nonneg(kr) * kr.standard[kr.element]) for kr in measured )
    return material(name, compute(iter.unmeasured, mfs))
end

"""
    computeKs(iter::Iteration, est::Material)::Dict{Element, Float64}

Given an estimate of the composition compute the corresponding k-ratios.
"""
function computeKs(iter::Iteration, est::Material, measured::Vector{KRatio})::Dict{Element, Float64}
    estkrs = Dict{Element,Float64}()
    for kr in measured
        # Build ZAF for std
        coating = get(kr.stdProps, :Coating, NullCoating())
        stdZaf = ZAF(iter.mctype, iter.fctype, kr.standard, kr.lines, kr.stdProps[:BeamEnergy], coating)
        # Build ZAF for unk
        coating = get(kr.unkProps, :Coating, NullCoating())
        unkZaf = ZAF(iter.mctype, iter.fctype, est, kr.lines, kr.unkProps[:BeamEnergy], coating)
        # Compute the total correction and the resulting k-ratio
        gzafc =  gZAFc(unkZaf, stdZaf, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle])
        estkrs[kr.element] = gzafc * est[kr.element] / kr.standard[kr.element]
    end
    return estkrs
end

"""
    iterate(iter::Iteration, name::String, measured::Vector{KRatio})

Iterate to find the composition that produces the measured k-ratios.
"""
function Base.iterate(iter::Iteration, name::String, measured::Vector{KRatio})::Material
    estcomp = firstEstimate(iter, name, measured)
    estkrs = computeKs(iter, estcomp, measured)
    while !converged(iter.converged, measured, estkrs)
        # println("$(estcomp) for $(estkrs)")
        estcomp = material(name, compute(iter.unmeasured, update(iter.updater, estcomp, measured, estkrs)))
        estkrs = computeKs(iter, estcomp, measured)
    end
    return estcomp
end
