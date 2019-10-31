using DataFrames
using TimerOutputs

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
    timer::TimerOutput

    Iteration(
        mct::Type{<:MatrixCorrection},
        fct::Type{<:FluorescenceCorrection};
        updater=NaiveUpdateRule(),
        converged=RMSBelowTolerance(0.001),
        unmeasured=NullUnmeasuredRule()) =
        new(mct,fct,updater,converged,unmeasured,TimerOutput())
end

function ZAF(iter::Iteration, mat::Material, kr::KRatio)::MultiZAF
    coating = get(kr.stdProps, :Coating, NullCoating())
    return ZAF(iter.mctype, iter.fctype, mat, kr.lines, kr.stdProps[:BeamEnergy], coating)
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
    # Precompute ZAFs for std
    nc, stdZafs = NullCoating(), Dict{KRatio, MultiZAF}()
    for kr in measured
        coating = get(kr.stdProps, :Coating, nc)
        @timeit iter.timer "ZAF[std]" stdZafs[kr] = ZAF(iter, kr.standard, kr)
    end
    for kr in measured
        # Build ZAF for unk
        @timeit iter.timer "ZAF[unk]" unkZaf = ZAF(iter, est, kr)
        # Compute the total correction and the resulting k-ratio
        @timeit iter.timer "gZAFc" gzafc =  gZAFc(unkZaf, stdZafs[kr], kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle])
        estkrs[kr.element] = gzafc * est[kr.element] / kr.standard[kr.element]
    end
    return estkrs
end


struct IterationResult
    comp::Material
    kratios::Vector{KRatio}
    computed::Dict{Element, Float64}
    converged::Bool
    iterations::Int
end

function Base.show(io::IO, itres::IterationResult)
    print(itres.converged ? "Converged in $(itres.iterations) to $(itres.comp)\n" : "Failed to converge after $(itres.iterations).")
end

NeXLCore.compare(itres::IterationResult, known::Material)::DataFrame =
    compare(itres.comp, known)

NeXLCore.compare(itress::AbstractVector{IterationResult}, known::Material)::DataFrame =
    mapreduce(itres->compare(itres, known),append!,itress)

NeXLCore.material(itres::IterationResult) = itres.comp

mutable struct Counter
    count::Int
    terminate::Int
    Counter(terminate::Int) = new(0,terminate)
end

update(it::Counter)::Bool =
    (it.count+=1) <= it.terminate

terminated(it::Counter) =
    it.count > it.terminate

"""
    iterateks(iter::Iteration, name::String, measured::Vector{KRatio})

Iterate to find the composition that produces the measured k-ratios.
"""
function iterateks(iter::Iteration, name::String, measured::Vector{KRatio})::IterationResult
    @timeit iter.timer "FirstComp" estcomp = firstEstimate(iter, name, measured)
    @timeit iter.timer "FirstKs" estkrs = computeKs(iter, estcomp, measured)
    iters = Counter(100)
    while !converged(iter.converged, measured, estkrs) && update(iters)
        # println("$(estcomp) for $(estkrs)")
        @timeit iter.timer "NextEst" upd = update(iter.updater, estcomp, measured, estkrs)
        @timeit iter.timer "Unmeasured" unmeas = compute(iter.unmeasured, upd)
        @timeit iter.timer "Material" estcomp = material(name, unmeas)
        @timeit iter.timer "ComputeKs" estkrs = computeKs(iter, estcomp, measured)
    end
    return IterationResult(estcomp, measured, estkrs, !terminated(iters), iters.count)
end
