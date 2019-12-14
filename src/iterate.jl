using DataFrames
using TimerOutputs

abstract type UnmeasuredElementRule end

struct NullUnmeasuredRule <: UnmeasuredElementRule end
"""
    compute(::Type{UnmeasuredElementRule}, inp::Dict{Element,Float64})::Dict{Element,Float64}

A null UnmeasuredElementRule.  Just returns the inputs.
"""
compute(::NullUnmeasuredRule, inp::Dict{Element,Float64})::Dict{Element,Float64} = inp

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
function update(::NaiveUpdateRule, #
    prevcomp::Material, #
    measured::Vector{KRatio}, #
    estkrs::Dict{Element,Float64}, #
)::Dict{Element, Float64}
    cnp1 = Dict{Element,Float64}()
    for mkr in measured
        cnp1[mkr.element] = (nonnegk(mkr) / estkrs[mkr.element]) * prevcomp[mkr.element]
    end
    return cnp1
end

function reset(::NaiveUpdateRule) end

struct WegsteinUpdateRule <: UpdateRule
    prevc::Vector{Material}
    prevk::Vector{Dict{Element,Float64}}
    factor::Float64
    norm::Vector{Float64}
    WegsteinUpdateRule(f::Float64=0.3) = new(Vector{Material}(), Vector{Dict{Element,Float64}}(),f,Float64[])
end

function update( #
    weg::WegsteinUpdateRule,
    prevcomp::Material,
    measured::Vector{KRatio},
    estkrs::Dict{Element,Float64},
)::Dict{Element,Float64}
    bound(x, min, max) = x < min ? min : (x > max ? max : x)
    cnp1, itercx = Dict{Element,Float64}(), length(weg.prevc)
    if itercx < 1
        push!(weg.norm, 1.0)
        for mkr in measured
            cnp1[mkr.element] = estkrs[mkr.element] > 0.0 ?
                                (nonnegk(mkr) / estkrs[mkr.element]) * prevcomp[mkr.element] : 0.0
        end
    else
        safediv(a, b) = isapprox(b, 0.0, atol = 1.0e-12) ? 1.0 : a / b
        push!(weg.norm, itercx % 5 ≠ 0 ? weg.norm[end] : #
            weg.norm[end] / sum(safediv(nonnegk(mkr), estkr[mkr.element]) * prevcomp[mkr.element] for mkr in measured))
        cn, cnm1, kn, knm1 = NeXLCore.asnormalized(prevcomp, weg.norm[end]), weg.prevc[end], estkrs, weg.prevk[end]
        for mkr in measured
            elm, km = mkr.element, nonnegk(mkr)
            if (kn[elm] > 0) && (knm1[elm] > 0)
                fcn, fcnm1 = cn[elm] / kn[elm], cnm1[elm] / knm1[elm] # c = k*f
                dfdk = ((fcn - fcnm1) / (cn[elm] - cnm1[elm]))
                cnp1[elm] = cn[elm] + (km * fcn - cn[elm]) / (1.0 - bound(km * dfdk, -weg.factor, weg.factor))
            else
                cnp1[elm] = 0.0
            end
        end
    end
    push!(weg.prevc, prevcomp)
    push!(weg.prevk, estkrs)
    return cnp1
end

function reset(weg::WegsteinUpdateRule)
    resize!(weg.prevc, 0)
    resize!(weg.prevk, 0)
end

struct RecordingUpdateRule <: UpdateRule
    base::UpdateRule
    estkrs::Vector{Dict{Element,Float64}}
    comps::Vector{Dict{Element,Float64}}
    prev::Vector{Material}
    meas::Dict{Element, KRatio}
    RecordingUpdateRule(ur::UpdateRule) =
        new(ur,Vector{Dict{Element,Float64}}(),Vector{Dict{Element,Float64}}(),Vector{Material}(),Dict{Element,KRatio}())
end

function NeXLUncertainties.asa(::Type{DataFrame}, rur::RecordingUpdateRule)
    dkrs, dcs, prev = Dict{Element, Vector{Float64}}(), Dict{Element,Vector{Float64}}(), Dict{Element,Vector{Float64}}()
    allelms = union(keys(rur.estkrs[1]), keys(rur.comps[1]), keys(rur.prev[1]), keys(rur.meas))
    for elm in allelms
        dkrs[elm], dcs[elm], prev[elm], meas[elm] = [], [], [], []
    end
    for i in eachindex(rur.estkrs)
        for elm in allelms
            push!(dkrs[elm], get(rur.estkrs[i], elm, 0.0))
            push!(dcs[elm], get(rur.comps[i], elm, 0.0))
            push!(prev[elm], rur.prev[i][elm])
            push!(meas[elm], rur.meas[elm])
        end
    end
    df = DataFrame(Iter=collect(eachindex(rur.estkrs)))
    for elm in allelms
        df[!, Symbol("Prev($(elm.symbol))")]=prev[elm]
        df[!, Symbol("Next($(elm.symbol))")]=dcs[elm]
        df[!, Symbol("kest($(elm.symbol))")]=dkrs[elm]
        df[!, Symbol("meas($(elm.symbol))")]=meas[elm]
    end
    return df
end

function NeXLMatrixCorrection.update( #
    rur::RecordingUpdateRule,
    prevcomp::Material,
    measured::Vector{KRatio},
    estkrs::Dict{Element,Float64},
)::Dict{Element,Float64}
    res = NeXLMatrixCorrection.update(rur.base, prevcomp, measured, estkrs)
    if isempty(rur.meas)
        append!(rur.meas, measured)
    end
    push!(rur.prev, prevcomp)
    push!(rur.estkrs, estkrs)
    push!(rur.comps, res)
    return res
end

function NeXLMatrixCorrection.reset(rur::RecordingUpdateRule)
    resize!(rur.estkrs, 0)
    resize!(rur.comps, 0)
    resize!(rur.prev, 0)
    resize!(rur.meas, 0)
    NeXLMatrixCorrection.reset(rur.base)
end

abstract type ConvergenceTest end

struct RMSBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(rbt::RMSBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool =
    sum((nonnegk(kr) - computed[kr.element])^2 for kr in meas) < rbt.tolerance^2

struct AllBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(abt::AllBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool =
    all(abs(nonnegk(kr) - computed[kr.element]) < abt.tolerance for kr in meas)

struct IsApproximate <: ConvergenceTest
    atol::Float64
    rtol::Float64
end

converged(ia::IsApproximate, meas::Vector{KRatio}, computed::Dict{Element,Float64}) =
    all((abs(1.0 - nonnegk(kr) / computed[kr.element]) < rtol) || (abs(nonnegk(kr) - computed[kr.element]) < atol) for kr in meas)

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
        updater = WegsteinUpdateRule(),
        converged = RMSBelowTolerance(0.0001),
        unmeasured = NullUnmeasuredRule(),
    ) = new(mct, fct, updater, converged, unmeasured, TimerOutput())
end

function _ZAF(iter::Iteration, mat::Material, props::Dict{Symbol,Any}, lines::Vector{CharXRay})::MultiZAF
    coating = get(props, :Coating, NullCoating())
    return ZAF(iter.mctype, iter.fctype, mat, lines, props[:BeamEnergy], coating)
end


"""
    firstEstimate(iter::Iteration)::Material

Make a first estimate at the composition.
"""
function firstEstimate(iter::Iteration, name::String, measured::Vector{KRatio})::Material
    mfs = Dict{Element,Float64}((kr.element, nonnegk(kr) * kr.standard[kr.element]) for kr in measured)
    return material(name, compute(iter.unmeasured, mfs))
end

struct IterationResult
    comp::Material
    kratios::Vector{KRatio}
    computed::Dict{Element,Float64}
    converged::Bool
    iterations::Int
end

function Base.show(io::IO, itres::IterationResult)
    print(itres.converged ? "Converged in $(itres.iterations) to $(itres.comp)\n" :
          "Failed to converge after $(itres.iterations).")
end

NeXLCore.compare(itres::IterationResult, known::Material)::DataFrame = compare(itres.comp, known)

NeXLCore.compare(itress::AbstractVector{IterationResult}, known::Material)::DataFrame =
    mapreduce(itres -> compare(itres, known), append!, itress)

NeXLCore.material(itres::IterationResult) = itres.comp

mutable struct Counter
    count::Int
    terminate::Int
    Counter(terminate::Int) = new(0, terminate)
end

update(it::Counter)::Bool = (it.count += 1) <= it.terminate

terminated(it::Counter) = it.count > it.terminate

value(it::Counter) = it.count

"""
    computeKs(iter::Iteration, est::Material)::Dict{Element, Float64}

Given an estimate of the composition compute the corresponding k-ratios.
"""
function computeKs(iter::Iteration, est::Material, measured::Vector{KRatio}, stdZafs::Dict{KRatio,MultiZAF})::Dict{Element,Float64}
    estkrs = Dict{Element,Float64}()
    for kr in measured
        if nonnegk(kr) > 0.0
            # Build ZAF for unk
            @timeit iter.timer "ZAF[unk]" unkZaf = _ZAF(iter, est, kr.unkProps, kr.lines)
            # Compute the total correction and the resulting k-ratio
            @timeit iter.timer "gZAFc" gzafc = gZAFc(
                unkZaf,
                stdZafs[kr],
                kr.unkProps[:TakeOffAngle],
                kr.stdProps[:TakeOffAngle],
            )
            estkrs[kr.element] = max(1.0e-7, gzafc * est[kr.element] / kr.standard[kr.element])
        else
            estkrs[kr.element] = 0.0
        end
    end
    return estkrs
end

"""
    iterateks(iter::Iteration, name::String, measured::Vector{KRatio})

Iterate to find the composition that produces the measured k-ratios.
"""
function iterateks(iter::Iteration, name::String, measured::Vector{KRatio})::IterationResult
    eval(computed) = sum((nonnegk(kr) - computed[kr.element])^2 for kr in measured)
    stdZafs = Dict( ( kr, _ZAF(iter, kr.standard, kr.stdProps, kr.lines) ) #
                        for kr in filter(k->k.kratio > 0.0, measured) )
    @timeit iter.timer "FirstComp" estcomp = firstEstimate(iter, name, measured)
    @timeit iter.timer "FirstKs" estkrs = computeKs(iter, estcomp, measured, stdZafs)
    # If no convergence report it but return closest result...
    bestComp, bestKrs, bestEval = estcomp, estkrs, eval(estkrs)
    iters = Counter(100)
    reset(iter.updater)
    norm = 1.0
    while !converged(iter.converged, measured, estkrs) && update(iters)
        @timeit iter.timer "NextEst" upd = update(iter.updater, estcomp, measured, estkrs)
        @timeit iter.timer "Unmeasured" unmeas = compute(iter.unmeasured, upd)
        estcomp = asnormalized(material(name, unmeas), norm)
        #estcomp = material(name, unmeas)
        @timeit iter.timer "ComputeKs" estkrs = computeKs(iter, estcomp, measured, stdZafs)
        if (iters.count > 4) && (iters.count % 5 == 0)
            norm = 1.0 / mapreduce(mkr->unmeas[mkr.element], +, measured)
        end
        if eval(estkrs)<bestEval
            bestComp, bestKrs, bestEval = estcomp, estkrs, eval(estkrs)
        end
    end
    if terminated(iters)
        @warn "$(name) did not converge - using best non-converged result."
        return IterationResult(bestComp, measured, bestKrs, false, iters.count)
    else
        return IterationResult(estcomp, measured, estkrs, true, iters.count)
    end
end
