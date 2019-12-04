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
    WegsteinUpdateRule(f::Float64=0.3) = new(Vector{Material}(), Vector{Dict{Element,Float64}}(),f)
end

function update( #
    weg::WegsteinUpdateRule,
    prevcomp::Material,
    measured::Vector{KRatio},
    estkrs::Dict{Element,Float64},
)::Dict{Element,Float64}
    bound(x,min,max) = x < min ? min : (x > max ? max : x)
    cnp1, itercx = Dict{Element,Float64}(), length(weg.prevc)
    if itercx < 1
        for mkr in measured
            cnp1[mkr.element] = estkrs[mkr.element] > 0.0 ? (nonnegk(mkr) / estkrs[mkr.element]) * prevcomp[mkr.element] : 0.0
        end
    else
        cn = itercx <= 5 ? NeXLCore.asnormalized(prevcomp) : prevcomp
        cnm1, kn, knm1 = weg.prevc[end], estkrs, weg.prevk[end]
        for mkr in measured
            elm, km = mkr.element, nonnegk(mkr)
            if (kn[elm] > 0) && (knm1[elm] > 0)
                fcn, fcnm1 = cn[elm]/kn[elm], cnm1[elm]/knm1[elm] # c = k*f
                dfdk = ((fcn - fcnm1) / (cn[elm] - cnm1[elm]))
                cnp1[elm] = cn[elm] + (km * fcn - cn[elm]) / (1.0 - bound(km*dfdk, -weg.factor, weg.factor))
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
    RecordingUpdateRule(ur::UpdateRule) = new(ur,Vector{Dict{Element,Float64}}(),Vector{Dict{Element,Float64}}())
end

function NeXLUncertainties.asa(::Type{DataFrame}, rur::RecordingUpdateRule)
    dkrs, dcs = Dict{Element, Vector{Float64}}(), Dict{Element,Vector{Float64}}()
    allelms = union(keys(rur.estkrs[1]), keys(rur.comps[1]))
    for elm in allelms
        dkrs[elm], dcs[elm] = [], []
    end
    for i in eachindex(rur.estkrs)
        for elm in allelms
            push!(dkrs[elm], get(rur.estkrs[i], elm, 0.0))
            push!(dcs[elm], get(rur.comps[i], elm, 0.0))
        end
    end
    df = DataFrame(Iter=collect(eachindex(rur.estkrs)))
    for elm in allelms
        df[!, Symbol("C($(elm.symbol))")]=dcs[elm]
        df[!, Symbol("k($(elm.symbol))")]=dkrs[elm]
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
    push!(rur.estkrs, estkrs)
    push!(rur.comps, res)
    return res
end

function NeXLMatrixCorrection.reset(rur::RecordingUpdateRule)
    resize!(rur.estkrs,0)
    resize!(rur.comps,0)
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
    stdZafs = Dict( ( kr, _ZAF(iter, kr.standard, kr.stdProps, kr.lines) ) #
                        for kr in filter(k->k.kratio > 0.0, measured) )
    @timeit iter.timer "FirstComp" estcomp = firstEstimate(iter, name, measured)
    @timeit iter.timer "FirstKs" estkrs = computeKs(iter, estcomp, measured, stdZafs)
    iters = Counter(100)
    reset(iter.updater)
    while !converged(iter.converged, measured, estkrs) && update(iters)
        @timeit iter.timer "NextEst" upd = update(iter.updater, estcomp, measured, estkrs)
        @timeit iter.timer "Unmeasured" unmeas = compute(iter.unmeasured, upd)
        @timeit iter.timer "Material" estcomp = material(name, unmeas)
        @timeit iter.timer "ComputeKs" estkrs = computeKs(iter, estcomp, measured, stdZafs)
    end
    return IterationResult(estcomp, measured, estkrs, !terminated(iters), iters.count)
end
