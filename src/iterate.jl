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
        zafs::Dict{Element, Float64}
    )::Dict{Element,Float64}

Determine the next estimate of the composition that brings the estimate k-ratios closer to measured.
"""
update(
    ::NaiveUpdateRule, #
    prevcomp::Material, #
    measured::Vector{KRatio}, #
    zafs::Dict{Element,Float64}, #
)::Dict{Element,Float64} =
    Dict(kr.element => value(nonnegk(kr)) * kr.standard[kr.element] / zafs[kr.element] for kr in measured)

function reset(::NaiveUpdateRule) end

mutable struct WegsteinUpdateRule <: UpdateRule
    fallback::UpdateRule
    prevc::Union{Material,Missing}
    prevfn::Union{Dict{Element,Float64},Missing}
    WegsteinUpdateRule() = new(NaiveUpdateRule(), missing, missing)
end

function update( #
    weg::WegsteinUpdateRule,
    prevcomp::Material,
    measured::Vector{KRatio},
    zafs::Dict{Element,Float64},
)::Dict{Element,Float64}
    # Default to naive
    cnp1 = update(weg.fallback, prevcomp, measured, zafs)
    fn = Dict(kr.element => kr.standard[kr.element] / zafs[kr.element] for kr in measured)
    if !ismissing(weg.prevc)
        cn, cnm1, fnm1 = prevcomp, weg.prevc, weg.prevfn
        for mkr in measured
            elm, km = mkr.element, value(nonnegk(mkr))
            if km > 0
                δfδc = (fn[elm] - fnm1[elm]) / (cn[elm] - cnm1[elm])
                den = 1.0 - km * δfδc
                if (abs(δfδc) < 10.0) && (abs(den) > 0.2)
                    Δc = (km * fn[elm] - cn[elm]) / den # Wegstein
                    cnp1[elm] = cn[elm] + Δc
                end
            end
        end
    end
    weg.prevc, weg.prevfn = prevcomp, fn
    return cnp1
end

function reset(weg::WegsteinUpdateRule)
    weg.prevc = missing
    weg.prevfn = missing
end

struct RecordingUpdateRule <: UpdateRule
    base::UpdateRule
    zafs::Vector{Dict{Element,Float64}}
    comps::Vector{Dict{Element,Float64}}
    prev::Vector{Material}
    meas::Dict{Element,KRatio}
    RecordingUpdateRule(ur::UpdateRule) = new(
    ur,
    Vector{Dict{Element,Float64}}(),
    Vector{Dict{Element,Float64}}(),
    Vector{Material}(),
    Dict{Element,KRatio}(),
    )
end

function NeXLUncertainties.asa(::Type{DataFrame}, rur::RecordingUpdateRule)
    dzafs, dcs = Dict{Element,Vector{Float64}}(), Dict{Element,Vector{Float64}}()
    prev, meas = Dict{Element,Vector{Float64}}(), Dict{Element,Vector{Float64}}()
    allelms = union(
        keys(rur.zafs[1]),
        keys(rur.comps[1]),
        keys(rur.prev[1]),
        keys(rur.meas),
    )
    for elm in allelms
        dzafs[elm], dcs[elm], prev[elm], meas[elm] = [], [], [], []
    end
    for i in eachindex(rur.zafs)
        for elm in allelms
            push!(dzafs[elm], get(rur.zafs[i], elm, 0.0))
            push!(dcs[elm], get(rur.comps[i], elm, 0.0))
            push!(prev[elm], rur.prev[i][elm])
            push!(meas[elm], value(rur.meas[elm].kratio))
        end
    end
    df = DataFrame(Iter = collect(eachindex(rur.zafs)))
    for elm in allelms
        df[!, Symbol("Prev($(elm.symbol))")] = prev[elm]
        df[!, Symbol("Next($(elm.symbol))")] = dcs[elm]
        df[!, Symbol("ZAF($(elm.symbol))")] = dzafs[elm]
        df[!, Symbol("meas($(elm.symbol))")] = meas[elm]
    end
    return df
end

function NeXLMatrixCorrection.update( #
    rur::RecordingUpdateRule,
    prevcomp::Material,
    measured::Vector{KRatio},
    zafs::Dict{Element,Float64},
)::Dict{Element,Float64}
    res = NeXLMatrixCorrection.update(rur.base, prevcomp, measured, zafs)
    if isempty(rur.meas)
        for kr in measured
            rur.meas[kr.element] = kr
        end
    end
    push!(rur.prev, prevcomp)
    push!(rur.zafs, zafs)
    push!(rur.comps, res)
    return res
end

function NeXLMatrixCorrection.reset(rur::RecordingUpdateRule)
    resize!(rur.zafs, 0)
    resize!(rur.comps, 0)
    resize!(rur.prev, 0)
    empty!(rur.meas)
    NeXLMatrixCorrection.reset(rur.base)
end

abstract type ConvergenceTest end

struct RMSBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(rbt::RMSBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool =
    sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in meas) < rbt.tolerance^2

struct AllBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(abt::AllBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool =
    all(abs(value(nonnegk(kr)) - computed[kr.element]) < abt.tolerance for kr in meas)

struct IsApproximate <: ConvergenceTest
    atol::Float64
    rtol::Float64
end

converged(ia::IsApproximate, meas::Vector{KRatio}, computed::Dict{Element,Float64}) =
    all((abs(1.0 - value(nonnegk(kr)) / computed[kr.element]) < rtol) || (abs(value(nonnegk(kr)) -
                                                                              computed[kr.element]) < atol) for kr in meas)

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
        converged = RMSBelowTolerance(0.00001),
        unmeasured = NullUnmeasuredRule(),
    ) = new(mct, fct, updater, converged, unmeasured, TimerOutput())
end

struct IterationResult
    comp::Material
    kratios::Vector{KRatio}
    computed::Dict{Element,Float64}
    converged::Bool
    iterations::Int
end

function NeXLUncertainties.asa(::Type{DataFrame}, ir::IterationResult)
    elms, mfs, ks, cks = Element[],
        Float64[],
        Union{Missing,Float64}[],
        Union{Missing,Float64}[]
    for elm in keys(ir.comp)
        push!(elms, elm)
        push!(mfs, ir.comp[elm])
        added = false
        for kr in ir.kratios
            if elm == kr.element
                push!(ks, value(kr.kratio))
                push!(cks, ir.computed[elm])
                added = true
                break
            end
        end
        if !added
            push!(ks, missing)
            push!(cks, missing)
        end
    end
    return DataFrame(
        :Element => elms,
        :Converged => [ir.converged for elm in keys(ir.comp)],
        :Iterations => [ir.iterations for elm in keys(ir.comp)],
        :Composition => mfs,
        :Measured => ks,
        :Computed => cks,
    )
end

function Base.show(io::IO, itres::IterationResult)
    print(
        io,
        itres.converged ? "Converged in $(itres.iterations) to $(itres.comp)\n" :
        "Failed to converge after $(itres.iterations) as $(itres.comp).",
    )
end

NeXLCore.compare(itres::IterationResult, known::Material)::DataFrame =
    compare(itres.comp, known)

NeXLCore.compare(itress::AbstractVector{IterationResult}, known::Material)::DataFrame =
    mapreduce(itres -> compare(itres, known), append!, itress)

NeXLCore.material(itres::IterationResult) = itres.comp

_ZAF(iter::Iteration, mat::Material, props::Dict{Symbol,Any}, lines::Vector{CharXRay})::MultiZAF =
    ZAF(
        iter.mctype,
        iter.fctype,
        mat,
        lines,
        props[:BeamEnergy],
        get(props, :Coating, NullCoating()),
    )

"""
    computeZAFs(
        iter::Iteration,
        est::Material,
        stdZafs::Dict{KRatio,MultiZAF}
    )::Dict{Element, Float64}

Given an estimate of the composition compute the corresponding k-ratios.
"""
function computeZAFs(iter::Iteration, est::Material, stdZafs::Dict{KRatio,MultiZAF})::Dict{
    Element,
    Float64,
}
    zaf(kr, zafs) = gZAFc(
    _ZAF(iter, est, kr.unkProps, kr.lines),
    zafs,
    kr.unkProps[:TakeOffAngle],
    kr.stdProps[:TakeOffAngle],
    )
    return Dict(kr.element => zaf(kr, zafs) for (kr, zafs) in stdZafs)
end


"""
    iterateks(iter::Iteration, name::String, measured::Vector{KRatio})

Iterate to find the composition that produces the measured k-ratios.
"""
function iterateks(
    iter::Iteration,
    name::String,
    measured::Vector{KRatio},
    maxIter::Int = 100,
)::IterationResult
    # Compute the C = k*C_std estimate
    firstEstimate(meas::Vector{KRatio}) =
        Dict(kr.element => value(nonnegk(kr)) * kr.standard[kr.element] for kr in meas)
    # Compute the estimated k-ratios
    computeKs(
        estComp::Material,
        zafs::Dict{Element,Float64},
        stdComps::Dict{Element,Float64},
    ) = Dict(elm => estComp[elm] * zafs[elm] / stdComps[elm] for (elm, zaf) in zafs)
    # Compute the k-ratio difference metric
    eval(computed) = sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in measured)
    # Compute the standard matrix correction factors
    stdZafs = Dict(kr => _ZAF(iter, kr.standard, kr.stdProps, kr.lines) for kr in measured)
    stdComps = Dict(kr.element => value(kr.standard[kr.element]) for kr in measured)
    # First estimate c_unk = k*c_std
    estcomp = material(name, compute(iter.unmeasured, firstEstimate(measured)))
    # Compute the associated matrix corrections
    zafs = computeZAFs(iter, estcomp, stdZafs)
    bestComp, bestKrs = estcomp, computeKs(estcomp, zafs, stdComps)
    bestEval, bestIter = 1.0e300, 0
    reset(iter.updater)
    for iters = 1:maxIter
        # How close are the calculated k-ratios to the measured?
        estkrs = computeKs(estcomp, zafs, stdComps)
        if eval(estkrs) < bestEval
            # If no convergence report it but return closest result...
            bestComp, bestKrs, bestEval, bestIter = estcomp, estkrs, eval(estkrs), iters
            if converged(iter.converged, measured, bestKrs)
                return IterationResult(estcomp, measured, bestKrs, true, bestIter)
            end
        end
        # Compute the next estimated mass fractions
        upd = update(iter.updater, estcomp, measured, zafs)
        # Apply unmeasured element rules
        estcomp = material(name, compute(iter.unmeasured, upd))
        # calculated matrix correction for estcomp
        zafs = computeZAFs(iter, estcomp, stdZafs)
    end
    @warn "$(name) did not converge in $(maxIter)."
    @warn "   Using best non-converged result from step $(bestIter)."
    return IterationResult(bestComp, measured, bestKrs, false, bestIter)
end
