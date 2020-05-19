using DataFrames

"""
The `UnmeasuredElementRule` mechanism provides a method to implement rules for adding unmeasured elements to
the fitting process.  Examples include element-by-stoichiometry or element-by-difference.
"""
abstract type UnmeasuredElementRule end

"""
The NullUnmeasuredRule adds no additional elements in the iteration process.
"""
struct NullUnmeasuredRule <: UnmeasuredElementRule end

"""
    compute(::Type{UnmeasuredElementRule}, inp::Dict{Element,Float64})::Dict{Element,Float64}

A null UnmeasuredElementRule.  Just returns the inputs.
"""
compute(::NullUnmeasuredRule, inp::Dict{Element,Float64})::Dict{Element,Float64} = inp

"""
The `UpdateRule` abstract type defines mechanisms to update the best composition estimate between
iteration steps.
"""
abstract type UpdateRule end

"""
The `NaiveUpdateRule` implements the 'method of successive approximations' to update
the composition between iteration steps.
"""
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
    Dict(kr.element => value(nonnegk(kr)) * value(kr.standard[kr.element]) / zafs[kr.element] for kr in measured)

function reset(::NaiveUpdateRule) end

"""
The `WegsteinUpdateRule` implements the very effective method of Reed and Mason (S.J.B. Reed and P.K. Mason,
Transactions ofthe Second National Conference on Electron Microprobe Analysis, Boston, 1967.) for updating
the composition estimate between iteration steps.
"""
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
    fn = Dict(kr.element => value(kr.standard[kr.element]) / zafs[kr.element] for kr in measured)
    if !ismissing(weg.prevc)
        cn, cnm1, fnm1 = prevcomp, weg.prevc, weg.prevfn
        for mkr in measured
            elm, km = mkr.element, value(nonnegk(mkr))
            if km > 0
                δfδc = (fn[elm] - fnm1[elm]) / (value(cn[elm]) - value(cnm1[elm]))
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
"""
    reset(weg::WegsteinUpdateRule)

Restart the WegsteinUpdateRule accumulators.
"""
function reset(weg::WegsteinUpdateRule)
    weg.prevc = missing
    weg.prevfn = missing
end

"""
The `RecordingUpdateRule` wraps other `UpdateRule` instances to provide diagnostics which may be tabulated
using asa(DataFrame, rur::RecordingUpdateRule) or plotted using Gadfly's plot(rur::RecordingUpdateRule).
"""
struct RecordingUpdateRule <: UpdateRule
    base::UpdateRule
    zafs::Vector{Dict{Element,Float64}}
    comps::Vector{Dict{Element,Float64}}
    prev::Vector{Material}
    meas::Dict{Element,KRatio}
    """
        RecordingUpdateRule(ur::UpdateRule)

    Wrap an UpdateRule instance with diagnostic recorders.
    """
    RecordingUpdateRule(ur::UpdateRule) = new(
        ur,
        Vector{Dict{Element,Float64}}(),
        Vector{Dict{Element,Float64}}(),
        Vector{Material}(),
        Dict{Element,KRatio}(),
    )
end
"""
    NeXLUncertainties.asa(::Type{DataFrame}, rur::RecordingUpdateRule)::DataFrame

Tabulate the iteration steps in a DataFrame.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, rur::RecordingUpdateRule)
    dzafs, dcs = Dict{Element,Vector{Float64}}(), Dict{Element,Vector{Float64}}()
    prev, meas = Dict{Element,Vector{Float64}}(), Dict{Element,Vector{Float64}}()
    allelms = union(keys(rur.zafs[1]), keys(rur.comps[1]), keys(rur.prev[1]), keys(rur.meas))
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

"""
The `ConvergenceTest` abstract type represents mechanisms to decide when the iteration has converged.

    converged(ct::ConvergenceTest, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool
"""
abstract type ConvergenceTest end

"""
The `RMSBelowTolerance` `ConvergenceTest` ensures that the root-mean-squared difference between
measured and computed is below a threshold.
"""
struct RMSBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(rbt::RMSBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool =
    sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in meas) < rbt.tolerance^2

"""
The `AllBelowTolerance` `ConvergenceTest` ensures that the difference between
measured and computed is below a threshold for each k-ratio.
"""
struct AllBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(abt::AllBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,Float64})::Bool =
    all(abs(value(nonnegk(kr)) - computed[kr.element]) < abt.tolerance for kr in meas)


"""
The `IsApproximate` `ConvergenceTest` checks that the k-ratio differences are either below an absolute threshold
or a relative tolerance.
"""
struct IsApproximate <: ConvergenceTest
    atol::Float64
    rtol::Float64
end

converged(ia::IsApproximate, meas::Vector{KRatio}, computed::Dict{Element,Float64}) = all(
    (abs(1.0 - value(nonnegk(kr)) / computed[kr.element]) < rtol) ||
    (abs(value(nonnegk(kr)) - computed[kr.element]) < atol) for kr in meas
)
"""
Collects the information necessary to define the iteration process including the `MatrixCorrection` and
`FLuorescenceCorrection` algorithms, the iteration `UpdateRule`, the `ConvergenceTest`, and an
`UnmeasuredElementRule`.
"""
struct Iteration
    mctype::Type{<:MatrixCorrection}
    fctype::Type{<:FluorescenceCorrection}
    cctype::Type{<:CoatingCorrection}
    updater::UpdateRule
    converged::ConvergenceTest
    unmeasured::UnmeasuredElementRule

    Iteration(
        mct::Type{<:MatrixCorrection},
        fct::Type{<:FluorescenceCorrection},
        cct::Type{<:CoatingCorrection};
        updater = WegsteinUpdateRule(),
        converged = RMSBelowTolerance(0.00001),
        unmeasured = NullUnmeasuredRule(),
    ) = new(mct, fct, cct, updater, converged, unmeasured)
end


"""
`IterationResult` contains the results of the iteration process including a Label identifying the source of
the k-ratios, the resulting Material, the initial and final k-ratios, whether the iteration converged and the
number of steps.  The results can be output using `asa(DataFrame, ir::IterationResult)`.
"""
struct IterationResult
    label::Label
    comp::Material
    kratios::Vector{KRatio}
    computed::Dict{Element,Float64}
    converged::Bool
    iterations::Int
    iterate::Iteration
end

"""
The source of the k-ratio data as a Label (often a CharXRayLabel).
"""
source(ir::IterationResult)::Label = ir.label

function NeXLUncertainties.asa(::Type{DataFrame}, ir::IterationResult; withZAF::Bool = true)
    elms, mfs, ks, cks, labels = String[], Float64[], Union{Missing,Float64}[], Union{Missing,Float64}[], Label[]
    g, z, a, f = Union{Missing,Float64}[], Union{Missing,Float64}[], Union{Missing,Float64}[], Union{Missing,Float64}[]
    c, gzafc, stds, dmfs = Union{Missing,Float64}[], Union{Missing,Float64}[], String[], Float64[]
    for elm in keys(ir.comp)
        push!(labels, ir.label)
        push!(elms, elm.symbol)
        rc = round(ir.comp[elm])
        push!(mfs, value(rc))
        push!(dmfs, σ(rc))
        added = false
        for kr in ir.kratios
            if elm == kr.element
                push!(ks, value(kr.kratio))
                if withZAF
                    zafs = _ZAF(ir.iterate, kr.standard, kr.stdProps, kr.lines)
                    zafu = _ZAF(ir.iterate, ir.comp, kr.unkProps, kr.lines)
                    push!(cks, σ(kr.kratio))
                    push!(stds, name(kr.standard))
                    push!(c, coating(zafu, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle]))
                    push!(z, Z(zafu, zafs))
                    push!(a, A(zafu, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle]))
                    push!(f, F(zafu, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle]))
                    push!(g, generation(zafu, zafs))
                    push!(gzafc, g[end] * z[end] * a[end] * f[end])
                else
                    push!(cks, ir.computed[elm])
                end
                added = true
                break
            end
        end
        if !added
            push!(ks, missing), push!(cks, missing), push!(stds, missing), push!(c, missing), push!(z, missing),
            push!(a, missing), push!(f, missing), push!(g, missing), push!(gzafc, missing)
        end
    end
    return withZAF ?
           DataFrame(
        :Label => labels,
        :Element => elms,
        :Standard => stds,
        #:Converged => [ir.converged for elm in keys(ir.comp)],
        #:Iterations => [ir.iterations for elm in keys(ir.comp)],
        Symbol("Mass Frac.") => mfs,
        Symbol("Δ[Mass Frac.]") => dmfs,
        Symbol("k[Meas]") => ks,
        Symbol("Δk[Meas]") => cks,
        :Generation => g,
        :Z => z,
        :A => a,
        :F => f,
        :Coating => c,
        :gZAFc => gzafc,
    ) :
           DataFrame(
        :Label => labels,
        :Element => elms,
        :Converged => [ir.converged for elm in keys(ir.comp)],
        :Iterations => [ir.iterations for elm in keys(ir.comp)],
        Symbol("Mass Frac.") => mfs,
        Symbol("Δ[Mass Frac.]") => dmfs,
        Symbol("k[Meas]") => ks,
        Symbol("k[Comp]") => cks,
    )
    return res
end

NeXLUncertainties.asa(::Type{DataFrame}, irs::AbstractVector{IterationResult}; withZAF::Bool = true)::DataFrame =
    mapreduce(ir->asa(DataFrame,ir,withZAF=withZAF), vcat, irs)

DataFrames.describe(irs::AbstractVector{IterationResult}) =
    describe(asa(DataFrame, collect(map(ir->ir.comp, irs)))[:, 2:end], :mean, :std, :min, :q25, :median, :q75, :max)

Base.show(io::IO, itres::IterationResult) = print(
    io,
    itres.converged ? "Converged to $(itres.comp) in $(itres.iterations) steps\n" :
        "Failed to converge after $(itres.iterations): Best estimate = $(itres.comp).",
)

NeXLCore.compare(itres::IterationResult, known::Material)::DataFrame = compare(itres.comp, known)

NeXLCore.compare(itress::AbstractVector{IterationResult}, known::Material)::DataFrame =
    mapreduce(itres -> compare(itres, known), append!, itress)
"""
    NeXLCore.material(itres::IterationResult)::Material
"""
NeXLCore.material(itres::IterationResult) = itres.comp

_ZAF(iter::Iteration, mat::Material, props::Dict{Symbol,Any}, lines::Vector{CharXRay})::MultiZAF =
    zafcorrection(iter.mctype, iter.fctype, iter.cctype, mat, lines, props[:BeamEnergy], get(props, :Coating, missing))

"""
    computeZAFs(
        iter::Iteration,
        est::Material,
        stdZafs::Dict{KRatio,MultiZAF}
    )::Dict{Element, Float64}

Given an estimate of the composition compute the corresponding k-ratios.
"""
function computeZAFs(iter::Iteration, est::Material, stdZafs::Dict{KRatio,MultiZAF})::Dict{Element,Float64}
    zaf(kr, zafs) =
        gZAFc(_ZAF(iter, est, kr.unkProps, kr.lines), zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle])
    return Dict(kr.element => zaf(kr, zafs) for (kr, zafs) in stdZafs)
end
"""
    quantify(iter::Iteration, label::Label, measured::Vector{KRatio}, maxIter::Int = 100)::IterationResult
    quantify(iter::Iteration, name::String, measured::Vector{KRatio})::IterationResult
    quantify(ffr::FitResult; strip::AbstractVector{Element} = Element[], mc::Type{<:MatrixCorrection} = XPP, fl::Type{<:FluorescenceCorrection} = ReedFluorescence, cc::Type{<:CoatingCorrection} = Coating)::IterationResult

Perform the iteration procedurer as described in `iter` using the `measured` k-ratios to produce the best
estimate `Material` in an `IterationResult` object.  The third form makes it easier to quantify the
k-ratios from filter fit spectra.
"""
quantify(iter::Iteration, name::String, measured::Vector{KRatio}) = quantify(iter, label(name), measured)

function quantify(iter::Iteration, label::Label, measured::Vector{KRatio}, maxIter::Int = 100)::IterationResult
    # Compute the C = k*C_std estimate
    firstEstimate(meas::Vector{KRatio})::Dict{Element,Float64} =
        Dict(kr.element => value(nonnegk(kr)) * value(kr.standard[kr.element]) for kr in meas)
    # Compute the estimated k-ratios
    computeKs(estComp::Material, zafs::Dict{Element,Float64}, stdComps::Dict{Element,Float64})::Dict{Element,Float64} =
        Dict(elm => estComp[elm] * zafs[elm] / stdComps[elm] for (elm, zaf) in zafs)
    function computefinal(estcomp, measured)
        final = Dict{Element,UncertainValue}(elm=>convert(UncertainValue, estcomp[elm]) for elm in keys(estcomp))
        for kr in measured
            elm = element(kr)
            final[elm] = uv(value(final[elm]), (value(final[elm])/value(kr.kratio))*σ(kr.kratio))
        end
        return material(name(estcomp), final)
    end
    # Compute the k-ratio difference metric
    eval(computed) = sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in measured)
    # Compute the standard matrix correction factors
    stdZafs = Dict(kr => _ZAF(iter, kr.standard, kr.stdProps, kr.lines) for kr in measured)
    stdComps = Dict(kr.element => value(kr.standard[kr.element]) for kr in measured)
    # First estimate c_unk = k*c_std
    estcomp = material(repr(label), compute(iter.unmeasured, firstEstimate(measured)))
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
                fc = computefinal(estcomp, measured)
                return IterationResult(label, fc, measured, bestKrs, true, bestIter, iter)
            end
        end
        # Compute the next estimated mass fractions
        upd = update(iter.updater, estcomp, measured, zafs)
        # Apply unmeasured element rules
        estcomp = material(repr(label), compute(iter.unmeasured, upd))
        # calculated matrix correction for estcomp
        zafs = computeZAFs(iter, estcomp, stdZafs)
    end
    @warn "$label did not converge in $(maxIter)."
    @warn "   Using best non-converged result from step $(bestIter)."
    return IterationResult(label, bestComp, measured, bestKrs, false, bestIter, iter)
end

quantify(
    sampleName::String,
    measured::Vector{KRatio};
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
    cc::Type{<:CoatingCorrection} = Coating,
) = quantify(Iteration(mc, fc, cc), label(sampleName), measured)

quantify(
    lbl::Label,
    measured::Vector{KRatio};
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
    cc::Type{<:CoatingCorrection} = Coating,
) = quantify(Iteration(mc, fc, cc), lbl, measured)
