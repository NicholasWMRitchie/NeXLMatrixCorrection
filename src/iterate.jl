using DataFrames
using Statistics
using ThreadsX
using Statistics
using LinearAlgebra

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
    )::Dict{Element, Float64}

Determine the next estimate of the composition that brings the estimate k-ratios closer to measured.
"""
function update(
    ::NaiveUpdateRule, #
    prevcomp::Material, #
    measured::Vector{KRatio}, #
    zafs::Dict{Element,Float64},
    state::Union{Nothing, Tuple} #
)::Tuple
    cnp1 = Dict(kr.element => value(nonnegk(kr)) * value(kr.standard[kr.element]) / zafs[kr.element] for kr in measured)
    return ( cnp1, nothing )

end

"""
The `WegsteinUpdateRule` implements the very effective method of Reed and Mason (S.J.B. Reed and P.K. Mason,
Transactions ofthe Second National Conference on Electron Microprobe Analysis, Boston, 1967.) for updating
the composition estimate between iteration steps.
"""
struct WegsteinUpdateRule <: UpdateRule end

function update( #
    ::WegsteinUpdateRule,
    prevcomp::Material,
    measured::Vector{KRatio},
    zafs::Dict{Element,Float64},
    state::Union{Nothing, Tuple}
)::Tuple
    cnp1 = Dict(kr.element => value(nonnegk(kr)) * value(kr.standard[kr.element]) / zafs[kr.element] for kr in measured)
    fn = Dict(kr.element => value(kr.standard[kr.element]) / zafs[kr.element] for kr in measured)
    if !isnothing(state)
        cn = prevcomp
        cnm1, fnm1 = state
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
    return ( cnp1, ( prevcomp, fn ) )
end

"""
The `ConvergenceTest` abstract type represents mechanisms to decide when the iteration has converged.

    converged(ct::ConvergenceTest, meas::Vector{KRatio}, computed::Dict{Element,<:AbstractFloat})::Bool
"""
abstract type ConvergenceTest end

"""
The `RMSBelowTolerance` `ConvergenceTest` ensures that the root-mean-squared difference between
measured and computed is below a threshold.
"""
struct RMSBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(rbt::RMSBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,<:AbstractFloat})::Bool =
    sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in meas) < rbt.tolerance^2

"""
The `AllBelowTolerance` `ConvergenceTest` ensures that the difference between
measured and computed is below a threshold for each k-ratio.
"""
struct AllBelowTolerance <: ConvergenceTest
    tolerance::Float64
end

converged(abt::AllBelowTolerance, meas::Vector{KRatio}, computed::Dict{Element,<:AbstractFloat})::Bool =
    all(abs(value(nonnegk(kr)) - computed[kr.element]) < abt.tolerance for kr in meas)


"""
The `IsApproximate` `ConvergenceTest` checks that the k-ratio differences are either below an absolute threshold
or a relative tolerance.
"""
struct IsApproximate <: ConvergenceTest
    atol::Float64
    rtol::Float64
end

converged(ia::IsApproximate, meas::Vector{KRatio}, computed::Dict{Element,<:AbstractFloat}) = all(
    (abs(1.0 - value(nonnegk(kr)) / computed[kr.element]) < ia.rtol) ||
    (abs(value(nonnegk(kr)) - computed[kr.element]) < ia.atol) for kr in meas
)
"""
    Iteration(;
        mc::Type{<:MatrixCorrection} = XPP,
        fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
        cc::Type{<:CoatingCorrection} = Coating,
        updater = WegsteinUpdateRule(),
        converged = RMSBelowTolerance(0.00001),
        unmeasured = NullUnmeasuredRule(),
    )

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

    Iteration(;
        mc::Type{<:MatrixCorrection} = XPP,
        fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
        cc::Type{<:CoatingCorrection} = Coating,
        updater = WegsteinUpdateRule(),
        converged = RMSBelowTolerance(0.00001),
        unmeasured = NullUnmeasuredRule(),
    )  = new(mc, fc, cc, updater, converged, unmeasured)
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
    elms, mfs, ks, cks, labels = String[], Float64[], Union{Missing,Float64}[], Union{Missing,Float64}[], String[]
    g, z, a, f = Union{Missing,Float64}[], Union{Missing,Float64}[], Union{Missing,Float64}[], Union{Missing,Float64}[]
    c, gzafc, stds, dmfs, xrays = Union{Missing,Float64}[], Union{Missing,Float64}[], String[], Float64[], String[]
    for elm in keys(ir.comp)
        push!(labels, repr(ir.label))
        push!(elms, elm.symbol)
        rc = round(ir.comp[elm])
        push!(mfs, value(rc))
        push!(dmfs, σ(rc))
        added = false  # computed element
        for kr in ir.kratios
            if elm == kr.element
                push!(ks, value(kr.kratio))
                push!(xrays, repr(kr.xrays))
                if withZAF
                    zafs = _ZAF(ir.iterate, kr.standard, kr.stdProps, kr.xrays)
                    zafu = _ZAF(ir.iterate, ir.comp, kr.unkProps, kr.xrays)
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
            push!(a, missing), push!(f, missing), push!(g, missing), push!(gzafc, missing), push!(xrays, missing)
        end
    end
    return withZAF ?
           DataFrame(
        :Label => labels,
        :Element => elms,
        :Standard => stds,
        :Xrays => xrays,
        #:Converged => [ir.converged for elm in keys(ir.comp)],
        #:Iterations => [ir.iterations for elm in keys(ir.comp)],
        Symbol("C") => mfs,
        Symbol("ΔC") => dmfs,
        Symbol("k") => ks,
        Symbol("Δk") => cks,
        :g => g,
        :Z => z,
        :A => a,
        :F => f,
        :c => c,
        :gZAFc => gzafc,
    ) :
           DataFrame(
        :Label => labels,
        :Element => elms,
        :Converged => [ir.converged for elm in keys(ir.comp)],
        :Iterations => [ir.iterations for elm in keys(ir.comp)],
        Symbol("C") => mfs,
        Symbol("ΔC") => dmfs,
        Symbol("k") => ks,
        Symbol("Δk") => cks,
    )
    return res
end

function NeXLUncertainties.asa(::Type{DataFrame}, irs::AbstractVector{IterationResult}; mode = :MassFraction, nominal = nothing)::DataFrame
    if isnothing(nominal) 
        asa(DataFrame, map(ir->ir.comp,irs), mode)
    else
        asa(DataFrame, [ ( ir.comp for ir in irs)..., nominal ], mode )
    end
end


DataFrames.describe(irs::AbstractVector{IterationResult}) =
    describe(asa(DataFrame, irs)[:, 2:end], :mean, :std, :min, :q25, :median, :q75, :max)

Base.show(io::IO, itres::IterationResult) = print(
    io,
    itres.converged ? "Converged to $(itres.comp) in $(itres.iterations) steps." :
        "Failed to converge in $(itres.iterations) iterations: Best estimate = $(itres.comp).",
)

NeXLCore.compare(itres::IterationResult, known::Material)::DataFrame = compare(itres.comp, known)

NeXLCore.compare(itress::AbstractVector{IterationResult}, known::Material)::DataFrame =
    mapreduce(itres -> compare(itres, known), append!, itress)

NeXLCore.compare(iter1::IterationResult, iter2::IterationResult) = compare(iter1.comp, iter2.comp)
"""
    NeXLCore.material(itres::IterationResult)::Material
"""
NeXLCore.material(itres::IterationResult) = itres.comp
NeXLCore.material(itress::AbstractVector{IterationResult})  = mean(material.(itress))

_ZAF(iter::Iteration, mat::Material, props::Dict{Symbol,Any}, xrays::Vector{CharXRay}) =
    zafcorrection(iter.mctype, iter.fctype, iter.cctype, convert(Material{Float64,Float64}, mat), xrays, props[:BeamEnergy], get(props, :Coating, missing))

"""
    computeZAFs(
        iter::Iteration,
        mat::Material{Float64, Float64},
        stdZafs::Dict{KRatio,MultiZAF}
    )::Dict{Element, Float64}

Given an estimate of the composition compute the corresponding gZAFc matrix correction.
"""
function computeZAFs(iter::Iteration, mat::Material, stdZafs::Dict{<:NeXLCore.KRatioBase, MultiZAF})
    mat64 = convert(Material{Float64,Float64}, mat)
    function gzafc(kr, zafs)
        zaf = _ZAF(iter, mat64, kr.unkProps, kr.xrays)
        gZAFc(zaf, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle])
    end
    return Dict(kr.element => gzafc(kr, zafs) for (kr, zafs) in stdZafs)
end


"""
    estimatecoating(substrate::Material, coating::Material, kcoating::KRatio, mc::Type{<:MatrixCorrection}=XPP)::Film

Use the measured k-ratio to estimate the mass-thickness of the coating material on the specified substrate.
Return the result as a Film which may be assigned to the :Coating property for k-ratios associated with the
substrate.

Assumptions:
  * The coating is the same for all measurements of the unknown
  * The coating element k-ratio is assumed to be assigned to the brightest line
"""
function estimatecoating(substrate::Material, coating::Material, kcoating::Union{KRatio,KRatios}, mc::Type{<:MatrixCorrection}=XPP)::Film
    coatingasfilm(mc, convert(Material{Float64, Float64}, substrate), convert(Material{Float64, Float64}, coating), #
            brightest(kcoating.xrays), kcoating.unkProps[:BeamEnergy], # 
            kcoating.unkProps[:TakeOffAngle], value(kcoating.kratio))
end

"""
    quantify(
        name::Union{String, Label},
        measured::Vector{KRatio}, 
        iteration::Iteration = Iteration(mc=XPP, fc=ReedFluorescence, cc=Coating);
        maxIter::Int = 100, 
        initalEstComp::Union{Nothing,Material}=nothing, 
        coating::Union{Nothing, Pair{CharXRay, <:Material}}=nothing
    )::IterationResult

    quantify(
        measured::AbstractVector{KRatios{T}},
        iteration::Iteration = Iteration(mc=XPP, fc=NullFluorescence, cc=NullCoating);
        kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
        maxIter::Int = 100, 
        initalEstComp::Union{Nothing,Material}=nothing, 
        coating::Union{Nothing, Pair{CharXRay, Material}}=nothing
    )

Perform the iteration procedurer as described in `iter` using the `measured` k-ratios to produce the best
estimate `Material` in an `IterationResult` object.  The third form makes it easier to quantify the
k-ratios from filter fit spectra.  `initalEstComp` is an optional first estimate of the composition.  This can
be a useful optimization when quantifying many similar k-ratios (like, for example, points on a 
hyper-spectrum.)

`coating` defines a CharXRay not present in the sample that is used to estimate the thickness of the coating
layer (of the paired `Material`).
"""
quantify(
    name::Union{Label, String},
    measured::Vector{KRatio},
    iteration::Iteration = Iteration();
    maxIter::Int = 100, 
    initalEstComp::Union{Nothing,Material}=nothing, 
    coating::Union{Nothing, Pair{CharXRay, Material}}=nothing
)::IterationResult = quantify(name, measured, iteration, maxIter, initalEstComp, coating)

function quantify(
    name::Union{Label, String},
    measured::Vector{KRatio},
    iteration::Iteration,
    maxIter::Int, 
    initalEstComp::Union{Nothing,Material}, 
    coating::Union{Nothing, Pair{CharXRay, Material}}
)::IterationResult
    aslbl(lbl::Label) = lbl
    aslbl(lbl) = label(lbl)
    initalEstComp = convert(Material{Float64,Float64}, initalEstComp)
    coating = convert(Material{Float64,Float64}, coating)
    lbl = aslbl(name)
    # Compute the C = k*C_std estimate
    firstEstimate(meas::Vector{KRatio}) =
        Dict(kr.element => value(nonnegk(kr)) * value(kr.standard[kr.element]) for kr in meas)
    # Compute the estimated k-ratios
    computeKs(comp, zafs, stdComps) =
        Dict(elm => comp[elm] * zafs[elm] / stdComps[elm] for (elm, zaf) in zafs)
    function computefinal(comp, meas::Vector{KRatio})
        final = Dict(map(meas) do kr
            elm = element(kr)
            elm => if comp[elm]>0.0 && value(kr.kratio) > 0.0 && σ(kr.kratio) > 0.0 
                uv(comp[elm], comp[elm]*fractional(kr.kratio))
            else
                uv(0.0, σ(kr.kratio))
            end
        end)
        return material(NeXLCore.name(comp), final)
    end
    @assert isnothing(coating) || (get(last(coating), :Density, -1.0) > 0.0) "You must provide a positive density for the coating material."
    # Is this k-ratio due to the coating?
    iscoating(k, coatmat) =  (!isnothing(coatmat)) && (first(coatmat) in k.xrays) && (element(k) in keys(last(coatmat)))
    # k-ratios from measured elements in the unknown - Remove k-ratios for unmeasured and coating elements
    kunk = filter(measured) do kr
        !(isunmeasured(iteration.unmeasured, element(kr)) || iscoating(kr, coating))
    end
    # k-ratios associated with the coating
    kcoat = filter(kr->iscoating(kr, coating), measured)
    # Compute the k-ratio difference metric
    eval(computed) = sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in kunk)
    # Compute the standard matrix correction factors
    stdZafs = Dict(kr => _ZAF(iteration, kr.standard, kr.stdProps, kr.xrays) for kr in kunk)
    stdComps = Dict(kr.element => value(kr.standard[kr.element]) for kr in kunk)
    # First estimate c_unk = k*c_std
    currComp = something(initalEstComp, material(repr(lbl), compute(iteration.unmeasured, firstEstimate(kunk))))
    # Compute the associated matrix corrections
    zafs = computeZAFs(iteration, currComp, stdZafs)
    bestComp, bestKrs = currComp, computeKs(currComp, zafs, stdComps)
    bestEval, bestIter, iter_state = 1.0e300, 0, nothing
    for iters in Base.OneTo(maxIter)
        if length(kcoat) >= 1
            coatings = estimatecoating(currComp, last(coating), first(kcoat), iteration.mctype)
            # Previous coatings are replaced on all the unknown's k-ratios
            foreach(k->k.unkProps[:Coating] = coatings, kunk)
        end
        # How close are the calculated k-ratios to the measured version of the k-ratios?
        estkrs = computeKs(currComp, zafs, stdComps)
        if eval(estkrs) < bestEval
            # If no convergence report it but return closest result...
            bestComp, bestKrs, bestEval, bestIter = currComp, estkrs, eval(estkrs), iters
            if converged(iteration.converged, kunk, bestKrs)
                fc = computefinal(currComp, kunk)
                return IterationResult(lbl, fc, measured, bestKrs, true, bestIter, iteration)
            end
        end
        # Compute the next estimated mass fractions
        upd, iter_state = update(iteration.updater, currComp, kunk, zafs, iter_state)
        # Apply unmeasured element rules
        currComp = material(repr(lbl), compute(iteration.unmeasured, upd))
        # calculated matrix correction for currComp
        zafs = computeZAFs(iteration, currComp, stdZafs)
    end
    @warn "$lbl did not converge in $(maxIter)."
    @warn "   Using best non-converged result from step $(bestIter)."
    return IterationResult(lbl, bestComp, measured, bestKrs, false, bestIter, iteration)
end

NeXLMatrixCorrection.quantify(
    measured::Vector{KRatios},
    iteration::Iteration = Iteration(mc=XPP, fc=NullFluorescence, cc=NullCoating);
    kro::KRatioOptimizer = SimpleKRatioOptimizer(1.5),
    maxIter::Int = 100, 
    coating::Union{Nothing, Pair{CharXRay, Material}}=nothing,
    name::AbstractString = "Map",
    ty::Union{Type{UncertainValue}, Type{Float64}, Type{Float32}}=Float32,
    maxErrors::Int = 5
) = quantify(measured, iteration, kro, maxIter, coating, name, ty, maxErrors)

function NeXLMatrixCorrection.quantify(
    measured::Vector{KRatios},
    iteration::Iteration,
    kro::KRatioOptimizer,
    maxIter::Int, 
    coating::Union{Nothing, Pair{CharXRay, Material}},
    name::AbstractString,
    ty::Union{Type{UncertainValue}, Type{Float64}, Type{Float32}},
    maxErrors::Int
)
    if (!isnothing(coating)) && (typeof(coating.second) != Material{Float64, Float64})
        coating = coating.first => convert(Material{Float64, Float64}, coating.second)
    end
    @assert all(size(measured[1])==size(krsi) for krsi in measured[2:end]) "All the KRatios need to be the same dimensions."
    # Compute the C = k*C_std estimate
    firstEstimate(meas::Vector{KRatio}) =
        Dict(kr.element => value(nonnegk(kr)) * value(kr.standard[kr.element]) for kr in meas)
    # Compute the estimated k-ratios
    computeKs(comp, zafs, stdComps) =
        Dict(elm => comp[elm] * zafs[elm] / stdComps[elm] for (elm, zaf) in zafs)
    @assert isnothing(coating) || (get(last(coating), :Density, -1.0) > 0.0) "You must provide a positive density for the coating material."
    # Is this k-ratio due to the coating?
    iscoating(k, coatmat) =  (!isnothing(coatmat)) && (first(coatmat) in k.xrays) && (element(k) in keys(last(coatmat)))
    # Compute the k-ratio difference metric
    eval(computed, kunk) = sum((value(nonnegk(kr)) - computed[kr.element])^2 for kr in kunk)
    # Pick the best sub-selection of `measured` to quantify
    optmeasured = brightest.(optimizeks(kro, measured))
    # Compute the standard matrix correction factors
    stdZafs = Dict(kr => _ZAF(iteration, kr.standard, kr.stdProps, kr.xrays) for kr in optmeasured)
    stdComps = Dict{Element,Float64}(kr.element => value(kr.standard[kr.element]) for kr in optmeasured)
    mats = Materials(name, [ element(kr) for kr in optmeasured ], ty, size(measured[1]))
    nerrors = Threads.Atomic{Int}(0)
    ThreadsX.foreach(CartesianIndices(mats)) do ci
        bestComp = NeXLCore.NULL_MATERIAL
        if nerrors[] < maxErrors
            try
                measured = KRatio[ kr[ci] for kr in optmeasured ]
                # k-ratios from measured elements in the unknown - Remove k-ratios for unmeasured and coating elements
                kunk = filter(measured) do kr
                    !(isunmeasured(iteration.unmeasured, element(kr)) || iscoating(kr, coating))
                end
                # k-ratios associated with the coating
                kcoat = filter(kr->iscoating(kr, coating), measured)
                # First estimate c_unk = k*c_std
                currComp = Material("first", compute(iteration.unmeasured, firstEstimate(kunk)))
                # Compute the associated matrix corrections
                zafs = computeZAFs(iteration, currComp, stdZafs)
                bestComp, bestKrs = currComp, computeKs(currComp, zafs, stdComps)
                bestEval, bestIter, iter_state = 1.0e300, 0, nothing
                for iters in Base.OneTo(maxIter)
                    if length(kcoat) >= 1
                        coatings = estimatecoating(currComp, last(coating), first(kcoat), iteration.mctype)
                        # Previous coatings are replaced on all the unknown's k-ratios
                        foreach(k->k.unkProps[:Coating] = coatings, kunk)
                    end
                    # How close are the calculated k-ratios to the measured version of the k-ratios?
                    estKrs = computeKs(currComp, zafs, stdComps)
                    if (ev = eval(estKrs, kunk)) < bestEval
                        # If no convergence report it but return closest result...
                        bestComp, bestKrs, bestEval, bestIter = currComp, estKrs, ev, iters
                        if converged(iteration.converged, kunk, bestKrs)
                            break
                        end
                    end
                    # Compute the next estimated mass fractions
                    upd, iter_state = update(iteration.updater, currComp, kunk, zafs, iter_state)
                    # Apply unmeasured element rules
                    currComp = Material("bogus", compute(iteration.unmeasured, upd))
                    # calculated matrix correction for currComp
                    zafs = computeZAFs(iteration, currComp, stdZafs)
                end
            catch ex
                Threads.atomic_add!(nerrors, 1)
                bestComp = NeXLCore.NULL_MATERIAL
                @error ex
                show(stderr, ex)
            end
        end
        mats[ci] = bestComp
    end
    if nerrors[] > maxErrors
        @error "Exceeded $maxErrors errors - terminating early."
    end
    return mats
end