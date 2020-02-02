using DataFrames
using TimerOutputs

struct Optimize
    mctype::Type{<:MatrixCorrection}
    fctype::Type{<:FluorescenceCorrection}
    unmeasured::UnmeasuredElementRule
    timer::TimerOutput

    Optimize(
        mct::Type{<:MatrixCorrection},
        fct::Type{<:FluorescenceCorrection};
        unmeasured = NullUnmeasuredRule(),
    ) = new(mct, fct, unmeasured, TimerOutput())
end

function _ZAF(iter::Iteration, mat::Material, props::Dict{Symbol,Any}, lines::Vector{CharXRay})::MultiZAF
    coating = get(props, :Coating, NullCoating())
    return ZAF(iter.mctype, iter.fctype, mat, lines, props[:BeamEnergy], coating)
end


"""
    firstEstimate(iter::Iteration)::Material

Make a first estimate at the composition.
"""
function firstEstimate(optim::Optimize, name::String, measured::Vector{KRatio})::Material
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
function computeKs(optim::Optimize, est::Material, measured::Vector{KRatio}, stdZafs::Dict{KRatio,MultiZAF})::Dict{Element,Float64}
    estkrs = Dict{Element,Float64}()
    for kr in measured
        if nonnegk(kr) > 0.0
            # Build ZAF for unk
            unkZaf = _ZAF(iter, est, kr.unkProps, kr.lines)
            # Compute the total correction and the resulting k-ratio
            gzafc = gZAFc(
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

function iterateks(optim::Optimize, name::String, measured::Vector{KRatio})::IterationResult
    est = firstEstimate(optim, name, measured)
    elms = sort!([ keys(est)... ])
    function f(c::Vector{Float64})
        est = toMat(c)
        ks = computeKs(optim, est, measured, stdZafs)
    end
    toVec(mat::Material) = Float64[ mat[elm] for elm in elms ]
    toMat(vec::Vector{Float64}) = material(name,Dict{Element,Float64}(elm=>c for (elm, c) in zip(elms,vec)))
    res = Optim.minimizer(optimize(f, toVec(est), ConjugateGradient(); autodiff = :forward)
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
