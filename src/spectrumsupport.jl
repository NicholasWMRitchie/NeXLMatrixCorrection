using .NeXLSpectrum
"""
   quantify(ffr::FilterFitResult, strip::AbstractVector{Element} = [], mc::Type{<:MatrixCorrection} = XPP, fl::Type{<:FluorescenceCorrection} = ReedFluorescence)::IterationResult

Facilitates quantifying `FilterFitResult` objects from extracting k-ratios from measured spectra. 
"""
function quantify(
    ffr::FilterFitResult,
    strip::AbstractVector{Element} = [],
    mc::Type{<:MatrixCorrection} = XPP,
    fl::Type{<:FluorescenceCorrection} = ReedFluorescence,
)::IterationResult
    iter = Iteration(mc, fl, updater = NeXLMatrixCorrection.WegsteinUpdateRule())
    krs = filter(kr -> !(element(kr) in strip), kratios(ffr))
    skro = SimpleKRatioOptimizer(1.5)
    return iterateks(iter, ffr.label, optimizeks(skro, krs))
end
