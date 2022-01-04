# NeXLMatrixCorrection

The matrix correction package in the **NeXL** microanalysis library for Julia.  `NeXLMatrixCorrection` depends upon
`NeXLUncertainties` and `NeXLCore`.

```julia
> using Pkg;
> Pkg.add("NeXLMatrixCorrection")
```
or 

```julia
> ]add NeXLMatrixCorrection
```

Currently `NeXLMatrixCorrection` implements the XPP, Riveros, CITZAF and XPhi matrix correction and Reed fluorescence correction algorithms for bulk and coated samples.  The library is designed to make it easy to add additional algorithms.

Primarily the algorithms in `NeXLMatrixCorrection` are designed to take a `Vector{NeXLCore.KRatio}` and return a `NeXLCore.Material`.  Since they are intended for both WDS and EDS, the k-ratio can represent one or more characteristic X-ray lines from a single element.  K-ratios compare a measured intensity with the intensity from a reference (standard) material. Typically, these two materials are measured at the same beam energy but multiple beam energy measurements are also supported.

The primary method is
```julia
function quantify(
    iter::Iteration,  # The Iteration object providing algorithmic details
    label::Label,     # A label for the unknown
    measured::Vector{KRatio}, # A complete set of k-ratios (one per element)
    maxIter::Int = 100, # The maximum number of iterations to try before
)::IterationResult
```
which is wrapped as
```julia
quantify(iter::Iteration, sampleName::String, measured::Vector{KRatio})
```
and
```julia
quantify(ffr::FilterFitResult,...)
```
to simplify usage.


The `KRatio` structure is defined in `NeXLCore`.
```Julia
KRatio(
    lines::Vector{CharXRay},
    unkProps::Dict{Symbol,<:Any},
    stdProps::Dict{Symbol,<:Any},
    standard::Material,
    kratio::AbstractFloat,
)
```
Usually it is sufficient to define the `unkProps` and `stdProps` corresponding to the `:BeamEnergy`, the `:TakeOffAngle`
which are, of course, in eV and radians.

```julia
KRatio(
  characteristic(n"F",ktransitions), # Builds a vector containing all the K shell characteristic x-rays for F
  # [ n"F K-L3" ], # An alternative with only one transition
  Dict(:BeamEnergy=>15.0e3, :TakeOffAngle=>deg2rad(40.0)),
  Dict(:BeamEnergy=>15.0, :TakeOffAngle=>deg2rad(40.0)),
  mat"CaF2",  # The standard material
  0.324 # The k-ratio
)
```

# Structure and Method Documentation

```@autodocs
Modules = [NeXLMatrixCorrection]
```
