# ![](docs/src/NeXL_sm.png)MatrixCorrection
| **Documentation**                        | **Build Status**                  |
|:----------------------------------------:|:---------------------------------:|
| [![][docs-stable-img]][docs-stable-url]  | [![][travis-img]][travis-url]     |


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://pages.nist.gov/NeXLMatrixCorrection.jl
[travis-img]: https://travis-ci.com/usnistgov/NeXLMatrixCorrection.jl.svg?branch=master
[travis-url]: https://travis-ci.com/usnistgov/NeXLMatrixCorrection.jl


The matrix correction package in the NeXL microanalysis library for Julia.  `NeXLMatrixCorrection` depends upon
`NeXLUncertainties` and `NeXLCore`.

`NeXLMatrixCorrection` is available through the Julia registry.  You can install it using the package manager from the command prompt.

```julia
]add NeXLMatrixCorrection
```

Currently `NeXLMatrixCorrection` implements the XPP, CitZAF, XPhi and Riveros matrix correction and Reed fluorescence correction algorithms for bulk and coated samples.  The library is designed to make it easy to add additional algorithms.

Primarily the algorithms in `NeXLMatrixCorrection` are designed to take a `Vector{NeXLCore.KRatio}` and return a `NeXLCore.Material`.  Since they are intended for both WDS and EDS, the k-ratio can represent one or more characteristic X-ray lines from a single element.  K-ratios compare a measured intensity with the intensity from a reference (standard) material. Typically, these two materials are measured at the same beam energy but multiple beam energy measurements are also supported.

The primary methods are
```julia
quantify(  # Generic method for EDS, WDS or mixed data
  sample::Union{String, Label}, # Sample name or Label
  measured::Vector{KRatio};     # The k-ratios
  mc::Type{<:MatrixCorrection}=XPP,  # Default algorithm choices
  fc::Type{<:FluorescenceCorrection}=ReedFluorescence,
  cc::Type{<:CoatingCorrection}=Coating)::IterationResult
quantify(ffr::FilterFitResult)::IterationResult  # Specialized for the results from fitted EDS spectra

# where

KRatio(
    lines::AbstractVector{CharXRay},  # CharXRay or X-rays measured
    unkProps::Dict{Symbol,<:Any},     # Properties of the measurement ( :BeamEnery, :TakeOffAngle )
    stdProps::Dict{Symbol,<:Any},     # Properties of the standard ( :BeamEnery, :TakeOffAngle )
    standard::Material,               # Composition of the standard
    kratio::AbstractFloat,            # The k-ratio (can be an UncertainValue)
)
```

### An example
```julia
julia> lbl = label("K458")
julia> unkProps = Dict(:BeamEnergy=>15.0e3, :TakeOffAngle=>deg2rad(40.0))
julia> stdProps = unkProps # Same for both (in this case...)
julia> krs = [
    KRatio([n"O K-L3"], unkProps, stdProps, mat"SiO2", 0.746227 ),
    KRatio([n"Si K-L3"], unkProps, stdProps, mat"SiO2", 0.441263 ),
    KRatio([n"Zn K-L3"], unkProps, stdProps, mat"Zn", 0.027776 ),
    KRatio([n"Ba L3-M5"], unkProps, stdProps, mat"BaCl", 0.447794 )
]
julia> res = quantify(lbl, krs)
```
Converged to K458[Si=0.2311,Ba=0.4212,O=0.3192,Zn=0.0307] in 7 steps
