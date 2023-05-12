# ![](docs/src/NeXL_sm.png)MatrixCorrection
| **Documentation**                        | **Build Status**                  |
|:----------------------------------------:|:---------------------------------:|
| [![][docs-stable-img]][docs-stable-url]  | [![][travis-img]][travis-url]     |


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://pages.nist.gov/NeXLMatrixCorrection.jl
[travis-img]: https://travis-ci.com/usnistgov/NeXLMatrixCorrection.jl.svg?branch=master
[travis-url]: https://travis-ci.com/usnistgov/NeXLMatrixCorrection.jl


#### Installation
Install NeXLMatrixCorrection using the Julia package manager

    > ]add NeXLMatrixCorrection

or

    > using Pkg
    > Pkg.add("NeXLMatrixCorrection")


The matrix correction package in the NeXL microanalysis library for Julia.  `NeXLMatrixCorrection` depends upon
`NeXLUncertainties` and `NeXLCore`. 

#### Quick Doc
Currently `NeXLMatrixCorrection` implements the XPP, CitZAF, XPhi and Riveros matrix correction and Reed fluorescence correction algorithms for bulk and coated samples.  The library is designed to make it easy to add additional algorithms.

Primarily the algorithms in `NeXLMatrixCorrection` are designed to take a `Vector{NeXLCore.KRatio}` and return a `NeXLCore.Material`.  Using broadcast syntax, it is possible to apply these algorithms directly to hyperspectrum data sets.  However, there are versions optimized for hyper-spectral data sets, that take a `Vector{KRatios}` and return a `NeXLCore.Materials`.   When used with (`NeXLSpectrum`)[https://github.com/usnistgov/NeXLSpectrum.jl], `quantify(...)` is extended to chain the `fit_spectrum(...)` operation to both fit standards to extract k-ratios and matrix correct both `Spectrum` and `HyperSpectrum` objects.

The k-ratio is a comparison of a measured intensity with the intensity from a reference (standard) material measured under the same conditions.  Since they are intended for both WDS and EDS, the k-ratio can represent one or more characteristic X-ray lines from a single element.  You can mix and match EDS and WDS k-ratios within a single matrix correction.  You can also mix k-ratios measured at different beam energies or, hypothetically, even different instruments.

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
lbl = label("K458")
unkProps = Dict(:BeamEnergy=>15.0e3, :TakeOffAngle=>deg2rad(40.0))
stdProps = unkProps # Same for both (in this case...)
krs = [
    KRatio([n"O K-L3"], unkProps, stdProps, mat"SiO2", 0.746227 ),
    KRatio([n"Si K-L3"], unkProps, stdProps, mat"SiO2", 0.441263 ),
    KRatio([n"Zn K-L3"], unkProps, stdProps, mat"Zn", 0.027776 ),
    KRatio([n"Ba L3-M5"], unkProps, stdProps, mat"BaCl", 0.447794 )
]
res = quantify(lbl, krs)
# Converged to K458[Si=0.2311,Ba=0.4212,O=0.3192,Zn=0.0307] in 7 steps
```

### An example using NeXLSpectrum
```julia
using NeXLSpectrum
path=joinpath(datadir(), "exp_raw", "ADM6005a spectra")
refs=references( [
        reference(n"C", joinpath(path, "C std.msa"), mat"C"),
        reference(n"O", joinpath(path, "SiO2 std.msa"), mat"SiO2"),
        reference(n"Si", joinpath(path, "SiO2 std.msa"), mat"SiO2"),
        reference(n"Al", joinpath(path, "Al std.msa"), mat"Al"),
        reference(n"Ca", joinpath(path, "CaF2 std.msa"), mat"CaF2"),
        reference(n"Ti", joinpath(path, "Ti trimmed.msa"), mat"Ti"),
        reference(n"Zn", joinpath(path, "Zn std.msa"), mat"Zn"),
        reference(n"Ge", joinpath(path, "Ge std.msa"), mat"Ge"),
    ], 132.0
)
spec = sp = loadspectrum(joinpath(datadir(), "exp_raw", "ADM6005a spectra", "ADM-6005a_1.msa"))
quantify(spec, refs, strip = [ n"C" ], coating = n"C K-L2"=>parse(Material, "C", density=1.9))
```