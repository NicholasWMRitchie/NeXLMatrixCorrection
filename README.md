# NeXLMatrixCorrection

The matrix correction package in the NeXL microanalysis library for Julia.  `NeXLMatrixCorrection` depends upon 
`NeXLUncertainties` and `NeXLCore`.

`NeXLMatrixCorrection` and some of its dependences are not yet available in the Julia registry so it must be installed using the URL.

```julia
using Pkg; 
Pkg.add.([
  "https://github.com/usnistgov/BoteSalvatICX.jl.git",
  "https://github.com/usnistgov/FFAST.jl.git",
  "https://github.com/NicholasWMRitchie/NeXLUncertainties.jl.git",
  "https://github.com/NicholasWMRitchie/NeXLCore.jl.git",
  "https://github.com/NicholasWMRitchie/NeXLMatrixCorrection.jl.git"
])
```

Currently `NeXLMatrixCorrection` implements the XPP matrix correction and Reed fluorescence correction algorithms for bulk and coated samples.  The library is designed to make it easy to add additional algorithms.

Primarily the algorithms in `NeXLMatrixCorrection` are designed to take a `Vector{NeXLCore.KRatio}` and return a `NeXLCore.Material`.  Since they are intended for both WDS and EDS, the k-ratio can represent one or more characteristic X-ray lines from a single element.  K-ratios compare a measured intensity with the intensity from a reference (standard) material. Typically, these two materials are measured at the same beam energy but multiple beam energy measurements are also supported.
