using Documenter
using Weave
using NeXLMatrixCorrection

include("../weave/buildweave.jl")

makedocs(modules = [NeXLMatrixCorrection], sitename = "NeXLMatrixCorrection.jl")

deploydocs(repo = "github.com/NeXLMatrixCorrection.jl.git")
