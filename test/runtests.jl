using Test
using NeXLMatrixCorrection

@testset "NeXLMatrixCorrection" begin
    include("xpp.jl")
    include("citzaf.jl")
    include("riveros.jl")
    include("iterate.jl")
    include("xppu.jl")
    include("aspure.jl")
    include("defaultstandards.jl")
end