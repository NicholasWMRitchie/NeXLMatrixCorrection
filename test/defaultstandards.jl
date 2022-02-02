using Test

@testset "Default Standards" begin
    @test length(getstandards(n"Rb"))>=8 # Check that the standard parse etc.
end