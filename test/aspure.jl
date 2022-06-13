using NeXLMatrixCorrection
using Test

@testset "aspure" begin
    krp = aspure(KRatio(
        [n"O K-L3"],
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        mat"SiO2",
        uv(0.605, 0.01),
    ))
    @test isapprox(value(krp), 0.13009, atol=0.0001)

    krp = aspure(KRatio(
        [n"Si K-L3"],
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        mat"SiO2",
        uv(0.139188, 0.01),
    ))
    @test isapprox(value(krp), 0.0529, atol=0.0001)
    @test isapprox(σ(krp), 0.0038, atol=0.0001)

    krp = aspure(KRatio(
        [n"Ca K-L3"],
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        mat"CaF2",
        uv(0.193347, 0.001),
    ))
    @test isapprox(value(krp), 0.0962, atol=0.0001)
    @test isapprox(σ(krp), 0.0005, atol=0.0001)

    krps = aspure(KRatios(
        [n"O K-L3"],
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        Dict(:BeamEnergy => 20.0e3, :TakeOffAngle => deg2rad(40.0)),
        mat"SiO2",
        [ 0.605 , 0.547978 ],
    ))
    @test isapprox(value(krps[1]), 0.13009, atol = 0.0001) 
    @test isapprox(value(krps[2]), 0.1179, atol = 0.0001) 
end