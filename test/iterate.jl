using Test
using NeXLMatrixCorrection

@testset "Iteration" begin
    # Measurement conditions
    props = Dict{Symbol,Any}(:BeamEnergy=>15.0e3,:TakeOffAngle=>deg2rad(40.0))
    # Unknown material
    unk = material("0.6Fe+0.4Cr",Dict(n"Fe"=>0.6, n"Cr"=>0.4))
    # Standard materials
    fe = pure(n"Fe")
    cr = pure(n"Cr")

    # Characteristic x-rays
    fek = characteristic(n"Fe",kalpha,0.01)
    crk = characteristic(n"Cr",kalpha,0.01)

    toa, e0 = props[:TakeOffAngle], props[:BeamEnergy]
    # Computed k-ratios
    kfek = NeXLMatrixCorrection.k(ZAF(XPPCorrection, ReedFluorescence, unk, fe, fek, e0)...,toa,toa)
    kcrk =  NeXLMatrixCorrection.k(ZAF(XPPCorrection, ReedFluorescence, unk, cr, crk, e0)...,toa,toa)

    # Use these as the "measured" k-ratios
    krs = [
        KRatio(fek,props,props,fe,kfek),
        KRatio(crk,props,props,cr,kcrk),
    ]
    print(krs)
    ENV["Columns"]=160
    up = RecordingUpdateRule(NeXLMatrixCorrection.WegsteinUpdateRule())
    iter=Iteration(XPPCorrection,ReedFluorescence, updater=up)
    println(analyticaltotal(unk))
    res=iterateks(iter, "Result", krs)
    print(res)
    @test res.converged
end
