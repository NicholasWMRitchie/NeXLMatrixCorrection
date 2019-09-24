using NeXLCore
using NeXLMatrixCorrection
using Test

@testset "K240" begin # Test XPP against DTSA-II using K240

    toa = deg2rad(40.0)
    k240 = NeXLCore.material("K240",Dict(n"O"=>0.340023, n"Mg"=>0.030154, n"Si"=>0.186986, n"Ti"=>0.059950, n"Zn"=>0.040168, n"Zr"=>0.074030, n"Ba"=>0.268689),missing)

    cxr = n"O K-L3"
    xpp = NeXLMatrixCorrection.XPP(k240,inner(cxr),20.0e3)
    @test isapprox(xpp.A, 393.01, atol=0.01)
    @test isapprox(xpp.a, 6754, atol=1.0)
    @test isapprox(xpp.B, -378725, atol=1.0)
    @test isapprox(xpp.b, 7786, atol=1.0)
    @test isapprox(xpp.ϕ0, 1.63011, atol = 0.0001)
    @test isapprox(xpp.F, 0.0016754, atol=0.000001)

    xppO = NeXLMatrixCorrection.XPP(pure(element(cxr)),inner(cxr),20.0e3)
    ZA = Fχ(xpp,cxr,toa)/Fχ(xppO,cxr,toa)
    Z = F(xpp)/F(xppO)
    @test isapprox(Z, 1.1759, atol = 0.0001)
    @test isapprox(ZA/Z, 0.2951, atol = 0.0001)


    cxr = n"Zr L3-M5"
    xpp = NeXLMatrixCorrection.XPP(k240,inner(cxr),20.0e3)
    @test isapprox(xpp.A, 129.53, atol=0.01)
    @test isapprox(xpp.a, 6668, atol=1.0)
    @test isapprox(xpp.B, -177201, atol=1.0)
    @test isapprox(xpp.b, 8239, atol=1.0)
    @test isapprox(xpp.ϕ0, 1.60571, atol = 0.0001)
    @test isapprox(xpp.F, 0.0012887, atol = 0.00001)

    xppO = NeXLMatrixCorrection.XPP(pure(element(cxr)),inner(cxr),20.0e3)
    ZA = Fχ(xpp,cxr,toa)/Fχ(xppO,cxr,toa)
    Z = F(xpp)/F(xppO)
    @test isapprox(Z, 0.8502, atol = 0.0001)
    @test isapprox(ZA/Z, 0.7463, atol = 0.0001)

    cxr = n"Mg K-L3"
    xpp = NeXLMatrixCorrection.XPP(k240,inner(cxr),20.0e3)
    @test isapprox(xpp.A, 215.96, atol=0.01)
    @test isapprox(xpp.a, 6681, atol=1.0)
    @test isapprox(xpp.B, -252614, atol=1.0)
    @test isapprox(xpp.b, 7973, atol=1.0)
    @test isapprox(xpp.ϕ0, 1.620603, atol = 0.0001)
    @test isapprox(xpp.F, 0.001468275, atol = 0.00001)

    xppO = NeXLMatrixCorrection.XPP(pure(element(cxr)),inner(cxr),20.0e3)
    ZA = Fχ(xpp,cxr,toa)/Fχ(xppO,cxr,toa)
    Z = F(xpp)/F(xppO)
    @test isapprox(Z, 1.0890, atol = 0.0001)
    @test isapprox(ZA/Z, 0.3889, atol = 0.0001)
end


sio2 = atomicfraction("Quartz",Dict(n"Si"=>1,n"O"=>2))
mgo = atomicfraction("MgO",Dict(n"Mg"=>1,n"O"=>1))
baf2 = atomicfraction("Barium Fluoride",Dict(n"Ba"=>1,n"F"=>2))
ti, zn, zr = pure(n"Ti"), pure(n"Zn"), pure(n"Zr")

e0, θ = 20.0e3, deg2rad(40.0)

cc, cc2 = carbonCoating(10.0), carbonCoating(20.0)

zafSi = xppZAF(k240, sio2, n"Si K", e0, stdCoating = cc2, unkCoating=cc)
zafMg = xppZAF(k240, mgo, n"Mg K", e0, unkCoating=cc)
zafBa = xppZAF(k240, baf2, n"Ba L3", e0, stdCoating = cc2, unkCoating=cc)
zafTi = xppZAF(k240, ti, n"Ti K", e0, unkCoating=cc)
zafZn = xppZAF(k240, zn, n"Zn K", e0, unkCoating=cc)
zafZr = xppZAF(k240, zr, n"Zr K", e0, unkCoating=cc)
zafO = xppZAF(k240, sio2, n"O K", e0, stdCoating = cc2, unkCoating=cc)

print(NeXLCore.summarize(Dict( [ zafSi, zafMg, zafBa, zafTi, zafZn, zafZr, zafO ]), θ))
