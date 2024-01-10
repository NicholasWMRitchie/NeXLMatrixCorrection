using NeXLMatrixCorrection
using Test

@testset "XPP" begin
    @testset "K240" begin # Test XPP against DTSA-II using K240

        toa = deg2rad(40.0)
        k240 = NeXLCore.material(
            "K240",
            n"O" => 0.340023,
            n"Mg" => 0.030154,
            n"Si" => 0.186986,
            n"Ti" => 0.059950,
            n"Zn" => 0.040168,
            n"Zr" => 0.074030,
            n"Ba" => 0.268689,
        )

        cxr = n"O K-L3"
        xpp = NeXLMatrixCorrection.XPP(k240, inner(cxr), 20.0e3)
        @test isapprox(xpp.A, 390.68, atol = 0.01)
        @test isapprox(xpp.a, 6754, atol = 1.0)
        @test isapprox(xpp.B, -377242.5, atol = 1.0)
        @test isapprox(xpp.b, 7788, atol = 1.0)
        @test isapprox(xpp.ϕ0, 1.63011, atol = 0.0001)
        @test isapprox(xpp.F, 0.0016734, atol = 0.000001)

        xppO = NeXLMatrixCorrection.XPP(pure(element(cxr)), inner(cxr), 20.0e3)
        ZA = ℱχ(xpp, cxr, toa) / ℱχ(xppO, cxr, toa)
        Z = ℱ(xpp) / ℱ(xppO)
        @test isapprox(Z, 1.1759, atol = 0.0001)
        @test isapprox(ZA / Z, 0.2950, atol = 0.0001)

        cxr = n"Zr L3-M5"
        xpp = NeXLMatrixCorrection.XPP(k240, inner(cxr), 20.0e3)
        @test isapprox(xpp.A, 129.23, atol = 0.01)
        @test isapprox(xpp.a, 6668, atol = 1.0)
        @test isapprox(xpp.B, -176883.3, atol = 1.0)
        @test isapprox(xpp.b, 8239, atol = 1.0)
        @test isapprox(xpp.ϕ0, 1.60571, atol = 0.0001)
        @test isapprox(xpp.F, 0.0012887, atol = 0.00001)

        xppO = NeXLMatrixCorrection.XPP(pure(element(cxr)), inner(cxr), 20.0e3)
        ZA = ℱχ(xpp, cxr, toa) / ℱχ(xppO, cxr, toa)
        Z = ℱ(xpp) / ℱ(xppO)
        @test isapprox(Z, 0.8502, atol = 0.0001)
        @test isapprox(ZA / Z, 0.7463, atol = 0.0001)

        cxr = n"Mg K-L3"
        xpp = NeXLMatrixCorrection.XPP(k240, inner(cxr), 20.0e3)
        @test isapprox(xpp.A, 215.55, atol = 0.01)
        @test isapprox(xpp.a, 6681, atol = 1.0)
        @test isapprox(xpp.B, -252271.2, atol = 1.0)
        @test isapprox(xpp.b, 7973, atol = 1.0)
        @test isapprox(xpp.ϕ0, 1.620603, atol = 0.0001)
        @test isapprox(xpp.F, 0.001468275, atol = 0.00001)

        xppO = NeXLMatrixCorrection.XPP(pure(element(cxr)), inner(cxr), 20.0e3)
        ZA = ℱχ(xpp, cxr, toa) / ℱχ(xppO, cxr, toa)
        Z = ℱ(xpp) / ℱ(xppO)
        @test isapprox(Z, 1.0890, atol = 0.0001)
        @test isapprox(ZA / Z, 0.38915, atol = 0.0001)
    end

    @testset "K240 ZAF" begin

        k240 = NeXLCore.material(
            "K240",
            n"O" => 0.340023,
            n"Mg" => 0.030154,
            n"Si" => 0.186986,
            n"Ti" => 0.059950,
            n"Zn" => 0.040168,
            n"Zr" => 0.074030,
            n"Ba" => 0.268689,
        )

        sio2 = atomicfraction("Quartz", n"Si" => 1, n"O" => 2)
        mgo = atomicfraction("MgO", n"Mg" => 1, n"O" => 1)
        baf2 = atomicfraction("Barium Fluoride", n"Ba" => 1, n"F" => 2)
        ti, zn, zr = pure(n"Ti"), pure(n"Zn"), pure(n"Zr")

        e0, θ = 17.0e3, deg2rad(40.0)

        zafSi = zafcorrection(XPP, ReedFluorescence, Coating, k240, sio2, n"Si K", e0)
        zafMg = zafcorrection(XPP, ReedFluorescence, Coating, k240, mgo, n"Mg K", e0)
        zafBa = zafcorrection(XPP, ReedFluorescence, Coating, k240, baf2, n"Ba L3", e0)
        zafTi = zafcorrection(XPP, ReedFluorescence, Coating, k240, ti, n"Ti K", e0)
        zafZn = zafcorrection(XPP, ReedFluorescence, Coating, k240, zn, n"Zn K", e0)
        zafZr = zafcorrection(XPP, ReedFluorescence, Coating, k240, zr, n"Zr L3", e0)
        zafO = zafcorrection(XPP, ReedFluorescence, Coating, k240, sio2, n"O K", e0)

        @test isapprox(Z(zafSi...), 1.1345, atol = 0.001) # 1.1411
        @test isapprox(Z(zafMg...), 1.1280, atol = 0.001) # 1.1279
        @test isapprox(Z(zafBa...), 0.8259, atol = 0.001) # 0.8261
        @test isapprox(Z(zafTi...), 0.9446, atol = 0.001) # 0.9439
        @test isapprox(Z(zafZn...), 0.8976, atol = 0.001) # 0.8978
        @test isapprox(Z(zafZr...), 0.8442, atol = 0.001) # 0.8443
        @test isapprox(Z(zafO...), 1.1316, atol = 0.001) # 1.1382

        @test isapprox(A(zafSi..., n"Si K-L3", θ, θ), 0.7652, atol = 0.001) # 0.7407
        @test isapprox(A(zafMg..., n"Mg K-L3", θ, θ), 0.5996, atol = 0.001) # 0.6021
        @test isapprox(A(zafBa..., n"Ba L3-M5", θ, θ), 1.0124, atol = 0.001) # 1.0128
        @test isapprox(A(zafTi..., n"Ti K-L3", θ, θ), 0.9607, atol = 0.001) # 0.9567
        @test isapprox(A(zafZn..., n"Zn K-L3", θ, θ), 0.9842, atol = 0.001) # 0.9837
        @test isapprox(A(zafZr..., n"Zr L3-M5", θ, θ), 0.7929, atol = 0.001) # 0.6826
        @test isapprox(A(zafO..., n"O K-L3", θ, θ), 0.7750, atol = 0.001) # 0.8310

        @test isapprox(ZA(zafSi..., n"Si K-L3", θ, θ), 1.1345 * 0.7652, atol = 0.001)
        @test isapprox(ZA(zafMg..., n"Mg K-L3", θ, θ), 1.1280 * 0.5996, atol = 0.001)
        @test isapprox(ZA(zafBa..., n"Ba L3-M5", θ, θ), 0.8259 * 1.0124, atol = 0.001)
        @test isapprox(ZA(zafTi..., n"Ti K-L3", θ, θ), 0.9446 * 0.9607, atol = 0.001)
        @test isapprox(ZA(zafZn..., n"Zn K-L3", θ, θ), 0.8973 * 0.9842, atol = 0.001)
        @test isapprox(ZA(zafZr..., n"Zr L3-M5", θ, θ), 0.8442 * 0.7929, atol = 0.001)
        @test isapprox(ZA(zafO..., n"O K-L3", θ, θ), 1.1316 * 0.7750, atol = 0.001)

        @test isapprox(F(zafSi..., n"Si K-L3", θ, θ), 1.0015, atol = 0.001) # 1.0029
        @test isapprox(F(zafMg..., n"Mg K-L3", θ, θ), 1.0041, atol = 0.001) # 1.0038
        @test isapprox(F(zafBa..., n"Ba L3-M5", θ, θ), 0.9987, atol = 0.001) # 0.9997
        @test isapprox(F(zafTi..., n"Ti K-L3", θ, θ), 1.0042, atol = 0.001) # 1.0007
        @test isapprox(F(zafZn..., n"Zn K-L3", θ, θ), 1.000, atol = 0.001) # 1.0000
        @test isapprox(F(zafZr..., n"Zr L3-M5", θ, θ), 1.0020, atol = 0.001) # 1.0106
        @test isapprox(F(zafO..., n"O K-L3", θ, θ), 0.9996, atol = 0.001) # 0.9997
        # TryZAF results
        # K240	C	    Z	    A	    F	    k		        C`	    Z`	    A`	    F`	    k`
        # O	    0.3400	1.1829	0.3373	1.0003	0.1357	SiO2	0.6348	1.1382	0.8310	0.9997	0.6002
        # Mg	0.0302	1.0927	0.4579	1.0038	0.0152	MgO	    0.0501	1.1279	0.6021	1.0038	0.0341
        # Si	0.1870	1.0806	0.6562	1.0029	0.1330	SiO2	0.3984	1.1411	0.7407	1.0029	0.3377
        # Ti	0.0600	0.9439	0.9567	1.0007	0.0542	Ti					
        # Zn	0.4020	0.8978	0.9837	1.0000	0.3550	Zn					
        # Zr	0.0740	0.8443	0.6826	1.0106	0.0431	Zr					
        # Ba	0.2687	0.7457	1.0361	1.0038	0.2084	BaF2	0.3430	0.8261	1.0128	1.0038	0.2881
        # SiO2											
        # Si	0.4694	0.9470	0.8859	1.0000	0.3938						
        # O	    0.5356	1.0393	0.4059	1.0006	0.2261						
        # MgO											
        # Mg	0.6030	0.9688	0.7605	1.0000	0.4443						
        # O	    0.3970	1.0500	0.5357	1.0018	0.2237						
        # BaF2											
        # Ba	0.7833	0.9027	1.0230	1.0000	0.7233						
        # F	    0.2167	1.3080	0.6009	1.0001	0.1703						
        #    print(asa(DataFrame, Dict( [ zafSi, zafMg, zafBa, zafTi, zafZn, zafZr, zafO ]), θ))
    end

    @testset "U3O8 at 25 keV" begin

        e0, toa = 25.0e3, deg2rad(40.0)
        u3o8 = atomicfraction("U3O8", n"U" => 3, n"O" => 8)
        k227 = material("K227", n"O" => 0.1639, n"Si" => 0.0935, n"Pb" => 0.7427)
        u = pure(n"U")

        zafO = zafcorrection(XPP, ReedFluorescence, Coating, u3o8, k227, n"O K", e0)
        zafU_L = zafcorrection(XPP, ReedFluorescence, Coating, u3o8, u, n"U L3", e0)
        zafU_M = zafcorrection(XPP, ReedFluorescence, Coating, u3o8, u, n"U M5", e0)

        @test isapprox(Z(zafO...), 1.0735, atol = 0.001)
        @test isapprox(Z(zafU_L...), 0.8965, atol = 0.001)
        @test isapprox(Z(zafU_M...), 0.9197, atol = 0.001)

        @test isapprox(A(zafO..., n"O K-L3", toa, toa), 1.5847, atol = 0.001)
        @test isapprox(A(zafU_L..., n"U L3-M5", toa, toa), 1.0089, atol = 0.001)
        @test isapprox(A(zafU_M..., n"U M5-N7", toa, toa), 1.0533, atol = 0.001)

        @test isapprox(F(zafO..., n"O K-L3", toa, toa), 0.9999, atol = 0.001)
        @test isapprox(F(zafU_L..., n"U L3-M5", toa, toa), 0.9999, atol = 0.001)
        @test isapprox(F(zafU_M..., n"U M5-N7", toa, toa), 0.9999, atol = 0.001)

    end


    @testset "K240 multiZAF" begin

        k240 = NeXLCore.material(
            "K240",
            n"O" => 0.340023,
            n"Mg" => 0.030154,
            n"Si" => 0.186986,
            n"Ti" => 0.059950,
            n"Zn" => 0.040168,
            n"Zr" => 0.074030,
            n"Ba" => 0.268689,
        )

        sio2 = atomicfraction("Quartz", n"Si" => 1, n"O" => 2)
        mgo = atomicfraction("MgO", n"Mg" => 1, n"O" => 1)
        baf2 = atomicfraction("Barium Fluoride", n"Ba" => 1, n"F" => 2)
        ti, zn, zr = pure(n"Ti"), pure(n"Zn"), pure(n"Zr")

        e0, θ = 17.0e3, deg2rad(40.0)

        cxrSi = characteristic(n"Si", ktransitions, 0.001, e0)
        zafSi = zafcorrection(XPP, ReedFluorescence, Coating, k240, sio2, cxrSi, e0)
        cxrMg = characteristic(n"Mg", ktransitions, 0.001, e0)
        zafMg = zafcorrection(XPP, ReedFluorescence, Coating, k240, mgo, cxrMg, e0)
        cxrBa = characteristic(n"Ba", ltransitions, 0.001, e0)
        zafBa = zafcorrection(XPP, ReedFluorescence, Coating, k240, baf2, cxrBa, e0)
        cxrTi = characteristic(n"Ti", ktransitions, 0.001, e0)
        zafTi = zafcorrection(XPP, ReedFluorescence, Coating, k240, ti, cxrTi, e0)
        cxrZn = characteristic(n"Zn", kalpha, 0.001, e0)
        zafZn = zafcorrection(XPP, ReedFluorescence, Coating, k240, zn, cxrZn, e0)
        cxrZr = characteristic(n"Zr", ltransitions, 0.001, e0)
        zafZr = zafcorrection(XPP, ReedFluorescence, Coating, k240, zr, cxrZr, e0)
        cxrO = characteristic(n"O", ktransitions, 0.001, e0)
        zafO = zafcorrection(XPP, ReedFluorescence, Coating, k240, sio2, cxrO, e0)

        @test isapprox(Z(zafSi...), 1.1345, atol = 0.001)
        @test isapprox(Z(zafMg...), 1.1280, atol = 0.001)
        @test_broken isapprox(Z(zafBa...), 0.8330, atol = 0.001)
        @test isapprox(Z(zafTi...), 0.9446, atol = 0.001)
        @test isapprox(Z(zafZn...), 0.8973, atol = 0.001)
        @test isapprox(Z(zafZr...), 0.8442, atol = 0.001)
        @test isapprox(Z(zafO...), 1.1316, atol = 0.001)

        @test isapprox(A(zafMg..., θ, θ), 0.5996, atol = 0.001)
        @test_broken isapprox(A(zafBa..., θ, θ), 1.0172, atol = 0.001)
        @test isapprox(A(zafTi..., θ, θ), 0.9615, atol = 0.001)
        @test isapprox(A(zafZn..., θ, θ), 0.9841, atol = 0.001)
        @test_broken isapprox(A(zafZr..., θ, θ), 0.7596, atol = 0.001)
        @test isapprox(A(zafO..., θ, θ), 0.7750, atol = 0.001)

        @test isapprox(F(zafSi..., θ, θ), 1.0015, atol = 0.001) # 1.0026
        @test isapprox(F(zafMg..., θ, θ), 1.0041, atol = 0.001) # 1.0033
        @test isapprox(F(zafBa..., θ, θ), 0.9999, atol = 0.001) # 1.0043
        @test isapprox(F(zafTi..., θ, θ), 1.0042, atol = 0.002) # 1.0009
        @test isapprox(F(zafZn..., θ, θ), 1.000, atol = 0.001)  # 1.000
        @test isapprox(F(zafZr..., θ, θ), 1.0015, atol = 0.001) # 1.0106
        @test isapprox(F(zafO..., θ, θ), 0.9996, atol = 0.001)  # 1.0006
        #    print(asa(DataFrame, Dict( [ zafSi, zafMg, zafBa, zafTi, zafZn, zafZr, zafO ]), θ))
    end

    @testset "Quantify K458" begin
        # Material	K458 = [O(0.3193 mass frac),Si(0.2308 mass frac),Zn(0.0307 mass frac),Ba(0.4192 mass frac),Σ=1.0000,3.5 g/cc]
        # Detector	Bruker 5 eV/ch - FWHM[Mn Kα]=128.6 eV - 2019-02-27 00:00
        # Algorithm	XPP - Pouchou & Pichoir Simplified
        # MAC	NIST-Chantler 2005
        # E0	15 keV
        # Take-off	40°

        # IUPAC	Seigbahn	Standard	Energy	 ZAF	  Z	  A	  F	k-ratio
        # O K-L3	O Kα1	SiO2	0.5249	1.2446	1.1591	1.0743	0.9996	0.746227
        # Si K-L3	Si Kα1	SiO2	1.7397	0.8937	1.1635	0.7672	1.0012	0.441263
        # Si K-M3	Si Kβ1	SiO2	1.8290	0.9150	1.1635	0.7854	1.0012	0.451796
        # Zn K-L3	Zn Kα1	Pure Zn	8.6389	0.9048	0.9201	0.9833	1.0000	0.027776
        # Zn K-M3	Zn Kβ1	Pure Zn	9.5720	0.9083	0.9201	0.9872	1.0000	0.027884
        # Zn L3-M5	Zn Lα1	Pure Zn	1.0116	0.4506	0.9261	0.4860	1.0011	0.013833
        # Ba L3-M5	Ba Lα1	BaCl	4.4663	0.8490	0.8175	1.0384	1.0001	0.447794
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
        @test res.converged
        @test isapprox(value(res.comp[n"Si"]), 0.2308, atol=0.0004)
        @test isapprox(value(res.comp[n"Ba"]), 0.4192, atol=0.0021)
        @test isapprox(value(res.comp[n"O"]), 0.3193, atol=0.0001)
        @test isapprox(value(res.comp[n"Zn"]), 0.0307, atol=0.0001)
    end


end
