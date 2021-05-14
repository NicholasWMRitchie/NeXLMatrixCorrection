using Test
using NeXLMatrixCorrection
using Random

function testIterate(unk, stds, e0, θ)
    toa = deg2rad(θ)
    props = Dict{Symbol,Any}(:BeamEnergy => e0, :TakeOffAngle => toa)
    krs = KRatio[]
    for (xrays, std) in stds
        elm = element(xrays[1])
        zu = zafcorrection(XPP, ReedFluorescence, Coating, unk, xrays, e0)
        zs = zafcorrection(XPP, ReedFluorescence, Coating, std, xrays, e0)
        k = gZAFc(zu, zs, toa, toa) * unk[elm] / std[elm]
        push!(krs, KRatio(xrays, props, props, std, k))
    end
    up = RecordingUpdateRule(NeXLMatrixCorrection.WegsteinUpdateRule())
    iter = Iteration(XPP, ReedFluorescence, Coating, updater = up)
    return quantify(iter, "Result", krs)
end

randomize(mat::Material, qty::Float64)::Material =
    material(mat.name, (elm => mat[elm] * (1.0 + qty * randn()) for elm in keys(mat))...)

function evaluate(mat, stds, e0, θ, tol = 0.0001)
    res = testIterate(mat, stds, e0, θ)
    return res.converged && all(abs(value(res.comp[elm]) - value(mat[elm])) < tol for elm in keys(mat))
end

Base.eps(::Type{UncertainValue}) = eps(Float64)
Base.isapprox(uv1::UncertainValue, uv2::UncertainValue; atol::Real=0.0) =
    isapprox(value(uv1),value(uv2),atol=atol)
Base.isapprox(uv1::UncertainValue, v2::Real; atol::Real=0.0) =
    isapprox(value(uv1), v2, atol=atol)
Base.isapprox(v1::Real, uv2::UncertainValue; atol::Real=0.0) =
    isapprox(v1, value(uv2), atol=atol)

@testset "Iteration tests" begin
    @testset "FeCr series at 15 keV (normalized inputs)" begin
        mat = material("0.6Fe+0.4Cr", n"Fe" => 0.6, n"Cr" => 0.4)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.2Fe+0.8Cr", n"Fe" => 0.2, n"Cr" => 0.8)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.8Fe+0.2Cr", n"Fe" => 0.8, n"Cr" => 0.2)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.99Fe+0.01Cr", n"Fe" => 0.99, n"Cr" => 0.01)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.01Fe+0.99Cr", n"Fe" => 0.01, n"Cr" => 0.99)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00055)
    end

    @testset "FeCr series at 15 keV (non-normalized inputs)" begin
        mat = material("0.62Fe+0.4Cr", n"Fe" => 0.62, n"Cr" => 0.4)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.18Fe+0.84Cr", n"Fe" => 0.18, n"Cr" => 0.84)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.85Fe+0.22Cr", n"Fe" => 0.85, n"Cr" => 0.22)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.95Fe+0.01Cr", n"Fe" => 0.95, n"Cr" => 0.01)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)

        mat = material("0.03Fe+0.99Cr", n"Fe" => 0.03, n"Cr" => 0.99)
        res = testIterate(mat, Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr")), 15.0e3, 40.0)

        @test res.converged
        @test isapprox(res.comp[n"Fe"], mat[n"Fe"], atol = 0.00005)
        @test isapprox(res.comp[n"Cr"], mat[n"Cr"], atol = 0.00005)
    end

    @testset "K240 tests - Pure standards" begin
        mat = material(
            "K240",
            n"O" => 0.340023,
            n"Mg" => 0.030154,
            n"Si" => 0.186986,
            n"Ti" => 0.059950,
            n"Zn" => 0.040168,
            n"Zr" => 0.074030,
            n"Ba" => 0.268689,
        )
        stds = Dict(
            [n"O K-L3"] => pure(n"O"),
            [n"Si K-L3"] => pure(n"Si"),
            [n"Mg K-L3"] => pure(n"Mg"),
            [n"Ba L3-M5"] => pure(n"Ba"),
            [n"Ti K-L3"] => pure(n"Ti"),
            [n"Zn K-L3"] => pure(n"Zn"),
            [n"Zr L3-M5"] => pure(n"Zr"),
        )
        res = testIterate(mat, stds, 20.0e3, 40)
        @test res.converged
        @test isapprox(res.comp[n"O"], mat[n"O"], atol = 0.00005)
        @test isapprox(res.comp[n"Si"], mat[n"Si"], atol = 0.00005)
        @test isapprox(res.comp[n"Mg"], mat[n"Mg"], atol = 0.00005)
        @test isapprox(res.comp[n"Ba"], mat[n"Ba"], atol = 0.00005)
        @test isapprox(res.comp[n"Ti"], mat[n"Ti"], atol = 0.00005)
        @test isapprox(res.comp[n"Zn"], mat[n"Zn"], atol = 0.00005)
        @test isapprox(res.comp[n"Zr"], mat[n"Zr"], atol = 0.00005)
    end

    @testset "K240 tests - Simple standards" begin
        mat = material(
            "K240",
            n"O" => 0.340023,
            n"Mg" => 0.030154,
            n"Si" => 0.186986,
            n"Ti" => 0.059950,
            n"Zn" => 0.040168,
            n"Zr" => 0.074030,
            n"Ba" => 0.268689,
        )
        stds = Dict(
            [n"O K-L3"] => parse(Material, "SiO2"),
            [n"Si K-L3"] => parse(Material, "SiO2"),
            [n"Mg K-L3"] => parse(Material, "MgO"),
            [n"Ba L3-M5"] => parse(Material, "BaF2"),
            [n"Ti K-L3"] => pure(n"Ti"),
            [n"Zn K-L3"] => pure(n"Zn"),
            [n"Zr L3-M5"] => pure(n"Zr"),
        )
        res = testIterate(mat, stds, 20.0e3, 40)
        @test res.converged
        @test isapprox(res.comp[n"O"], mat[n"O"], atol = 0.00005)
        @test isapprox(res.comp[n"Si"], mat[n"Si"], atol = 0.00005)
        @test isapprox(res.comp[n"Mg"], mat[n"Mg"], atol = 0.00005)
        @test isapprox(res.comp[n"Ba"], mat[n"Ba"], atol = 0.00005)
        @test isapprox(res.comp[n"Ti"], mat[n"Ti"], atol = 0.00005)
        @test isapprox(res.comp[n"Zn"], mat[n"Zn"], atol = 0.00005)
        @test isapprox(res.comp[n"Zr"], mat[n"Zr"], atol = 0.00005)
    end
    @testset "K240 tests - Simple standards - EDS mode" begin
        mat = randomize(
            material(
                "K240",
                n"O" => 0.340023,
                n"Mg" => 0.030154,
                n"Si" => 0.186986,
                n"Ti" => 0.059950,
                n"Zn" => 0.040168,
                n"Zr" => 0.074030,
                n"Ba" => 0.268689,
            ),
            0.1,
        )
        stds = Dict(
            characteristic(n"O", ktransitions) => parse(Material, "SiO2"),
            characteristic(n"Si", ktransitions) => parse(Material, "SiO2"),
            characteristic(n"Mg", ktransitions) => parse(Material, "MgO"),
            characteristic(n"Ba", ltransitions) => parse(Material, "BaF2"),
            characteristic(n"Ti", ktransitions) => pure(n"Ti"),
            characteristic(n"Zn", ktransitions) => pure(n"Zn"),
            characteristic(n"Zr", kalpha) => pure(n"Zr"),
        )
        res = testIterate(mat, stds, 20.0e3, 40)
        @test res.converged
        @test isapprox(res.comp[n"O"], mat[n"O"], atol = 0.00005)
        @test isapprox(res.comp[n"Si"], mat[n"Si"], atol = 0.00005)
        @test isapprox(res.comp[n"Mg"], mat[n"Mg"], atol = 0.00005)
        @test isapprox(res.comp[n"Ba"], mat[n"Ba"], atol = 0.00005)
        @test isapprox(res.comp[n"Ti"], mat[n"Ti"], atol = 0.00005)
        @test isapprox(res.comp[n"Zn"], mat[n"Zn"], atol = 0.00005)
        @test isapprox(res.comp[n"Zr"], mat[n"Zr"], atol = 0.00005)
    end
    @testset "K240 tests - Unnorm" begin
        mat = material(
            "K240",
            n"O" => 0.340023,
            n"Mg" => 0.030154,
            n"Si" => 0.186986,
            n"Ti" => 0.059950,
            n"Zn" => 0.040168,
            n"Zr" => 0.074030,
            n"Ba" => 0.268689,
        )
        stds = Dict(
            [n"O K-L3"] => parse(Material, "SiO2"),
            [n"Si K-L3"] => parse(Material, "SiO2"),
            [n"Mg K-L3"] => parse(Material, "MgO"),
            [n"Ba L3-M5"] => parse(Material, "BaF2"),
            [n"Ti K-L3"] => pure(n"Ti"),
            [n"Zn K-L3"] => pure(n"Zn"),
            [n"Zr L3-M5"] => pure(n"Zr"),
        )
        Random.seed!(0xEA7BADF00D0)
        @test evaluate(randomize(mat, 0.01), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.1), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.2), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.3), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.4), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
    end
    @testset "FeCr tests - Unnorm" begin
        mat = material("FeCr", n"Fe" => 0.6, n"Cr" => 0.4)
        stds = Dict([n"Fe K-L3"] => pure(n"Fe"), [n"Cr K-L3"] => pure(n"Cr"))
        Random.seed!(0xEA7BADF00D)
        @test evaluate(randomize(mat, 0.01), stds, 15.0e3, 40)
        @test evaluate(randomize(mat, 0.1), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.2), stds, 17.0e3, 40)
        @test evaluate(randomize(mat, 0.3), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.4), stds, 14.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 16.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 22.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 19.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 17.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 20.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 14.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 26.0e3, 40)
        @test evaluate(randomize(mat, 0.5), stds, 30.0e3, 40)
    end
    @testset "K309 test" begin
        # From a measurement
        props = Dict(:BeamEnergy=>15.0e3, :TakeOffAngle=>deg2rad(40.0))
        fe = KRatio(characteristic(n"Fe",kalpha),props,props,mat"Fe",UncertainValue(0.08884060796139037,0.000363897))
        al = KRatio(characteristic(n"Al",ktransitions),props,props,mat"Al2O3",UncertainValue(0.12212761167159845,0.000257044))
        si = KRatio(characteristic(n"Si",ktransitions),props,props,mat"Si",UncertainValue(0.13931565459523768,0.000175283))
        ba = KRatio(characteristic(n"Ba",ltransitions),props,props,mat"BaF2",UncertainValue(0.14026640286207112,0.000494559))
        o = KRatio(characteristic(n"O",ktransitions),props,props,mat"Al2O3",UncertainValue(0.636480599484999,0.000823425))
        ca = KRatio(characteristic(n"Ca",ktransitions),props,props,mat"CaF2",UncertainValue(0.20597043348028954,0.000417253))
        res = quantify(nl"K309", [al,ba,ca,fe,o,si])
        @test isapprox(res.comp[n"Fe"],0.1038,atol=0.0001) # 0.1036 in DTSA-II
        @test isapprox(res.comp[n"Al"],0.0756,atol=0.0001) # 0.0756
        @test isapprox(res.comp[n"Si"],0.1829,atol=0.0001) # 0.1823
        @test isapprox(res.comp[n"Ba"],0.1410,atol=0.0001) # 0.1402
        @test isapprox(res.comp[n"O"],0.3885,atol=0.0001)  # 0.3899
        @test isapprox(res.comp[n"Ca"],0.1060,atol=0.0001) # 0.1036
        @test isapprox(analyticaltotal(res.comp),0.997797,atol=0.000001) # 0.9986
    end
    @testset "20 nm C on SiO2 at 5 keV" begin
        props = Dict(:BeamEnergy=>5.0e3, :TakeOffAngle=>deg2rad(40.0))
        si = KRatio([n"Si K-L3"],props,props,mat"SiO2",0.9219)
        o = KRatio([n"O K-L3"],props,props,mat"SiO2",0.8950)
        c = KRatio([n"C K-L2"],props,props,mat"C",0.0548)
        res = quantify(nl"20 nm C on SiO2 at 5 keV", [si, o, c], coating=n"C K-L2"=>parse(Material,"C",density=1.9))
        @test isapprox(thickness(si.unkProps[:Coating]), 2.41e-6, atol=0.01e-6)
    end
    @testset "4 nm of Al2O3 on Al at 3 keV" begin
        props = Dict(:BeamEnergy=>3.0e3, :TakeOffAngle=>deg2rad(40.0))
        al = KRatio([n"Al K-L3"], props, props, mat"Al", 4402. / 4657.)
        o = KRatio([n"O K-L3"], props, props, mat"SiO2", 629.7 / 10061.)
        res = quantify(nl"4 nm Al2O3 on Al at 3 keV", [al, o], coating=n"O K-L3"=>parse(Material,"Al2O3",density=3.987))
        @test isapprox(thickness(al.unkProps[:Coating]), 4.85e-7, atol=0.01e-6)
    end
end
