using Test
using NeXLMatrixCorrection
using Random

function testIterate(unk, stds, e0, θ)
    toa = deg2rad(θ)
    props = Dict{Symbol,Any}(:BeamEnergy => e0, :TakeOffAngle => toa)
    krs = KRatio[]
    for (lines, std) in stds
        elm = element(lines[1])
        zu = zafcorrection(XPP, ReedFluorescence, Coating, unk, lines, e0)
        zs = zafcorrection(XPP, ReedFluorescence, Coating, std, lines, e0)
        k = gZAFc(zu, zs, toa, toa) * unk[elm] / std[elm]
        push!(krs, KRatio(lines, props, props, std, k))
    end
    up = RecordingUpdateRule(NeXLMatrixCorrection.WegsteinUpdateRule())
    iter = Iteration(XPP, ReedFluorescence, Coating, updater = up)
    return quantify(iter, "Result", krs)
end

randomize(mat::Material, qty::Float64)::Material =
    material(mat.name, (elm => mat[elm] * (1.0 + qty * randn()) for elm in keys(mat))...)

function evaluate(mat, stds, e0, θ, tol = 0.0001)
    res = testIterate(mat, stds, e0, θ)
    return res.converged && all(abs(res.comp[elm] - mat[elm]) < tol for elm in keys(mat))
end


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
end
