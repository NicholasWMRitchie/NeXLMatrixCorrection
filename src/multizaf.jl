# Implements ZAF matrix correction based on a k-ratio from multiple simultaneously
# measured characteristic X-ray lines.  MultiZAF can hold the ZAF correction
using BoteSalvatICX

struct MultiZAF
    xrays::Vector{CharXRay}
    zafs::Dict{AtomicShell,ZAFCorrection}
    function MultiZAF(xrays, zafs)
        elm = element(xrays[1])
        mat = material(first(values(zafs)))
        e0 = beamEnergy(first(values(zafs)))
        @assert(
            all(f -> isequal(element(f), elm), xrays),
            "MultiZAF constructor: All the characteristic X-rays must be from the same element.",
        )
        @assert(
            all(f -> isequal(element(f), elm), keys(zafs)),
            "MultiZAF constructor: All the shells must be from the same element as the X-rays.",
        )
        @assert(
            all(f -> haskey(zafs, inner(f)), xrays),
            "MultiZAF constructor: There must be a ZAF correction for each characteristic X-ray.",
        )
        @assert(
            all(f -> isequal(material(f), mat), values(zafs)),
            "MultiZAF constructor: All the materials must match.",
        )
        @assert(
            all(f -> isequal(beamEnergy(f), e0), values(zafs)),
            "MultiZAF constructor: All the beam energies must match.",
        )
        return new(xrays, zafs)
    end
end

"""
    shells(mz::MultiZAF)
A set of all shells supported by this MultiZAF
"""
shells(mz::MultiZAF) = keys(zafs)

NeXLCore.element(mz::MultiZAF) = element(xray[1])

NeXLCore.characteristic(mz::MultiZAF) = mz.xrays

NeXLCore.material(mz::MultiZAF) = material(first(values(mz.zafs)))

beamEnergy(mz::MultiZAF) = beamEnergy(first(values(mz.zafs)))

NeXLCore.name(mz::MultiZAF) =
    repr(brightest(mz.xrays)) * "+" * string(length(mz.xrays) - 1) * " others"

commonXrays(mz1::MultiZAF, mz2::MultiZAF) =
    union(characteristic(unk), characteristic(std))

function Z(unk::MultiZAF, std::MultiZAF)
    z, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        zz = Z(zafU, zafS)
        for cxr in cxrs2
            z += weight(cxr) * zz
            n += weight(cxr)
        end
    end
    return z / n
end

function A(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)
    a, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            a += weight(cxr) * A(zafU, zafS, cxr, θtoa)
            n += weight(cxr)
        end
    end
    return a / n
end

function F(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)
    f, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            f += weight(cxr) * F(zafU, zafS, cxr, θtoa)
            n += weight(cxr)
        end
    end
    return f / n
end

function generation(unk::MultiZAF, std::MultiZAF)
    g, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        icxU, icxS = ionizationCrossSection(
                sh,
                beamEnergy(zafU),
            ),
            ionizationCrossSection(sh, beamEnergy(zafS))
        for cxr in cxrs2
            g += weight(cxr) * (icxU / icxS)
            n += weight(cxr)
        end
    end
    return g / n
end

function gZAFc(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)
    a, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        icxU, icxS = ionizationCrossSection(
                sh,
                beamEnergy(zafU),
            ),
            ionizationCrossSection(sh, beamEnergy(zafS))
        for cxr in cxrs2
            a += weight(cxr) * (icxU / icxS) * ZAFc(zafU, zafS, cxr, θtoa)
            n += weight(cxr)
        end
    end
    return a / n
end

function coating(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)
    c, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            c += weight(cxr) * coating(zafU, zafS, cxr, θtoa)
            n += weight(cxr)
        end
    end
    return c / n
end

"""
    summarize(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)::DataFrame

Summarize a matrix correction relative to the specified unknown and standard in
a DataFrame.
"""
function NeXLCore.summarize(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)::DataFrame
    tot = gZAFc(unk, std)
    return DataFrame(
        Unknown = [name(material(unk))],
        E₀ᵤ = [beamEnergy(unk)],
        Standard = [name(material(std))],
        E₀ₛ = [beamEnergy(std)],
        Xrays = [name(union(characteristic(unk), characteristic(std)))],
        Generation = [generation(unk, std)],
        Z = [Z(unk, std)],
        A = [A(unk, std, θtoa)],
        F = [F(unk, std, θtoa)],
        coating = [coating(unk, std, θtoa)],
        gZAFc = [tot],
        k = [tot * material(unk)[elm] / material(std)[elm]],
    )
end

"""
    detail(unk::MultiZAF, std::MultiZAF)::DataFrame

Tabulate each term in the MultiZAF matrix correction in a DataFrame.
"""
function detail(unk::MultiZAF, std::MultiZAF, θtoa::AbstractFloat)::DataFrame
    stds, stdE0, unks = Vector{String}(), Vector{Float64}(), Vector{String}()
    unkE0, xray, g = Vector{Float64}(), Vector{CharXRay}(), Vector{Float64}()
    z, a, f = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    c, wgt, zaf, k = Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}()
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        matU, matS = material(zafU), material(zafS)
        elm = element(sh)
        for cxr in cxrs2
            push!(unks, name(material(zafU)))
            push!(unkE0, beamEnergy(zafU))
            push!(stds, name(material(zafS)))
            push!(stdE0, beamEnergy(zafS))
            push!(xray, cxr)
            push!(wgt, weight(cxr))
            push!(g, generation(zafU, zafS, cxr))
            push!(z, Z(zafU, zafS))
            push!(a, A(zafU, zafS, cxr, θtoa))
            push!(f, F(zafU, zafS, cxr, θtoa))
            push!(c, coating(zafU, zafS, cxr))
            tot = ZAFc(zafU, zafS, cxr, θtoa)
            push!(zaf, tot)
            push!(k, tot * matU[elm] / matS[elm])
        end
    end
    return DataFrame(
        Unknown = unks,
        E0unk = 0.001*unkE0,
        Standard = stds,
        E0std = 0.001*stdE0,
        Xray = xray,
        Weight = w,
        Generation = g,
        Z = z,
        A = a,
        F = f,
        c = c,
        ZAF = zaf,
        k = k,
    )
end

"""
    detail(mzs::AbstractArray{Tuple{MultiZAF, MultiZAF}})::DataFrame

Summarize a matrix correction relative to the specified unknown and standard in
a DataFrame.
"""
detail(mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θtoa::AbstractFloat) =
    mapreduce((unk, std) -> detail(unk, std, θtoa), append!, mzs)

"""
    summarize(mzs::Dict{MultiZAF, MultiZAF}}, θtoa::AbstractFloat)::DataFrame

Summarize a matrix correction relative to a specified Dict of unknowns and
standards in a DataFrame.
"""
NeXLCore.summarize(mzs::Dict{MultiZAF,MultiZAF}, θtoa::AbstractFloat) =
    mapreduce((unk, std) -> summarize(unk, std, θtoa), append!, mzs)
