# Implements ZAF matrix correction based on a k-ratio from multiple simultaneously
# measured characteristic X-ray lines.  MultiZAF can hold the ZAF correction

struct MultiZAF
    xrays::Vector{CharXRay}
    zafs::Dict{AtomicSubShell,ZAFCorrection}
    function MultiZAF(xrays, zafs)
        elm = element(xrays[1])
        mat = material(first(values(zafs)))
        e0 = beamEnergy(first(values(zafs)))
        @assert all(f -> isequal(element(f), elm), xrays)
            "MultiZAF constructor: All the characteristic X-rays must be from the same element."
        @assert all(f -> isequal(element(f), elm), keys(zafs))
            "MultiZAF constructor: All the shells must be from the same element as the X-rays."
        @assert all(f -> haskey(zafs, inner(f)), xrays)
            "MultiZAF constructor: There must be a ZAF correction for each characteristic X-ray.",
        @assert all(f -> isequal(material(f), mat), values(zafs))
            "MultiZAF constructor: All the materials must match."
        @assert all(f -> isequal(beamEnergy(f), e0), values(zafs))
            "MultiZAF constructor: All the beam energies must match."
        return new(xrays, zafs)
    end
end

"""
    shells(mz::MultiZAF)

A set of all sub-shells supported by this MultiZAF
"""
NeXLCore.atomicsubshells(mz::MultiZAF) = keys(mz.zafs)

NeXLCore.element(mz::MultiZAF) = element(mz.xrays[1])

NeXLCore.characteristic(mz::MultiZAF) = mz.xrays

NeXLCore.material(mz::MultiZAF) = material(first(values(mz.zafs)))

beamEnergy(mz::MultiZAF) = beamEnergy(first(values(mz.zafs)))

NeXLCore.name(mz::MultiZAF) =
    repr(brightest(mz.xrays)) * "+" * string(length(mz.xrays) - 1) * " others"

commonXrays(unk::MultiZAF, std::MultiZAF) =
    union(characteristic(unk), characteristic(std))

function Z(unk::MultiZAF, std::MultiZAF)
    z, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        norm = sum(weight.(cxrs2))
        n += norm
        z += Z(unk.zafs[sh], std.zafs[sh])*norm
    end
    return z / n
end

function A(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    a, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            w = weight(cxr)
            a += w * A(zafU, zafS, cxr, θunk, θstd)
            n += w
        end
    end
    return a / n
end

function F(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    f, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            w = weight(cxr)
            f += w * F(zafU, zafS, cxr, θunk, θstd)
            n += w
        end
    end
    return f / n
end

function generation(unk::MultiZAF, std::MultiZAF)
    g, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        icxU = ionizationcrosssection(sh, beamEnergy(zafU))
        icxS = ionizationcrosssection(sh, beamEnergy(zafS))
        for cxr in cxrs2
            w = weight(cxr)
            g += w * (icxU / icxS)
            n += w
        end
    end
    return g / n
end

function coating(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    c, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            c += weight(cxr) * coating(zafU, zafS, cxr, θunk, θstd)
            n += weight(cxr)
        end
    end
    return c / n
end

function gZAFc(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    a, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        icxU = ionizationcrosssection(sh, beamEnergy(zafU))
        icxS = ionizationcrosssection(sh, beamEnergy(zafS))
        z = Z(zafU, zafS)
        for cxr in cxrs2
            w = weight(cxr)
            a += w * (icxU / icxS) * z * A(zafU, zafS, cxr, θunk, θstd) *
                F(zafU, zafS, cxr, θunk, θstd) * coating(zafU, zafS, cxr, θunk, θstd)
            n += w
        end
    end
    return a / n
end

"""
    k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat) =

Compute the k-ratio for the specified measurement conditions. k = C / ZAF
"""
k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat) =
    (material(unk)[element(unk)] / material(std)[element(unk)])*gZAFc(unk, std, θunk, θstd)


"""
    tabulate(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard in
a DataFrame.
"""
function tabulate(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame
    tot = gZAFc(unk, std, θunk, θstd)
    @assert isequal(element(unk), element(std)) "The unknown and standard's elements must match."
    return DataFrame(
        Unknown = [name(material(unk))],
        E₀ᵤ = [beamEnergy(unk)],
        Standard = [name(material(std))],
        E₀ₛ = [beamEnergy(std)],
        Xrays = [name(union(characteristic(unk), characteristic(std)))],
        Generation = [generation(unk, std)],
        Z = [Z(unk, std)],
        A = [A(unk, std, θunk, θstd)],
        F = [F(unk, std, θunk, θstd)],
        coating = [coating(unk, std, θunk, θstd)],
        gZAFc = [tot],
        k = [tot * material(unk)[element(unk)] / material(std)[element(unk)]],
    )
end

"""
    detail(unk::MultiZAF, std::MultiZAF)::DataFrame

Tabulate each term in the MultiZAF matrix correction in a DataFrame.
"""
function detail(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame
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
            push!(g, generation(zafU, zafS, inner(cxr)))
            push!(z, Z(zafU, zafS))
            push!(a, A(zafU, zafS, cxr, θunk, θstd))
            push!(f, F(zafU, zafS, cxr, θunk, θstd))
            push!(c, coating(zafU, zafS, cxr, θunk, θstd))
            tot = ZAFc(zafU, zafS, cxr, θunk, θstd)
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
        Weight = wgt,
        Generation = g,
        Z = z,
        A = a,
        F = f,
        c = c,
        ZAF = zaf,
        k = k
    )
end

"""
    detail(mzs::AbstractArray{Tuple{MultiZAF, MultiZAF}})::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard in
a DataFrame.
"""
detail(mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat) =
    mapreduce(tmm -> detail(tmm[1], tmm[2], θunk, θstd), append!, mzs)

"""
    tabulate(mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to a specified Dict of unknowns and
standards in a DataFrame.
"""
tabulate(mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat) =
    mapreduce(tmm -> tabulate(tmm[1], tmm[2], θunk, θstd), append!, mzs)


function tabulate(unk::Material,
                std::Material,
                lines::Vector{CharXRay},
                e0::Float64,
                toa::Float64;
                zacorr::Type{<:MatrixCorrection} = XPP,
                fcorr::Type{<:FluorescenceCorrection} = ReedFluorescence)
    zafs = ZAF(zacorr, fcorr, unk, std, lines, e0)
    tabulate(zafs..., toa, toa)
end
