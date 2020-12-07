"""
The `MultiZAF` structure holds the information necessary to perform matrix correction on a collection of
characteristic X-rays that were measured simultaneously from the same element.
Use `zafcorrection(...)` to construct these rather than the internal constructor.
"""
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
    zafcorrection(
      mctype::Type{<:MatrixCorrection},
      fctype::Type{<:FluorescenceCorrection},
      cctype::Type{<:CoatingCorrection},
      mat::Material,
      cxrs,
      e0,
      coating=missing
    )

Constructs a MultiZAF around the mctype and fctype algorithms for a collection of CharXRay `cxrs`.
"""
function zafcorrection(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    cctype::Type{<:CoatingCorrection},
    mat::Material,
    cxrs,
    e0::Real,
    coating::Union{Film,AbstractVector{Film},Missing} = missing,
)
    mat = asnormalized(mat)
    shells = union(filter(sh->energy(sh)<e0, inner.(cxrs)))
    zafs = Dict((sh, zafcorrection(mctype, fctype, cctype, mat, sh, e0, coating)) for sh in shells)
    return MultiZAF(cxrs, zafs)
end

"""
    zafcorrection(
      mctype::Type{<:MatrixCorrection},
      fctype::Type{<:FluorescenceCorrection},
      cctype::Type{<:CoatingCorrection},
      unk::Material,
      std::Material,
      cxrs,
      e0;
      unkCoating::Union{Film,AbstractVector{Film},Missing} = missing,
      stdCoating::Union{Film,AbstractVector{Film},Missing} = missing,
    )

Constructs a tuple of MultiZAF around the mctype and fctype correction algorithms for the unknown and standard for a
collection of CharXRay `cxrs`.
"""
function zafcorrection(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    cctype::Type{<:CoatingCorrection},
    unk::Material,
    std::Material,
    cxrs,
    e0::Real;
    unkCoating::Union{Film,AbstractVector{Film},Missing} = missing,
    stdCoating::Union{Film,AbstractVector{Film},Missing} = missing,
)
    lines = filter(cxr->energy(inner(cxr))<e0, cxrs)
    return (
        zafcorrection(mctype, fctype, cctype, unk, lines, e0, unkCoating),
        zafcorrection(mctype, fctype, cctype, std, lines, e0, stdCoating),
    )
end


"""
    shells(mz::MultiZAF)

A set of all sub-shells supported by this MultiZAF
"""
NeXLCore.atomicsubshells(mz::MultiZAF) = keys(mz.zafs)
"""
    NeXLCore.element(mz::MultiZAF)

The element associate with this `MultiZAF`
"""
NeXLCore.element(mz::MultiZAF) = element(mz.xrays[1])

"""
    NeXLCore.characteristic(mz::MultiZAF)

The X-rays associated with this `MultiZAF`.
"""
NeXLCore.characteristic(mz::MultiZAF) = mz.xrays

NeXLCore.material(mz::MultiZAF) = material(first(values(mz.zafs)))

beamEnergy(mz::MultiZAF) = beamEnergy(first(values(mz.zafs)))

NeXLCore.name(mz::MultiZAF) = repr(brightest(mz.xrays)) * "+" * string(length(mz.xrays) - 1) * " others"

commonXrays(unk::MultiZAF, std::MultiZAF) = union(characteristic(unk), characteristic(std))

"""
    Z(unk::MultiZAF, std::MultiZAF)

The Z (atomic number) correction for `unk` relative to `std`.
"""
function Z(unk::MultiZAF, std::MultiZAF)
    z, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        norm = sum(weight.(cxrs2))
        n += norm
        z += Z(unk.zafs[sh], std.zafs[sh]) * norm
    end
    return z / n
end

"""
    A(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)

The A (absorption) correction for `unk` relative to `std`.
"""
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

"""
    F(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)

The F (fluoresence) correction for `unk` relative to `std`.
"""
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

"""
    generation(unk::MultiZAF, std::MultiZAF)

The generation correction for `unk` relative to `std`.  Usually, 1.0 unless the standard and unknown were collected
at different beam energies.
"""
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

"""
    coating(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)

The conductive (or other) coating correction factor.
"""
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
"""
    gZAFc(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    gZAFc(kr::KRatio, unkComp::Material; mc::Type{<:MatrixCorrection} = XPP, fc::Type{<:FluorescenceCorrection} = ReedFluorescence, cc::Type{<:CoatingCorrection} = Coating)

The combined generation, atomic number, absorption and generation corrections.
"""
function gZAFc(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    a, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        icxU = ionizationcrosssection(sh, beamEnergy(zafU))
        icxS = ionizationcrosssection(sh, beamEnergy(zafS))
        for cxr in cxrs2
            w = weight(cxr)
            a +=
                w *
                (icxU / icxS) *
                ZA(zafU, zafS, cxr, θunk, θstd) *
                F(zafU, zafS, cxr, θunk, θstd) *
                coating(zafU, zafS, cxr, θunk, θstd)
            n += w
        end
    end
    return a / n
end

function gZAFc(
    kr::KRatio,
    unkComp::Material;
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
    cc::Type{<:CoatingCorrection} = Coating,
)
    elm = kr.element
    lines = filter(cxr->energy(inner(cxr))<min(kr.unkProps[:BeamEnergy],kr.stdProps[:BeamEnergy]), kr.lines)
    zu = zafcorrection(mc, fc, cc, unkComp, lines, kr.unkProps[:BeamEnergy])
    zs = zafcorrection(mc, fc, cc, kr.standard, lines, kr.stdProps[:BeamEnergy])
    return gZAFc(zu, zs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle])
end

"""
    k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)

The computed k-ratio for the unknown relative to standard.
"""
function k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    elm = element(unk)
    return (nonneg(material(unk), elm) / nonneg(material(std), elm)) * gZAFc(unk, std, θunk, θstd)
end

"""
    asa(::Type{DataFrame}, unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard in a DataFrame.
"""
function NeXLUncertainties.asa(
    ::Type{DataFrame},
    unk::MultiZAF,
    std::MultiZAF,
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame
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
    detail(::Type{DataFrame}, unk::MultiZAF, std::MultiZAF)::DataFrame

Tabulate each term in the MultiZAF matrix correction in a DataFrame.
"""
function detail(::Type{DataFrame}, unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame
    stds, stdE0, unks = Vector{String}(), Vector{Float64}(), Vector{String}()
    unkE0, xray, g = Vector{Float64}(), Vector{CharXRay}(), Vector{Float64}()
    z, a, f = Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    c, wgt, zaf, k = Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
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
        E0unk = 0.001 * unkE0,
        Standard = stds,
        E0std = 0.001 * stdE0,
        Xray = xray,
        Weight = wgt,
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
    detail(::Type{DataFrame}, mzs::AbstractArray{Tuple{MultiZAF, MultiZAF}})::DataFrame

Tabulate the details of a matrix correction relative to the specified unknown and standard in a DataFrame.
"""
detail(::Type{DataFrame}, mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat) =
    mapreduce(tmm -> detail(tmm[1], tmm[2], θunk, θstd), append!, mzs)

"""
    asa(::Type{DataFrame}, mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to a specified Dict of unknowns and standards in a DataFrame.
"""
NeXLUncertainties.asa(
    ::Type{DataFrame},
    mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}},
    θunk::AbstractFloat,
    θstd::AbstractFloat,
) = mapreduce(tmm -> asa(DataFrame, tmm[1], tmm[2], θunk, θstd), append!, mzs)


function NeXLUncertainties.asa(
    ::Type{DataFrame},
    unk::Material,
    std::Material,
    lines::Vector{CharXRay},
    e0::Float64,
    toa::Float64;
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
    coating::Type{<:CoatingCorrection} = Coating,
)
    flines = filter(cxr->energy(inner(cxr)) < e0, lines)
    zafs = zafcorrection(mc, fc, coating, unk, std, flines, e0)
    asa(DataFrame, zafs..., toa, toa)
end

function aspure(c::Material, cxr::CharXRay, e0::Float64, toa::Float64)
    zs = zafcorrection(XPP, ReedFluorescence, NullCoating, unk, pure(element(cxr)), [cxr], e0)
    return k(zs..., toa, toa)
end
