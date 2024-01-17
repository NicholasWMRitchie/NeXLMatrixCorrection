"""
The `MultiZAF` structure optimizes the calculation of the matrix correction 
for mutliple lines measured simultaneously from a single contiguous peak as
represented by a single k-ratio (as is often the case with energy dispersive
measurements.)  For example, when the potassium K line is measured by EDS, 
one k-ratio is measured that corresponds to the K K-L2, K K-L3, K K-M2 & 
K K-M3 transitions.  Most of the matrix correction calculation depends upon
the shell the "K K" shell and only a small part depends on the exact 
characteristic X-ray.  To optimize the calculation, the `ZAFCorrection` 
object is calculated once for each `AtomicSubShell` and reused to calculated
all the `CharXRay`s.

DOn't construct MultiZAF objects directly; instead use `zafcorrection(...)`.
"""
struct MultiZAF
    xrays::Vector{CharXRay}
    zafs::Dict{AtomicSubShell,ZAFCorrection}
    function MultiZAF(xrays, zafs)
        elm = element(first(xrays))
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
    coating::Union{Film,AbstractVector{Film},Missing}=missing,
)
    mat = asnormalized(convert(Material{Float64, Float64}, mat))
    shells = union(filter(sh -> energy(sh) < e0, inner.(cxrs)))
    zafs = Dict(sh => zafcorrection(mctype, fctype, cctype, mat, sh, e0, coating) for sh in shells)
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
    unkCoating::Union{Film,AbstractVector{Film},Missing}=missing,
    stdCoating::Union{Film,AbstractVector{Film},Missing}=missing
)
    xrays = filter(cxr -> energy(inner(cxr)) < e0, cxrs)
    return (
        zafcorrection(mctype, fctype, cctype, unk, xrays, e0, unkCoating),
        zafcorrection(mctype, fctype, cctype, std, xrays, e0, stdCoating),
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

function commonXrays(unk::MultiZAF, std::MultiZAF)
    cunk = characteristic(unk)
    filter(cxr -> cxr in cunk, characteristic(std))
end
"""
    Z(unk::MultiZAF, std::MultiZAF)

The Z (atomic number) correction for `unk` relative to `std`.
"""
function Z(unk::MultiZAF, std::MultiZAF)
    z, n = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        norm = sum(weight.(NormalizeByShell, cxrs2))
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
            w = weight(NormalizeByShell, cxr)
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
    n, f = 0.0, 0.0
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        for cxr in cxrs2
            w = weight(NormalizeByShell, cxr)
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
            w = weight(NormalizeByShell, cxr)
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
            c += weight(NormalizeByShell, cxr) * coating(zafU, zafS, cxr, θunk, θstd)
            n += weight(NormalizeByShell, cxr)
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
            w = weight(NormalizeByShell, cxr)
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
    mc::Type{<:MatrixCorrection}=XPP,
    fc::Type{<:FluorescenceCorrection}=ReedFluorescence,
    cc::Type{<:CoatingCorrection}=Coating
)
    xrays = filter(cxr -> energy(inner(cxr)) < min(kr.unkProps[:BeamEnergy], kr.stdProps[:BeamEnergy]), kr.xrays)
    zu = zafcorrection(mc, fc, cc, unkComp, xrays, kr.unkProps[:BeamEnergy])
    zs = zafcorrection(mc, fc, cc, kr.standard, xrays, kr.stdProps[:BeamEnergy])
    return gZAFc(zu, zs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle])
end

"""
    k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)

Compute the k-ratio for the unknown relative to standard.
"""
function k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)
    elm = element(unk)
    return (nonneg(material(unk), elm) / nonneg(material(std), elm)) * gZAFc(unk, std, θunk, θstd)
end

function aspure(unk::Material, cxr::CharXRay, e0::Float64, toa::Float64)
    zs = zafcorrection(XPP, ReedFluorescence, NullCoating, unk, pure(element(cxr)), [cxr], e0)
    return k(zs..., toa, toa)
end


"""
    NeXLCore.KRatio(
        cxrs::Vector{CharXRay}, 
        unk_mat::Material, 
        unk_props::Dict{Symbol,Any}, 
        std_mat::Material, 
        std_props::Dict{Symbol,Any};
        mc = XPP,
        fc = ReedFluorescence,
        cc = Coating
    )

Construct a KRatio object for the specified measurement.
The :BeamEnergy and :TakeOffAngle properties are required.
The :Coating property is optional.
"""
function NeXLCore.KRatio(
    cxrs::Vector{CharXRay},
    unk_mat::Material,
    unk_props::Dict{Symbol,Any},
    std_mat::Material,
    std_props::Dict{Symbol,Any};
    mc::Type{<:MatrixCorrection}=XPP,
    fc::Type{<:FluorescenceCorrection}=ReedFluorescence,
    cc::Type{<:CoatingCorrection}=Coating
)
    KRatio(
        cxrs,
        unk_props,
        std_props,
        std_mat,
        k(
            zafcorrection(mc, fc, cc, unk_mat, cxrs, unk_props[:BeamEnergy], get(unk_props, :Coating, missing)),
            zafcorrection(mc, fc, cc, std_mat, cxrs, std_props[:BeamEnergy], get(std_props, :Coating, missing)),
            unk_props[:TakeOffAngle], std_props[:TakeOffAngle]
        )
    )
end
