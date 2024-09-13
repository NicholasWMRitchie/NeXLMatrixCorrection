"""
    ZAFCorrection

Pulls together the ZA, F and coating corrections into a single structure.
"""
struct ZAFCorrection
    za::MatrixCorrection
    f::FluorescenceCorrection
    coating::CoatingCorrection
    """
        ZAFCorrection(za::MatrixCorrection, f::FluorescenceCorrection, coating::CoatingCorrection)
    """
    ZAFCorrection(za::MatrixCorrection, f::FluorescenceCorrection, coating::CoatingCorrection) = new(za, f, coating)
end

NeXLCore.material(zaf::ZAFCorrection) = material(zaf.za)

"""
    Z(unk::ZAFCorrection, std::ZAFCorrection)

Computes the atomic number correction.
"""
Z(unk::ZAFCorrection, std::ZAFCorrection) = Z(unk.za, std.za)

"""
    A(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)

Computes the absorption correction.
"""
A(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    A(unk.za, std.za, cxr, θunk, θstd)

"""
    ZA(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)

Computes the combined atomic number and fluorescence correction.
"""
ZA(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    ZA(unk.za, std.za, cxr, θunk, θstd)

"""
    coating(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)

Computes the coating correction.
"""
coating(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    transmission(unk.coating, cxr, θunk) / transmission(std.coating, cxr, θstd)

"""
    generation(unk::ZAFCorrection, std::ZAFCorrection, ass::AtomicSubShell)

Computes a correction factor for differences X-ray generation due to differences in beam energy.
"""
generation(unk::ZAFCorrection, std::ZAFCorrection, ass::AtomicSubShell) =
    ionizationcrosssection(ass, beamEnergy(unk)) / ionizationcrosssection(ass, beamEnergy(std))

"""
    F(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)

Computes the secondary fluorescence correction.
"""
F(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    F(unk.f, cxr, θunk) / F(std.f, cxr, θstd)

"""
    ZAFc(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)

Computes the combined correction for atomic number, absorption, secondary fluorescence and generation.
"""
ZAFc(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    ZA(unk, std, cxr, θunk, θstd) * F(unk, std, cxr, θunk, θstd) * coating(unk, std, cxr, θunk, θstd)


"""
    gZAFc(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)

Computes the combined correction for generation, atomic number, absorption, and secondary fluorescence.
"""
gZAFc(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    generation(unk, std, inner(cxr)) * ZA(unk, std, cxr, θunk, θstd) * F(unk, std, cxr, θunk, θstd) * coating(unk, std, cxr, θunk, θstd)


"""
    k(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)

The computed k-ratio for the unknown relative to standard.
"""
function k(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat)
    elm = element(cxr)
    return (nonneg(material(unk), elm) / nonneg(material(std), elm)) * ZAFc(unk, std, cxr, θunk, θstd)
end

beamEnergy(zaf::ZAFCorrection) = beamEnergy(zaf.za)

Base.show(io::IO, cc::ZAFCorrection) = print(io, "ZAF[", cc.za, ", ", cc.f, ", ", cc.coating, "]")


"""
    zafcorrection(
      mctype::Type{<:MatrixCorrection},
      fctype::Type{<:FluorescenceCorrection},
      cctype::Type{<:CoatingCorrection},
      mat::Material,
      ashell::AtomicSubShell,
      e0::Real,
      coating::Union{Film,AbstractVector{Film},Missing}
    )

Constructs an ZAFCorrection object using the mctype correction model with
the fluorescence model for the specified parameters.
"""
function zafcorrection(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    cctype::Type{<:CoatingCorrection},
    mat::Material,
    ashell::AtomicSubShell,
    e0::Real,
    coating::Union{Film,AbstractVector{Film},Missing},
)
    deds(coating::Missing) = 0.0
    deds(coating::Film) = dEds(JoyLuo, e0, coating.material)*coating.thickness
    deds(coating::AbstractVector{Film}) = mapreduce(f->dEds(JoyLuo, e0, f.mat)*f.thicknes, +, coating)
    nn = analyticaltotal(Float64, mat) > 0
    norm = asnormalized(mat)  # Ensures convergence of the interation algorithms...
    e0c = e0 + deds(coating)
    return ZAFCorrection(
        matrixcorrection(nn ? mctype : NullCorrection, norm, ashell, e0c),
        fluorescencecorrection(nn ? fctype : NullFluorescence, norm, ashell, e0c),
        coatingcorrection(nn ? cctype : NullCoating, coating),
    )
end

"""
    zafcorrection(
       mctype::Type{<:MatrixCorrection},
       fctype::Type{<:FluorescenceCorrection},
       cctype::Type{<:CoatingCorrection},
       unk::Material,
       std::Material,
       ashell::AtomicSubShell,
       e0::Real;
       unkCoating::Union{Film,AbstractVector{Film},Missing} = missing,
       stdCoating::Union{Film,AbstractVector{Film},Missing} = missing,
    )

Creates a matched pair of ZAFCorrection objects using the matrix correction algorithm
for the specified unknown and standard.
"""
function zafcorrection(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    cctype::Type{<:CoatingCorrection},
    unk::Material,
    std::Material,
    ashell::AtomicSubShell,
    e0::Real;
    unkCoating::Union{Film,AbstractVector{Film},Missing} = missing,
    stdCoating::Union{Film,AbstractVector{Film},Missing} = missing,
)
    @assert energy(ashell) < e0 "The shell energy must be less than the beam energy."
    return (
        zafcorrection(mctype, fctype, cctype, unk, ashell, e0, unkCoating),
        zafcorrection(mctype, fctype, cctype, std, ashell, e0, stdCoating),
    )
end
