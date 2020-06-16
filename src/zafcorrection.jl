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
    transmission(std.coating, cxr, θunk) / transmission(unk.coating, cxr, θstd)

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
    NeXLUncertainties.asa(::Type{DataFrame}, unk::ZAFCorrection, std::ZAFCorrection, trans::AbstractVector{Transition},
    θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard for the iterable of Transition, trans.
"""
function NeXLUncertainties.asa(#
    ::Type{DataFrame},
    unk::ZAFCorrection,
    std::ZAFCorrection,
    trans::AbstractVector{Transition},
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame
    @assert isequal(atomicsubshell(unk.za), atomicsubshell(std.za))
    "The atomic sub-shell for the standard and unknown don't match."
    cxrs = characteristic(
        element(atomicsubshell(unk.za)),
        trans,
        1.0e-9,
        0.999 * min(beamEnergy(unk.za), beamEnergy(std.za)),
    )
    stds, stdE0, unks, unkE0, xray =
        Vector{String}(), Vector{Float64}(), Vector{String}(), Vector{Float64}(), Vector{CharXRay}()
    z, a, f, c, zaf, k, unkToa, stdToa = Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Float64}()
    for cxr in cxrs
        if isequal(inner(cxr), atomicsubshell(std.za))
            elm = element(cxr)
            push!(unks, name(material(unk.za)))
            push!(unkE0, beamEnergy(unk.za))
            push!(unkToa, rad2deg(θunk))
            push!(stds, name(material(std.za)))
            push!(stdE0, beamEnergy(std.za))
            push!(stdToa, rad2deg(θstd))
            push!(xray, cxr)
            push!(z, Z(unk, std))
            push!(a, A(unk, std, cxr, θunk, θstd))
            push!(f, F(unk, std, cxr, θunk, θstd))
            push!(c, coating(unk, std, cxr, θunk, θstd))
            tot = ZAFc(unk, std, cxr, θunk, θstd)
            push!(zaf, tot)
            push!(k, tot * material(unk.za)[elm] / material(std.za)[elm])
        end
    end
    return DataFrame(
        Unknown = unks,
        E0unk = unkE0,
        TOAunk = unkToa,
        Standard = stds,
        E0std = stdE0,
        TOAstd = stdToa,
        Xray = xray,
        Z = z,
        A = a,
        F = f,
        c = c,
        ZAF = zaf,
        k = k,
    )
end

NeXLUncertainties.asa( #
    ::Type{DataFrame},
    unk::ZAFCorrection,
    std::ZAFCorrection,
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame = asa(DataFrame, unk, std, alltransitions, θunk, θstd)

function NeXLUncertainties.asa( #
    ::Type{DataFrame},
    zafs::Dict{ZAFCorrection,ZAFCorrection},
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame
    df = DataFrame()
    for (unk, std) in zafs
        append!(df, asa(DataFrame, unk, std, θunk, θstd))
    end
    return df
end

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
    norm = asnormalized(mat)  # Ensures convergence of the interation algorithms...
    return ZAFCorrection(
        matrixcorrection(mctype, norm, ashell, e0),
        fluorescencecorrection(fctype, norm, ashell, e0),
        coatingcorrection(cctype, coating),
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

function NeXLUncertainties.asa( #
    ::Type{DataFrame},
    unk::Material,
    stds::Dict{Element,Material},
    e0::Real,
    θunk::AbstractFloat,
    θstd::AbstractFloat;
    mctype::Type{<:MatrixCorrection} = XPP,
    fctype::Type{<:FluorescenceCorrection} = ReedFluorescence,
    cctype::Type{<:CoatingCorrection} = Coating,
)
    df = DataFrame()
    for (elm, std) in stds
        for ashell in atomicsubshells(elm, e0)
            append!(
                df,
                asa( #
                    DataFrame,
                    zafcorrection(mctype, fctype, cctype, unk, std, ashell, e0)...,
                    alltransitions,
                    θunk,
                    θstd,
                ),
            )
        end
    end
    return df
end
