# Generic methods for matrix correctio

using DataFrames
using PeriodicTable

"""
MatrixCorrection structs should implement

    F(mc::MatrixCorrection)
    Fχ(mc::MatrixCorrection, xray::CharXRay, θtoa::AbstractFloat)
    atomicsubshell(mc::MatrixCorrection)
    material(mc::MatrixCorrection)
    beamEnergy(mc::MatrixCorrection)
    ϕ(ρz)
    ϕabs(ρz, θtoa)
"""
abstract type MatrixCorrection end

"""
    NullCorrection
Implements Castaing's First Approximation
"""
struct NullCorrection <: MatrixCorrection
    material::Material
    subshell::AtomicSubShell
    E0::AbstractFloat

    NullCorrection(mat::Material, ashell::AtomicSubShell, e0) = new(mat, ashell, e0)
end

Base.show(io::IO, nc::NullCorrection) =
    print(io, "Unity[" + nc.material, ", ", subshell, ", ", 0.001 * e0, " keV]")

F(mc::NullCorrection) = 1.0
Fχ(mc::NullCorrection, xray::CharXRay, θtoa::AbstractFloat) = 1.0
NeXLCore.atomicsubshell(mc::NullCorrection) = mc.subshell
NeXLCore.material(mc::NullCorrection) = mc.material
beamEnergy(mc::NullCorrection) = mc.E0 # in eV

"""
    χ(mat::Material, xray::CharXRay, θtoa)
Angle adjusted mass absorption coefficient.
"""
χ(mat::Material, xray::CharXRay, θtoa::AbstractFloat) =
    mac(mat, xray) * csc(θtoa)

"""
    ZA(
      unk::MatrixCorrection,
      std::MatrixCorrection,
      xray::CharXRay,
      θunk::AbstractFloat,
      θstd::AbstractFloat
    )

The atomic number and absorption correction factors.
"""
function ZA(
    unk::MatrixCorrection,
    std::MatrixCorrection,
    xray::CharXRay,
    θunk::AbstractFloat,
    θstd::AbstractFloat
)
    @assert isequal(unk.subshell, inner(xray)) "Unknown and X-ray don't match in XPP"
    @assert isequal(std.subshell, inner(xray)) "Standard and X-ray don't match in XPP"
    return Fχ(unk, xray, θunk) / Fχ(std, xray, θstd)
end

"""
    Z(unk::MatrixCorrection, std::MatrixCorrection)
The atomic number correction factor.
"""
function Z(unk::MatrixCorrection, std::MatrixCorrection)
    @assert isequal(unk.subshell, std.subshell)
        "Unknown and standard matrix corrections don't apply to the same sub-shell."
    return F(unk) / F(std)
end

"""
    A(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay, χcunk=0.0, tcunk=0.0, χcstd=0.0, tcstd=0.0)
The absorption correction factors.
"""
function A(
    unk::MatrixCorrection,
    std::MatrixCorrection,
    xray::CharXRay,
    θunk::AbstractFloat,
    θstd::AbstractFloat
)
    @assert isequal(unk.subshell, inner(xray)) "Unknown and X-ray don't match in XPP"
    @assert isequal(std.subshell, inner(xray)) "Standard and X-ray don't match in XPP"
    return ZA(unk, std, xray, θunk, θstd) / Z(unk, std)
end

"""
Implements

    F(unk::FluorescenceCorrection, xray::CharXRay, θtoa::AbstractFloat)
"""
abstract type FluorescenceCorrection end

"""
    NullFluorescence

Implements Castaing's First Approximation
"""
struct NullFluorescence <: FluorescenceCorrection

    NullFluorescence(mat::Material, ashell::AtomicSubShell, e0::AbstractFloat) =
        new()
end

Base.show(io::IO, nc::NullFluorescence) =
    print(io, "Null[Fluor]")

F(nc::NullFluorescence, cxr::CharXRay, θtoa::AbstractFloat) = 1.0

"""
    CoatingCorrection
Implements
    transmission(zaf::CoatingCorrection, xray::CharXRay)
"""
abstract type CoatingCorrection end

"""
    Coating
Implements a simple single layer coating correction.
"""
struct Coating <: CoatingCorrection
    layer::Film
    Coating(mat::Material, thickness::AbstractFloat) =
        new(Film(mat, thickness))

end

"""
    carboncoating(nm)

Constructs a carbon coating of the specified thickness (in nanometers).
"""
carboncoating(nm) = Coating(pure(n"C"), nm * 1.0e-7)

"""
    transmission(zaf::Coating, xray::CharXRay, toa)
Calculate the transmission fraction for the specified X-ray through the coating
in the direction of the detector.
"""
transmission(cc::Coating, xray::CharXRay, θtoa::AbstractFloat) =
    transmission(layer, xray, θtoa)

Base.show(io::IO, coating::Coating) =
    Base.show(io, coating.layer)

"""
    NullCoating
No coating (100%) transmission.
"""
struct NullCoating <: CoatingCorrection end

"""
    transmission(zaf::NullCoating, xray::CharXRay, toa)
Calculate the transmission fraction for the specified X-ray through no coating (1.0).
"""
transmission(nc::NullCoating, xray::CharXRay, θtoa::AbstractFloat) = 1.0

Base.show(io::IO, nc::NullCoating) = print(io, "no coating")

"""
    ZAFCorrection
Pulls together the ZA, F and coating corrections into a single structure.
"""
struct ZAFCorrection
    za::MatrixCorrection
    f::FluorescenceCorrection
    coating::CoatingCorrection

    ZAFCorrection(
        za::MatrixCorrection,
        f::FluorescenceCorrection,
        coating::CoatingCorrection = NullCoating(),
    ) = new(za, f, coating)
end

NeXLCore.material(zaf::ZAFCorrection) = material(zaf.za)

Z(unk::ZAFCorrection, std::ZAFCorrection) = Z(unk.za, std.za)

A(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    A(unk.za, std.za, cxr, θunk, θstd)

ZA(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    ZA(unk.za, std.za, cxr, θunk, θstd)

coating(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    cxr::CharXRay,
    θunk::AbstractFloat,
    θstd::AbstractFloat
) = transmission(unk.coating, cxr, θunk) / transmission(std.coating, cxr, θstd)

generation(unk::ZAFCorrection, std::ZAFCorrection, ass::AtomicSubShell) =
    ionizationcrosssection(ass, beamEnergy(unk))/ionizationcrosssection(ass, beamEnergy(std))

F(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θunk::AbstractFloat, θstd::AbstractFloat) =
    F(unk.f, cxr, θunk) / F(std.f, cxr, θstd)

ZAFc(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    cxr::CharXRay,
    θunk::AbstractFloat,
    θstd::AbstractFloat
) =
    Z(unk, std) * A(unk, std, cxr, θunk, θstd) * F(unk, std, cxr, θunk, θstd) *
    coating(unk, std, cxr, θunk, θstd)

beamEnergy(zaf::ZAFCorrection) = beamEnergy(zaf.za)

Base.show(io::IO, cc::ZAFCorrection) =
    print(io, "ZAF[", cc.za, ", ", cc.f, ", ", cc.coating, "]")

"""
    zaf(za::MatrixCorrection, f::FluorescenceCorrection, coating::CoatingCorrection = NullCoating())
Construct a ZAFCorrection object
"""
zaf(
    za::MatrixCorrection,
    f::FluorescenceCorrection,
    coating::CoatingCorrection = NullCoating(),
) = ZAFCorrection(za, f, coating)

"""
    NeXLUncertainties.asa(::Type{DataFrame}, unk::ZAFCorrection, std::ZAFCorrection, trans)::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard
for the iterable of Transition, trans.
"""
function NeXLUncertainties.asa(#
    ::Type{DataFrame},
    unk::ZAFCorrection,
    std::ZAFCorrection,
    trans,
    θunk::AbstractFloat,
    θstd::AbstractFloat
)::DataFrame
    @assert isequal(atomicsubshell(unk.za), atomicsubshell(std.za))
        "The atomic sub-shell for the standard and unknown don't match."
    cxrs = characteristic(
        element(atomicsubshell(unk.za)),
        trans,
        1.0e-9,
        0.999 * min(beamEnergy(unk.za), beamEnergy(std.za)),
    )
    stds, stdE0, unks, unkE0, xray = Vector{String}(),
        Vector{Float64}(),
        Vector{String}(),
        Vector{Float64}(),
        Vector{CharXRay}()
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

tabulate( #
    unk::ZAFCorrection,
    std::ZAFCorrection,
    θunk::AbstractFloat,
    θstd::AbstractFloat
)::DataFrame = tabulate(unk, std, alltransitions, θunk, θstd)

function tabulate( #
    zafs::Dict{ZAFCorrection,ZAFCorrection},
    θunk::AbstractFloat,
    θstd::AbstractFloat
)::DataFrame
    df = DataFrame()
    for (unk, std) in zafs
        append!(df, tabulate(unk, std, θunk, θstd))
    end
    return df
end

"""
    ZAF(
      mctype::Type{<:MatrixCorrection},
      fctype::Type{<:FluorescenceCorrection},
      mat::Material,
      ashell::AtomicSubShell,
      e0,
      coating=NullCoating()
    )

Constructs an ZAFCorrection object using the mctype correction model with
the fluorescence model for the specified parameters.
"""
function ZAF(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    mat::Material,
    ashell::AtomicSubShell,
    e0::AbstractFloat,
    coating = NullCoating(),
)
    norm = asnormalized(mat)  # Ensures convergence of the interation algorithms...
    return ZAFCorrection(
        matrixcorrection(mctype, norm, ashell, e0),
        fluorescencecorrection(fctype, norm, ashell, e0),
        coating,
    )
end

"""
    ZAF(
       mctype::Type{<:MatrixCorrection},
       fctype::Type{<:FluorescenceCorrection},
       unk::Material,
       std::Material,
       ashell::AtomicSubShell,
       e0::AbstractFloat;
       unkCoating = NullCoating(),
       stdCoating = NullCoating(),
    )

Creates a matched pair of ZAFCorrection objects using the matrix correction algorithm
for the specified unknown and standard.
"""
ZAF(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    unk::Material,
    std::Material,
    ashell::AtomicSubShell,
    e0::AbstractFloat;
    unkCoating = NullCoating(),
    stdCoating = NullCoating()
) = (
    ZAF(mctype, fctype, unk, ashell, e0, unkCoating),
    ZAF(mctype, fctype, std, ashell, e0, stdCoating),
)

"""
    ZAF(
      mctype::Type{<:MatrixCorrection},
      fctype::Type{<:FluorescenceCorrection},
      mat::Material,
      cxrs,
      e0,
      coating=NeXLCore.NullCoating()
    )

Constructs a MultiZAF around the mctype and fctype algorithms for a collection of CharXRay <code>cxrs</code>.
"""
function ZAF(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    mat::Material,
    cxrs,
    e0::AbstractFloat,
    coating = NullCoating()
)
    zafs = Dict((sh, ZAF(mctype, fctype, mat, sh, e0, coating)) for sh in union(inner.(cxrs)))
    return MultiZAF(cxrs, zafs)
end

"""
    ZAF(
      mctype::Type{<:MatrixCorrection},
      fctype::Type{<:FluorescenceCorrection},
      unk::Material,
      std::Material,
      cxrs,
      e0;
      unkCoating = NullCoating(),
      stdCoating = NullCoating(),
    )

Constructs a tuple of MultiZAF around the mctype and fctype correction algorithms for the unknown and standard for a
collection of CharXRay <code>cxrs</code>.
"""
ZAF(
    mctype::Type{<:MatrixCorrection},
    fctype::Type{<:FluorescenceCorrection},
    unk::Material,
    std::Material,
    cxrs,
    e0::AbstractFloat;
    unkCoating = NullCoating(),
    stdCoating = NullCoating(),
) = ( ZAF(mctype, fctype, unk, cxrs, e0, unkCoating),  ZAF(mctype, fctype, std, cxrs, e0, stdCoating))
