# Generic methods for matrix correctio

using DataFrames
using PeriodicTable

"""
MatrixCorrection structs should implement

    F(mc::MatrixCorrection)
    Fχ(mc::MatrixCorrection, xray::CharXRay, θtoa::AbstractFloat)
    atomicshell(mc::MatrixCorrection)
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
    shell::AtomicShell
    E0::AbstractFloat

    NullCorrection(mat::Material, ashell::AtomicShell, e0) = new(mat, shell, e0)
end

Base.show(io::IO, nc::NullCorrection) =
    print(io, "Unity[" + nc.material, ", ", shell, ", ", 0.001 * e0, " keV]")

F(mc::NullCorrection) = 1.0
Fχ(mc::NullCorrection, xray::CharXRay, θtoa::AbstractFloat) = 1.0
NeXLCore.atomicshell(mc::NullCorrection) = mc.shell
NeXLCore.material(mc::NullCorrection) = mc.material
beamEnergy(mc::NullCorrection) = mc.E0 # in eV

"""
    χ(mat::Material, xray::CharXRay, θtoa)
Angle adjusted mass absorption coefficient.
"""
χ(mat::Material, xray::CharXRay, θtoa::AbstractFloat) =
    mac(mat, xray) * csc(θtoa)

"""
    ZA(unk::XPP, std::XPP, xray::CharXRay, χcunk=0.0, tcunk=0.0, χcstd=0.0, tcstd=0.0)
The atomic number and absorption correction factors.
"""
function ZA(
    unk::MatrixCorrection,
    std::MatrixCorrection,
    xray::CharXRay,
    θtoa::AbstractFloat,
)
    @assert(
        isequal(unk.shell, inner(xray)),
        "Unknown and X-ray don't match in XPP",
    )
    @assert(
        isequal(std.shell, inner(xray)),
        "Standard and X-ray don't match in XPP",
    )
    return Fχ(unk, xray, θtoa) / Fχ(std, xray, θtoa)
end

"""
    Z(unk::XPP, std::XPP)
The atomic number correction factor.
"""
function Z(unk::MatrixCorrection, std::MatrixCorrection)
    @assert(
        isequal(unk.shell, std.shell),
        "Unknown and standard matrix corrections don't apply to the same shell.",
    )
    return F(unk) / F(std)
end

"""
    A(unk::XPP, std::XPP, xray::CharXRay, χcunk=0.0, tcunk=0.0, χcstd=0.0, tcstd=0.0)
The absorption correction factors.
"""
function A(
    unk::MatrixCorrection,
    std::MatrixCorrection,
    xray::CharXRay,
    θtoa::AbstractFloat,
)
    @assert(
        isequal(unk.shell, inner(xray)),
        "Unknown and X-ray don't match in XPP",
    )
    @assert(
        isequal(std.shell, inner(xray)),
        "Standard and X-ray don't match in XPP",
    )
    return ZA(unk, std, xray, θtoa) / Z(unk, std)
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

    NullFluorescence(mat::Material, ashell::AtomicShell, e0::AbstractFloat) =
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
    coating::Material
    thickness
end

"""
    carbonCoating(nm)

Constructs a carbon coating of the specified thickness (in nanometers).
"""
carbonCoating(nm) = Coating(pure(n"C"), nm * 1.0e-7)

"""
    transmission(zaf::Coating, xray::CharXRay, toa)
Calculate the transmission fraction for the specified X-ray through the coating
in the direction of the detector.
"""
transmission(cc::Coating, xray::CharXRay, θtoa::AbstractFloat) =
    exp(-χ(cc.coating, xray, θtoa) * cc.thickness)

Base.show(io::IO, coating::Coating) =
    print(io, 1.0e7 * coating.thickness, " nm of ", name(coating.coating))

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

Z(unk::ZAFCorrection, std::ZAFCorrection) = Z(unk.za, std.za)

A(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θtoa::AbstractFloat) =
    A(unk.za, std.za, cxr, θtoa)

coating(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    cxr::CharXRay,
    θtoa::AbstractFloat,
) = transmission(unk.coating, cxr, θtoa) / transmission(std.coating, cxr, θtoa)

F(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay, θtoa::AbstractFloat) =
    F(unk.f, cxr, θtoa) / F(std.f, cxr, θtoa)

ZAFc(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    cxr::CharXRay,
    θtoa::AbstractFloat,
) =
    Z(unk, std) * A(unk, std, cxr, θtoa) * F(unk, std, cxr, θtoa) *
    coating(unk, std, cxr, θtoa)

NeXLCore.material(zaf::ZAFCorrection) = material(zaf.za)
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
    NeXLCore.summarize(unk::ZAFCorrection, std::ZAFCorrection, trans)::DataFrame

Summarize a matrix correction relative to the specified unknown and standard
for the iterable of Transition, trans.
"""
function NeXLCore.summarize(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    trans,
    θtoa::AbstractFloat,
)::DataFrame
    @assert(
        isequal(atomicshell(unk.za), atomicshell(std.za)),
        "The atomic shell for the standard and unknown don't match.",
    )
    cxrs = characteristic(
        element(atomicshell(unk.za)),
        trans,
        1.0e-9,
        0.999 * min(beamEnergy(unk.za), beamEnergy(std.za)),
    )
    stds, stdE0, unks, unkE0, xray = Vector{String}(),
        Vector{Float64}(),
        Vector{String}(),
        Vector{Float64}(),
        Vector{CharXRay}()
    z, a, f, c, zaf, k, toa = Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}(),
        Vector{Float64}()
    for cxr in cxrs
        if isequal(inner(cxr), atomicshell(std.za))
            elm = element(cxr)
            push!(unks, name(material(unk.za)))
            push!(unkE0, beamEnergy(unk.za))
            push!(stds, name(material(std.za)))
            push!(stdE0, beamEnergy(std.za))
            push!(toa, rad2deg(θtoa))
            push!(xray, cxr)
            push!(z, Z(unk, std))
            push!(a, A(unk, std, cxr, θtoa))
            push!(f, F(unk, std, cxr, θtoa))
            push!(c, coating(unk, std, cxr, θtoa))
            tot = ZAFc(unk, std, cxr, θtoa)
            push!(zaf, tot)
            push!(k, tot * material(unk.za)[elm] / material(std.za)[elm])
        end
    end
    return DataFrame(
        Unknown = unks,
        E0unk = unkE0,
        Standard = stds,
        E0std = stdE0,
        TOA = toa,
        Xray = xray,
        Z = z,
        A = a,
        F = f,
        c = c,
        ZAF = zaf,
        k = k,
    )
end

NeXLCore.summarize(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    θtoa::AbstractFloat,
)::DataFrame = summarize(unk, std, alltransitions, θtoa)

function NeXLCore.summarize(
    zafs::Dict{ZAFCorrection,ZAFCorrection},
    θtoa::AbstractFloat,
)::DataFrame
    df = DataFrame()
    for (unk, std) in zafs
        append!(df, summarize(unk, std, θtoa))
    end
    return df
end
