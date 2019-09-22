# Generic methods for matrix correctio

using DataFrames
using PeriodicTable

"""
MatrixCorrection structs should implement

    F(mc::MatrixCorrection)
    Fχ(mc::MatrixCorrection, xray::CharXRay)
    atomicshell(mc::MatrixCorrection)
    takeOffAngle(mc::MatrixCorrection)
    material(mc::MatrixCorrection)
    beamEnergy(mc::MatrixCorrection)
    ϕ(ρz)
    ϕabs(ρz)
"""
abstract type MatrixCorrection end

"""
    NullCorrection
Implements Castaing's First Approximation
"""
struct NullCorrection <: MatrixCorrection
    material::Material
    shell::AtomicShell
    θtoa
    E0
    NullCorrection(mat::Material, ashell::AtomicShell, toa, e0) =
        new(mat, shell, toa, e0)
end

Base.show(io::IO, nc::NullCorrection) =
    print(io,"Unity["+nc.material,", ",shell,", ",e0," keV, ",180.0*θtoa/π,"°]")

F(mc::NullCorrection) = 1.0
Fχ(mc::NullCorrection, xray::CharXRay) = 1.0
NeXLCore.atomicshell(mc::NullCorrection) = mc.shell
takeOffAngle(mc::NullCorrection) = mc.θtoa
NeXLCore.material(mc::NullCorrection) = mc.material
beamEnergy(mc::NullCorrection) = mc.E0

"""
    χ(mat::Material, xray::CharXRay, θtoa)
Angle adjusted mass absorption coefficient.
"""
χ(mat::Material, xray::CharXRay, θtoa) =
    mac(mat, xray) * csc(θtoa)

"""
    ZA(unk::XPP, std::XPP, xray::CharXRay, χcunk=0.0, tcunk=0.0, χcstd=0.0, tcstd=0.0)
The atomic number and absorption correction factors.
"""
function ZA(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay)
    @assert(isequal(unk.shell,inner(xray)),"Unknown and X-ray don't match in XPP")
    @assert(isequal(std.shell,inner(xray)),"Standard and X-ray don't match in XPP")
    return Fχ(unk, xray) / Fχ(std, xray)
end

"""
    Z(unk::XPP, std::XPP)
The atomic number correction factor.
"""
function Z(unk::MatrixCorrection, std::MatrixCorrection)
    @assert(isequal(unk.shell,std.shell),"Unknown and standard matrix corrections don't apply to the same shell.")
    return F(unk) / F(std)
end

"""
    A(unk::XPP, std::XPP, xray::CharXRay, χcunk=0.0, tcunk=0.0, χcstd=0.0, tcstd=0.0)
The absorption correction factors.
"""
function A(unk::MatrixCorrection, std::MatrixCorrection, xray::CharXRay)
    @assert(isequal(unk.shell,inner(xray)),"Unknown and X-ray don't match in XPP")
    @assert(isequal(std.shell,inner(xray)),"Standard and X-ray don't match in XPP")
    return ZA(unk, std, xray) / Z(unk, std)
end

"""
Implements
    F(unk::FluorescenceCorrection, xray::CharXRay)
"""
abstract type FluorescenceCorrection end

"""
    NullCorrection
Implements Castaing's First Approximation
"""
struct NullFluorescence <: FluorescenceCorrection
    material::Material
    shell::AtomicShell
    E0
    θtoa
    NullFluorescence(mat::Material, ashell::AtomicShell, e0, toa) =
        new(mat, ashell, e0, deg2rad(toa))
end

Base.show(io::IO, nc::NullFluorescence) =
    print(io,"Null[",name(nc.material),", ",nc.shell,", ",nc.E0," keV, ",rad2deg(nc.θtoa),"°]")

F(nc::NullFluorescence, cxr::CharXRay) = 1.0

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

carbonCoating(nm) = Coating(pure(n"C"),nm*1.0e-7)

"""
    transmission(zaf::Coating, xray::CharXRay, toa)
Calculate the transmission fraction for the specified X-ray through the coating
in the direction of the detector.
"""
transmission(cc::Coating, xray::CharXRay, toa) =
    exp(-χ(cc.coating, xray, toa)*cc.thickness)

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
transmission(nc::NullCoating, xray::CharXRay, toa) = 1.0

Base.show(io::IO, nc::NullCoating) =
    print(io, "no coating")

"""
    ZAFCorrection
Pulls together the ZA, F and coating corrections into a single structure.
"""
struct ZAFCorrection
    za::MatrixCorrection
    f::FluorescenceCorrection
    coating::CoatingCorrection
    
    ZAFCorrection(za::MatrixCorrection, f::FluorescenceCorrection, coating::CoatingCorrection = NullCoating()) =
        new(za, f, coating)
end

Z(unk::ZAFCorrection, std::ZAFCorrection) = Z(unk.za, std.za)
A(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay) =
    A(unk.za, std.za, cxr)
coating(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay) =
    transmission(unk.coating,cxr,takeOffAngle(unk.za))/
       transmission(std.coating,cxr,takeOffAngle(std.za))
F(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay) =
    F(unk.f, cxr) / F(std.f, cxr)
ZAFc(unk::ZAFCorrection, std::ZAFCorrection, cxr::CharXRay) =
    Z(unk, std)*A(unk, std, cxr)*F(unk, std, cxr)*coating(unk, std, cxr)
NeXLCore.material(zaf::ZAFCorrection) = material(zaf.za)
beamEnergy(zaf::ZAFCorrection) = beamEnergy(zaf.za)
takeOffAngle(zaf::ZAFCorrection) = takeOffAngle(zaf.za)

Base.show(io::IO, cc::ZAFCorrection) =
    print(io, "ZAF[",cc.za,", ",cc.f,", ",cc.coating, "]")

"""
    zaf(za::MatrixCorrection, f::FluorescenceCorrection, coating::CoatingCorrection = NullCoating())
Construct a ZAFCorrection object
"""
zaf(za::MatrixCorrection, f::FluorescenceCorrection, coating::CoatingCorrection = NullCoating()) =
    ZAFCorrection(za, f, coating)

"""
    NeXLCore.summarize(unk::ZAFCorrection, std::ZAFCorrection, trans)::DataFrame

Summarize a matrix correction relative to the specified unknown and standard
for the iterable of Transition, trans.
"""
function NeXLCore.summarize(unk::ZAFCorrection, std::ZAFCorrection, trans)::DataFrame
    @assert(isequal(atomicshell(unk.za), atomicshell(std.za)), "The atomic shell for the standard and unknown don't match.")
    cxrs = characteristic(element(atomicshell(unk.za)), trans, 1.0e-9, 0.999*1000.0*min(beamEnergy(unk.za),beamEnergy(std.za)))
    stds, stdE0, unks, unkE0, xray = Vector{String}(), Vector{Float64}(), Vector{String}(), Vector{Float64}(), Vector{CharXRay}()
    z, a, f, c, zaf, k = Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    for cxr in cxrs
        if isequal(inner(cxr), atomicshell(std.za))
            elm = element(cxr)
            push!(unks, name(material(unk.za)))
            push!(unkE0, beamEnergy(unk.za))
            push!(stds, name(material(std.za)))
            push!(stdE0, beamEnergy(std.za))
            push!(xray, cxr)
            push!(z, Z(unk, std))
            push!(a, A(unk, std, cxr))
            push!(f, F(unk, std, cxr))
            push!(c, coating(unk, std, cxr))
            tot = ZAFc(unk, std, cxr)
            push!(zaf, tot)
            push!(k, tot * material(unk.za)[elm] / material(std.za)[elm])
        end
    end
    return DataFrame(Unknown=unks, E0unk=unkE0, Standard=stds, E0std=stdE0, Xray=xray, Z=z, A=a, F=f, c=c, ZAF=zaf, k=k)
end

NeXLCore.summarize(unk::ZAFCorrection, std::ZAFCorrection)::DataFrame =
    summarize(unk,std,alltransitions)

function NeXLCore.summarize(zafs::Dict{ZAFCorrection, ZAFCorrection})::DataFrame
    df = DataFrame()
    for (unk, std) in zafs
        append!(df, summarize(unk, std))
    end
    return df
end

function NeXLCore.summarize(unk::Material, stds::Dict{Element, Material}, e0, θtoa)
    df = DataFrame()
    for (elm, std) in stds
        for ashell in atomicshells(elm, 1000.0*e0)
            unkXPP=xppZAF(unk, ashell, e0, θtoa)
            stdXPP=xppZAF(std, ashell, e0, θtoa)
            append!(df,summarize(unkXPP,stdXPP,alltransitions))
        end
    end
    return df
end
