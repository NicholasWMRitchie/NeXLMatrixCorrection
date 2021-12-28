"""
An abstract type for implementing secondary-fluorescence corrections.

    F(unk::FluorescenceCorrection, xray::CharXRay, θtoa::AbstractFloat)
"""
abstract type FluorescenceCorrection end

"""
    NullFluorescence

Implements Castaing's First Approximation (i.e. No correction F=1)
"""
struct NullFluorescence <: FluorescenceCorrection

    NullFluorescence(mat::Material, ashell::AtomicSubShell, e0::AbstractFloat) = new()
end

Base.show(io::IO, nc::NullFluorescence) = print(io, "Null[Fluor]")

F(nc::NullFluorescence, cxr::CharXRay, θtoa::AbstractFloat) = 1.0

"""
    fluorescence(fltype::Type{<:FluorescenceCorrection}, comp::Material, secondary::AtomicSubShell, e0::Float64)

Construct an instance of a fltype correction structure to compute the
secondary fluorescence in the specified material and beam energy.
"""
function fluorescencecorrection(
    fltype::Type{<:FluorescenceCorrection},
    comp::Material,
    secondary::AtomicSubShell,
    e0::Float64;
    eThresh = 2.5e3,
    wThresh = 0.01,
)
    test(cxr, ee, wt) =
        (energy(cxr) > ee) && (energy(cxr) < ee + eThresh) && (NeXLCore.edgeenergy(cxr) < e0) && (weight(NormalizeToUnity, cxr) > wt)
    char4elm(elm, wt) = characteristic(elm, alltransitions, cxr -> test(cxr, energy(secondary), wt / comp[elm]))
    primaries = mapreduce(elm -> char4elm(elm, wThresh), append!, keys(comp))
    return fluorescencecorrection(fltype, comp, primaries, secondary, e0)
end

fluorescencecorrection(
    fltype::Type{NullFluorescence},
    comp::Material,
    secondary::AtomicSubShell,
    e0::Float64;
    eThresh = 2.5e3,
    wThresh = 0.01, 
) = NullFluorescence(comp, secondary, e0)


fluorescencecorrection(
    ::Type{NullFluorescence},
    comp::Material,
    primarys::Vector{CharXRay},
    secondary::AtomicSubShell,
    e0::Float64,
) = NullFluorescence(comp, secondary, e0)

