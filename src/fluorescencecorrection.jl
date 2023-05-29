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

    NullFluorescence(_::Material, _::AtomicSubShell, _::AbstractFloat) = new()
end

Base.show(io::IO, _::NullFluorescence) = print(io, "NullFluorescence[]")

F(nc::NullFluorescence, cxr::CharXRay, θtoa::AbstractFloat) = 1.0

"""
    fluorescence(fltype::Type{<:FluorescenceCorrection}, comp::Material, secondary::AtomicSubShell, e0::AbstractFloat)

Construct an instance of a fltype correction structure to compute the
secondary fluorescence in the specified material and beam energy.
"""
function fluorescencecorrection(
    fltype::Type{<:FluorescenceCorrection},
    comp::Material,
    secondary::AtomicSubShell,
    e0::AbstractFloat;
    eThresh = 2.5e3,
    wThresh = 0.01,
)
    es = energy(secondary)
    primaries = mapreduce(append!, keys(comp)) do elm 
        characteristic(elm, alltransitions) do cxr
            (energy(cxr) > es) && (energy(cxr) < es + eThresh) && #
            (NeXLCore.edgeenergy(cxr) < 0.8*e0) && (weight(NormalizeToUnity, cxr) * value(comp[elm]) > wThresh)
        end
    end
    return fluorescencecorrection(fltype, comp, primaries, secondary, e0)
end

fluorescencecorrection(
    ::Type{NullFluorescence},
    comp::Material,
    secondary::AtomicSubShell,
    e0::AbstractFloat;
    vargs...
) = NullFluorescence(comp, secondary, e0)


fluorescencecorrection(
    ::Type{NullFluorescence},
    comp::Material,
    primarys::Vector{CharXRay},
    secondary::AtomicSubShell,
    e0::AbstractFloat,
) = NullFluorescence(comp, secondary, e0)

