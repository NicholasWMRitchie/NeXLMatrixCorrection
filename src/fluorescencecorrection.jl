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

fluorescencecorrection(
    ::Type{NullFluorescence},
    comp::Material,
    primarys::Vector{CharXRay},
    secondary::AtomicSubShell,
    e0::Float64,
) = NullFluorescence(comp, secondary, e0)
