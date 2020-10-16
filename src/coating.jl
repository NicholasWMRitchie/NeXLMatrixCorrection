"""
An abstract type for implementing coating correction algorithms.

    NeXLCore.transmission(zaf::CoatingCorrection, xray::CharXRay)
"""
abstract type CoatingCorrection end

NeXLCore.minproperties(::Type{<:CoatingCorrection}) = ( )  # None are necessary but :Coating is optional

"""
    Coating
Implements a simple multi-layer coating correction.
"""
struct Coating <: CoatingCorrection
    layers::Vector{Film}
    Coating(coatings::AbstractVector{Film}) = new(coatings)
end

coatingcorrection(::Type{Coating}, film::Film) = Coating([film])
coatingcorrection(::Type{Coating}, films::AbstractVector{Film}) = Coating(films)

NeXLCore.minproperties(::Type{Coating}) = ( :Coating )  # :Coating = Film() | Film[]

"""
    carboncoating(nm)

Constructs a carbon coating of the specified thickness (in nanometers).
"""
carboncoating(nm) = Film(pure(n"C"), nm * 1.0e-7)

"""
    NeXLCore.transmission(zaf::Coating, xray::CharXRay, toa)

Calculate the transmission fraction for the specified X-ray through the coating
in the direction of the detector.
"""
NeXLCore.transmission(cc::Coating, xray::CharXRay, θtoa::AbstractFloat) =
    mapreduce(lyr -> transmission(lyr, xray, θtoa), *, cc.layers, init = 1.0)

Base.show(io::IO, coating::Coating) = length(coating.layers) == 1 ? #
    (thickness(coating.layers[1]) <= 0.0 ? Base.show(io, "No coating") : Base.show(io, coating.layers[1])) : #
    Base.show(io, coating.layers)

"""
    NullCoating

No coating (100%) transmission.
"""
struct NullCoating <: CoatingCorrection
    NullCoating(films::AbstractVector{Film}) = new()
    NullCoating(film::Film) = new()
    NullCoating() = new()
end

coatingcorrection(::Type{<:CoatingCorrection}, mss::Missing) = NullCoating()
coatingcorrection(::Type{NullCoating}, film::Film) = NullCoating()
coatingcorrection(::Type{NullCoating}, film::AbstractVector{Film}) = NullCoating()

"""
    NeXLCore.transmission(zaf::NullCoating, xray::CharXRay, toa)

Calculate the transmission fraction for the specified X-ray through no coating (1.0).
"""
NeXLCore.transmission(nc::NullCoating, xray::CharXRay, θtoa::AbstractFloat) = 1.0

Base.show(io::IO, nc::NullCoating) = print(io, "No coating.")
