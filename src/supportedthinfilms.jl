"""
    SupportedThinFilms

A struct representing a supported thin film sample - a series of thin Film(s) over a substrate Material.
"""
struct SupportedThinFilms
    layer::Vector{Film} # Surface to interior
    substrate::Material # "Infinite" substrate

    function SupportedThinFilms(films::Vector{Film}, substrate::Material)
        @assert all(f->haskey(material(f), :Density), films) "Each layer must define the material :Density property."
        @assert haskey(substrate, :Density) "The substrate must define the material :Density property."
        return new(films, substrate)
    end
end

NeXLCore.massthickness(stf::SupportedThinFilms, i::Int) = NeXLCore.massthickness(stf.layer[i])

"""
    nlayers(stf::SupportedThinFilms)

Number of layers including the substrate.
"""
nlayers(stf::SupportedThinFilms) = length(stf.layer)+1

"""
    eachindex(stf::SupportedThinFilms)

Helpful to Iterate through each layer from surface to substrate.
"""
Base.eachindex(stf::SupportedThinFilms) = Base.OneTo(nlayers(stf))

"""
   outer(stf::SupportedThinFilms, i::Int)

Depth in cm of the side of the layer closer to the surface.
"""
function NeXLCore.outer(stf::SupportedThinFilms, i::Int)
    return i==1 ? 0 : sum(massthickness(l) for l in stf.layer[1:i-1])
end
"""
   inner(stf::SupportedThinFilms, i::Int)

Depth in cm of the side of the layer further from the surface.
"""
NeXLCore.inner(stf::SupportedThinFilms, i::Int) = i==nlayers(stf) ? Inf : outer(stf,i+1)

"""
    NeXLCore.material(stf::SupportedThinFilms, i::Int)

The Material from which the i-th layer (or substrate) is constructed.
"""
NeXLCore.material(stf::SupportedThinFilms, i::Int) = i <= length(stf.layer) ? material(stf.layer[i]) : stf.substrate

"""
    χs(stf::SupportedThinFilms, cxr::CharXRay, toa::Float64)::Vector{Float64}

Compute the reduced mass-absorption coefficient for each layer in the SupportedThinFilm.
"""
χs(stf::SupportedThinFilms, cxr::CharXRay, toa::Float64)::Vector{Float64} = #
    [ χ(material(stf,i), cxr, toa) for i in eachindex(stf) ]

"""
    Base.indexin(stf::SupportedThinFilms, ρz::Float64)

Returns the index of the layer in which the depth ρz is located.
"""
Base.indexin(stf::SupportedThinFilms, ρz::Float64) = findfirst(i->ρz < inner(stf,i), eachindex(stf))

"""
    transmission(stf::SupportedThinFilms, χs::Vector{Float64}, ρz::Float64)

The fraction of the intensity generated at ρz that will be emitted. `χs` are the reduced mass-absorption coefficients
which define the angle of exit and the mass-absorption coefficients.
"""
function NeXLCore.transmission(stf::SupportedThinFilms, χs::Vector{Float64}, ρz::Float64)
    @assert ρz >= 0.0
    k = indexin(stf, ρz)
    return (k <= 1 ? 1.0 : exp(-sum(χs[i]*massthickness(stf,i) for i in 1:k-1)))*exp(-(ρz-outer(stf, k))*χs[k])
end
