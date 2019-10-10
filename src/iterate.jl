

struct KRatio
    lines::Array{CharXRay}
    properties::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratio::AbstractFloat
end

abstract type UnmeasredElementRule end

struct Iteration
    mctype::Type{MatrixCorrection}
    fctype::Type{FluorescenceCorrection}
    kratios::Array{KRatio}
    unmeasuredElement::UnmeasredElementRule
    coating::Coating
end



struct NoUnmeasuredElement

end

compute(noume::NoUnmeasuredElement, comp::Material)::Material =
    comp






function iterate(iter::Iteration, tolerance = 1.0e-5)::Material
    first = firstEstimate(kratios)




end
