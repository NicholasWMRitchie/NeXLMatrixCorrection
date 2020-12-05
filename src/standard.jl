"""
`Standard` represents a mechanism for combining `KRatio`s measured from a `Material`
so that they can serve as similar standards in the quantification process.
The `KRatio`s are typically measured from a complex material that is chosen to be similar 
in composition to an unknown.  The complex material and the unknown are measured using
the same simple pure elemental or simple compound standard.  For example, you 
might construct a `Standard` using the material Ca₅(PO₄)₃F as a complex standard for Ca.  
The k-ratios in the standard might be from Ca₅(PO₄)₃F) measured using CaF₂ on the Ca Kα 
and Kβ lines. This Standard might then be used to quantify a natural Apatite sample which
was also measured using CaF₂.  
"""
struct Standard
    material::Material # The complex material to be used as a standard
    kratios::Vector{KRatio} # k-ratios measured from `material`

    function Standard(mat::Material, krs::AbstractVector{KRatio})
        @assert length(krs) > 0 "Specify at least one k-ratio for $mat."
        @assert krs[1].element in keys(mat) "The element $(krs[1].element) must exist in $mat."
        @assert all(kr->kr.element == krs[1].element, krs) "Not all the k-ratios are for the element $elm."
        return new(mat, krs)
    end
end

NeXLCore.element(std::Standard) = std.kratios[1].element
NeXLCore.material(std::Standard) = std.material
function Base.show(io::IO, std::Standard) 
    println(io, "Standard[$(material(std)) for $(element(std)),")
    for kr in std.kratios
        println(io,"\t$kr")
    end
    println(io, "]")
end

"""
    findmatch(kr::KRatio, std::Standard)::Boolean

Find the `KRatio` within the specified `Standard` that is appropriate to standardize `kr`.
The two `KRatio` objects must be collected on the same element using the same characteristic X-ray lines, beam 
energy and take-off angle.
"""
function findmatch(kr::KRatio, std::Standard)::Union{KRatio,Nothing}
    if element(kr) == element(std) 
        for skr in std.kratios
            if kr.stdProps[:BeamEnergy] == skr.stdProps[:BeamEnergy] && #
                kr.stdProps[:TakeOffAngle] == skr.stdProps[:TakeOffAngle] && #
                isequal(kr.standard, skr.standard) && kr.lines == skr.lines
                return skr
            end
        end    
    end
    return nothing
end

"""
    standardize(kr::KRatio, standards::AbstractVector{Material, KRatio})::KRatio

Used to convert measurements relative to a simple standard (`kr`) to measurements relative to a similar standards (`stds`).
Whenever possible, the KRatio `kr` is replaced by the ratio of `kr/std` where `std` is a KRatio that matches `kr`
in the sense that it was measured using the same simple standard, same characteristic lines, the same beam energy and take-off 
angle.
"""
function standardize(kr::KRatio, standards::AbstractVector{Standard})::KRatio
    function divide(n, d)
        f = value(n)/value(d)
        s = sqrt(f^2*((σ(n)/value(n))^2+(σ(d)/value(d))^2))
        return s > 0 ? UncertainValue(f, s) : f
    end
    for standard in standards
        stdk = findmatch(kr, standard)
        if !isnothing(stdk)
            @assert kr.element==stdk.element "The elements in $kr & $stdk don't match."
            @assert kr.lines == stdk.lines "The characteristic lines in $kr & $stdk don't match."
            @assert kr.stdProps[:BeamEnergy] == stdk.stdProps[:BeamEnergy] "The beam energies in $kr & $stdk don't match."
            @assert kr.stdProps[:TakeOffAngle] == stdk.stdProps[:TakeOffAngle] "The take-off angles in $kr & $stdk don't match."
            return KRatio(kr.lines, kr.unkProps, stdk.unkProps, standard.material, divide(kr.kratio,stdk.kratio))
        end
    end
    return kr
end