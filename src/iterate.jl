"""
    KRatio

The k-ratio is the result of two intensity measurements - one on a standard
with known composition and one on an unknown. Each measurement has properties
like :BeamEnergy, :TakeOffAngle, :Coating that characterize the measurement.

Properties: (These Symbols are intentionally the same used in NeXLSpectrum)

    :BeamEnergy incident beam energy in eV
    :TakeOffAngle in radians
    :Coating A NeXLCore.Layer object describing a conductive coating
"""
struct KRatio
    lines::Vector{CharXRay} # Which CharXRays were measured?
    unkProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    stdProps::Dict{Symbol,Any} # Beam energy, take-off angle, coating, ???
    standard::Material
    kratio::AbstractFloat
    function KRatio(
        lines::Vector{CharXRay},
        unkProps::Dict{Symbol,Any},
        stdProps::Dict{Symbol,Any},
        standard::Material,
        kratio::AbstractFloat
        )
            elm = element(lines[0])
            if length(lines)==0
                error("Must specify at least one characteristic X-ray.")
            if !all(l->element(l)==elm for l in lines)
                error("The characteristic X-rays must all be from the same element.")
            if standard[elm]==0.0
                error("The standard must contain the element $(elm).")
            if kratio<0.0
                error("The k-ratio must be non-negative.")
            return new(lines,unkProps,stdProps,standard,kratio)
        end
end

abstract type UnmeasuredElementRule end

"""
    compute(::Type{UnmeasuredElementRule}, inp::Dict{Element,AbstractFloat})::Dict{Element,AbstractFloat}

A null UnmeasuredElementRule.  Just returns the inputs.
"""
compute(::Type{<:UnmeasuredElementRule}, inp::Dict{Element,AbstractFloat})::Dict{Element,AbstractFloat} =
    inp

struct Iteration
    name::String # Sample name
    mctype::Type{<:MatrixCorrection}
    fctype::Type{<:FluorescenceCorrection}
    kratios::Vector{KRatio}
    zaf::Dict{KRatio, MultiZAF}
    unmeasuredElement::UnmeasuredElementRule
end

allElements(mat1::Material, mat2::Material) =
    union(keys(mat1),keys(mat2))

delta(mat1::Material, mat2::Material)::Dict{Element,AbstractFloat} =
    Dict((elm, mat1[elm] - mat2[elm]) for elm in allElements(mat1,mat2))

test1(mat1::Material, mat2::Material, tol::Float64)::Bool =
    sqrt(mapreduce(elm->(mat1[elm]-mat2[elm])^2,+,allElements(mat1,mat2)))<tol

test2(mat1::Material, mat2::Material, tol::Float64)::Bool =
    all(elm->mat1[elm]-mat2[elm]<tol for elm in allElements(mat1,mat2))

function firstEstimate(iter::Iteration)::Material
    mfs = Dict{Element,Float64}()
    for kr in iter.kratios
        elm = element(kr.lines[0])
        mfs[elm] = kr.kratio * kr.standard[elm]
    end
    return material(iter.name, mfs)
end

function computeKs(iter::Iteration, est::Material):Dict{Element, AbstractFloat}
    for kr in iter.kratios
        e0 = kr.stdProps[:BeamEnergy]
        coating = get(kr.stdProps, :Coating, NeXLCore.NullCoating())
        unkZaf = ZAF(iter.mctype,iter.fctype,kr.lines,e0,coating)
        
        k = (Cunk * unkZaf)/(Cstd * stdZaf)

        kr.stdZaf

        ZAF(
            mctype::Type{<:MatrixCorrection},
            fctype::Type{<:FluorescenceCorrection},
            mat::Material,
            cxrs,
            e0::AbstractFloat,
            coating = NeXLCore.NullCoating()
        )

    end


end





function iterate(iter::Iteration, tolerance = 1.0e-5)::Material
    first = firstEstimate(kratios)




end
