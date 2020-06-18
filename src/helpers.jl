using DataFrames

function zaf(
        mat::Material,
        e0::Float64,
        toa::Float64=deg2rad(40.0);
        mc::Type{<:MatrixCorrection} = XPP,
        fc::Type{<:FluorescenceCorrection}=ReedFluorescence,
        coatings::Dict{String,Film}=Dict{String,Film}(),
        stds::Dict{Element,Material}=Dict{Element,Material}())
    @info "Using the $mc Z⋅A correction and the $fc F correction."
    stdfor(elm) = haskey(stds, elm) ? stds[elm] : pure(elm)
    coatz(mat) = get(coatings,name(mat),missing)
    @assert toa>=0.0 && toa <= π/2 "The take-off angle is out-of-range [0, π/2]"
    @assert e0 > 1.0e3 "The beam energy less than 1 keV."
    res = DataFrame(Material=String[], Standard=String[], Line=CharXRay[], E₀=Float64[], Z=Float64[], A=Float64[], #
                    F=Float64[], c=Float64[], ZAFc=Float64[], k=Float64[])
    for elm in keys(mat)
        cxrs = characteristic(elm, alltransitions, 0.01, e0)
        std = stdfor(elm)
        for (sh, cxrs) in splitbyshell(cxrs)
            zaf = zafcorrection(mc, fc, Coating, mat, std, sh, e0, unkCoating=coatz(mat), stdCoating=coatz(std))
            for cxr in cxrs
                row = [
                    mat.name, std.name, cxr, e0, #
                    Z(zaf...),
                    A(zaf..., cxr, toa, toa),
                    F(zaf..., cxr, toa, toa),
                    NeXLMatrixCorrection.coating(zaf..., cxr, toa, toa),
                    ZAFc(zaf..., cxr, toa, toa),
                    k(zaf..., cxr, toa, toa)
                ]
                push!(res, row)
            end
        end
    end
    return res
end
