using DataFrames

function zaf(
        mat::Material,
        e0::Float64,
        toa::Float64=deg2rad(40.0);
        mc::Type{<:MatrixCorrection} = XPP,
        fc::Type{<:FluorescenceCorrection}=ReedFluorescence,
        coatings::Dict{String,Film}=Dict{String,Film}(),
        stds::Dict{Element,Material}=Dict{Element,Material}(),
        minweight=0.01)
    # @info "Using the $mc Z⋅A correction and the $fc F correction."
    stdfor(elm) = haskey(stds, elm) ? stds[elm] : pure(elm)
    coatz(mat) = get(coatings,name(mat),missing)
    @assert toa>=0.0 && toa <= π/2 "The take-off angle is out-of-range [0, π/2]"
    @assert e0 > 1.0e3 "The beam energy less than 1 keV."
    res = DataFrame(Material=String[], Standard=String[], Line=CharXRay[], weight=Float64[], E₀=Float64[], Z=Float64[], A=Float64[], #
                    F=Float64[], c=Float64[], ZAFc=Float64[], k=Float64[])
    for elm in keys(mat)
        cxrs = characteristic(elm, alltransitions, 0.01, e0)
        std = stdfor(elm)
        for (sh, cxrs) in splitbyshell(cxrs)
            zaf = zafcorrection(mc, fc, Coating, mat, std, sh, e0, unkCoating=coatz(mat), stdCoating=coatz(std))
            for cxr in filter(cxr->weight(cxr)>=minweight, cxrs)
                row = [
                    mat.name, std.name, cxr, weight(cxr), e0, #
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

function NeXLMatrixCorrection.zaf(ir::IterationResult; minweight=0.01)
    stds = Dict( kr.element=>kr.standard for kr in ir.kratios)
    spec = ir.label.spectrum
    NeXLMatrixCorrection.zaf(ir.comp, spec[:BeamEnergy], spec[:TakeOffAngle];
        mc=ir.iterate.mctype, fc=ir.iterate.fctype, stds=stds, minweight=minweight)
end

function NeXLUncertainties.asa(
    ::Type{DataFrame}, 
    krs::AbstractVector{KRatio}, 
    withComputedKs::Bool, 
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection}=ReedFluorescence
)
    fm, sm, cm = Union{Float64,Missing}, Union{String,Missing}, Union{Material,Missing}
    lines, mease0, meastoa, meascomp = String[], fm[], fm[], cm[]
    refe0, reftoa, refcomp, krv, dkrv, cks, ratio = fm[], fm[], cm[], Float64[], Float64[], fm[], fm[]
    for kr in krs
        push!(lines, repr(kr.lines))
        meas = kr.unkProps
        push!(mease0, get(meas, :BeamEnergy, missing))
        push!(meastoa, get(meas, :TakeOffAngle, missing))
        push!(meascomp, get(meas, :Composition, missing))
        ref = kr.stdProps
        push!(refe0, get(ref, :BeamEnergy, missing))
        push!(reftoa, get(ref, :TakeOffAngle, missing))
        push!(refcomp, get(ref, :Composition, missing))
        push!(krv, value(kr.kratio))
        push!(dkrv, σ(kr.kratio))
        if withComputedKs
            elm = element(kr.lines[1])
            prim = brightest(kr.lines)
            if any(ismissing.((meascomp[end], mease0[end], meastoa[end], refcomp[end], refe0[end], reftoa[end]))) ||
               (NeXLCore.nonneg(meascomp[end], elm) < 1.0e-6) ||
               (NeXLCore.nonneg(refcomp[end], elm) < 1.0e-6) || (energy(inner(prim))>0.95*min(mease0[end],refe0[end]))
                push!(cks, missing)
                push!(ratio, missing)
            else
                zs = zafcorrection(mc, fc, Coating, meascomp[end], kr.lines, mease0[end])
                zr = zafcorrection(mc, fc, Coating, refcomp[end], kr.lines, refe0[end])
                k =
                    gZAFc(zs, zr, meastoa[end], reftoa[end]) * NeXLCore.nonneg(meascomp[end], elm) /
                    NeXLCore.nonneg(refcomp[end], elm)
                push!(cks, k)
                push!(ratio, value(kr.kratio) / k)
            end
        end
    end
    res = DataFrame(
        Lines = lines,
        E0meas = mease0,
        TOAmeas = meastoa,
        Cmeas = meascomp,
        E0ref = refe0,
        TOAref = reftoa,
        Cref = refcomp,
        K = krv,
        ΔK = dkrv,
    )
    if withComputedKs
        res[:, :Kxpp] = cks
        res[:, :Ratio] = ratio
    end
    return res
end