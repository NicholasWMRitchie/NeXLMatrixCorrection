module NeXLMatrixCorrectionDataFramesExt

using NeXLCore
using NeXLMatrixCorrection
using DataFrames

NeXLCore.compare(itress::AbstractVector{IterationResult}, known::Material)::DataFrame =
    mapreduce(itres -> compare(itres, known), append!, itress)
NeXLCore.compare(itres::IterationResult, known::Material)::DataFrame = compare(itres.comp, known)
NeXLCore.compare(iter1::IterationResult, iter2::IterationResult) = compare(iter1.comp, iter2.comp)


function NeXLMatrixCorrection.zaf(
        mat::Material,
        e0::AbstractFloat,
        toa::AbstractFloat=deg2rad(40.0);
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
            for cxr in filter(cxr->weight(NormalizeToUnity, cxr)>=minweight, cxrs)
                row = [
                    mat.name, std.name, cxr, weight(NormalizeToUnity, cxr), e0, #
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

function DataFrames.DataFrame(
    krs::AbstractVector{KRatio}, 
    withComputedKs::Bool, 
    mc::Type{<:MatrixCorrection} = XPP,
    fc::Type{<:FluorescenceCorrection}=ReedFluorescence
)
    fm, cm = Union{Float64,Missing}, Union{Material,Missing}
    xrays, mease0, meastoa, meascomp = String[], fm[], fm[], cm[]
    refe0, reftoa, refcomp, krv, dkrv, cks, ratio = fm[], fm[], cm[], Float64[], Float64[], fm[], fm[]
    for kr in krs
        push!(xrays, repr(kr.xrays))
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
            elm = element(kr.xrays[1])
            prim = brightest(kr.xrays)
            if any(ismissing.((meascomp[end], mease0[end], meastoa[end], refcomp[end], refe0[end], reftoa[end]))) ||
               (NeXLCore.nonneg(meascomp[end], elm) < 1.0e-6) ||
               (NeXLCore.nonneg(refcomp[end], elm) < 1.0e-6) || (energy(inner(prim))>0.95*min(mease0[end],refe0[end]))
                push!(cks, missing)
                push!(ratio, missing)
            else
                zs = zafcorrection(mc, fc, Coating, meascomp[end], kr.xrays, mease0[end])
                zr = zafcorrection(mc, fc, Coating, refcomp[end], kr.xrays, refe0[end])
                k =
                    gZAFc(zs, zr, meastoa[end], reftoa[end]) * NeXLCore.nonneg(meascomp[end], elm) /
                    NeXLCore.nonneg(refcomp[end], elm)
                push!(cks, k)
                push!(ratio, value(kr.kratio) / k)
            end
        end
    end
    res = DataFrame(
        Xrays = xrays,
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

function DataFrames.DataFrame(ir::IterationResult; withZAF::Bool=true)
    elms, mfs, ks, cks, labels = String[], Float64[], Union{Missing,Float64}[], Union{Missing,Float64}[], String[]
    g, z, a, f = Union{Missing,Float64}[], Union{Missing,Float64}[], Union{Missing,Float64}[], Union{Missing,Float64}[]
    c, gzafc, stds, dmfs, xrays = Union{Missing,Float64}[], Union{Missing,Float64}[], String[], Float64[], String[]
    for elm in keys(ir.comp)
        push!(labels, repr(ir.label))
        push!(elms, elm.symbol)
        rc = round(ir.comp[elm])
        push!(mfs, value(rc))
        push!(dmfs, σ(rc))
        added = false  # computed element
        for kr in ir.kratios
            if elm == kr.element
                push!(ks, value(kr.kratio))
                push!(xrays, repr(kr.xrays))
                if withZAF
                    zafs = _ZAF(ir.iterate, convert(Material{Float64,Float64}, kr.standard), kr.stdProps, kr.xrays)
                    zafu = _ZAF(ir.iterate, convert(Material{Float64,Float64}, ir.comp), kr.unkProps, kr.xrays)
                    push!(cks, σ(kr.kratio))
                    push!(stds, name(kr.standard))
                    push!(c, coating(zafu, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle]))
                    push!(z, Z(zafu, zafs))
                    push!(a, A(zafu, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle]))
                    push!(f, F(zafu, zafs, kr.unkProps[:TakeOffAngle], kr.stdProps[:TakeOffAngle]))
                    push!(g, generation(zafu, zafs))
                    push!(gzafc, g[end] * z[end] * a[end] * f[end] * c[end])
                else
                    push!(cks, ir.computed[elm])
                end
                added = true
                break
            end
        end
        if !added
            push!(ks, missing), push!(cks, missing), push!(stds, missing), push!(c, missing), push!(z, missing),
            push!(a, missing), push!(f, missing), push!(g, missing), push!(gzafc, missing), push!(xrays, missing)
        end
    end
    return withZAF ?
           DataFrame(
        :Label => labels,
        :Element => elms,
        :Standard => stds,
        :Xrays => xrays,
        #:Converged => [ir.converged for elm in keys(ir.comp)],
        #:Iterations => [ir.iterations for elm in keys(ir.comp)],
        Symbol("C") => mfs,
        Symbol("ΔC") => dmfs,
        Symbol("k") => ks,
        Symbol("Δk") => cks,
        :g => g,
        :Z => z,
        :A => a,
        :F => f,
        :c => c,
        :gZAFc => gzafc,
    ) :
           DataFrame(
        :Label => labels,
        :Element => elms,
        :Converged => [ir.converged for elm in keys(ir.comp)],
        :Iterations => [ir.iterations for elm in keys(ir.comp)],
        Symbol("C") => mfs,
        Symbol("ΔC") => dmfs,
        Symbol("k") => ks,
        Symbol("Δk") => cks,
    )
    return res
end

function DataFrames.DataFrame(irs::AbstractVector{IterationResult}; mode=:MassFraction, nominal=nothing)::DataFrame
    if isnothing(nominal)
        DataFrame(map(ir -> ir.comp, irs), mode)
    else
        DataFrame([(ir.comp for ir in irs)..., nominal], mode)
    end
end


DataFrames.describe(irs::AbstractVector{IterationResult}) =
    describe(DataFrame(irs)[:, 2:end], :mean, :std, :min, :q25, :median, :q75, :max)

    """
    DataFrame(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard in a DataFrame.
"""
function DataFrames.DataFrame(
    unk::MultiZAF,
    std::MultiZAF,
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame
    tot = gZAFc(unk, std, θunk, θstd)
    @assert isequal(element(unk), element(std)) "The unknown and standard's elements must match."
    return DataFrame(
        Unknown=[name(material(unk))],
        E₀ᵤ=[beamEnergy(unk)],
        Standard=[name(material(std))],
        E₀ₛ=[beamEnergy(std)],
        Xrays=[name(union(characteristic(unk), characteristic(std)))],
        Generation=[generation(unk, std)],
        Z=[Z(unk, std)],
        A=[A(unk, std, θunk, θstd)],
        F=[F(unk, std, θunk, θstd)],
        coating=[coating(unk, std, θunk, θstd)],
        gZAFc=[tot],
        k=[tot * material(unk)[element(unk)] / material(std)[element(unk)]],
    )
end

"""
    detail(unk::MultiZAF, std::MultiZAF)::DataFrame

Tabulate each term in the MultiZAF matrix correction in a DataFrame.
"""
function NeXLMatrixCorrection.detail(unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame
    stds, stdE0, unks = String[], Float64[], String[]
    unkE0, xray, g = Float64[], CharXRay[], Float64[]
    z, a, f = Float64[], Float64[], Float64[]
    c, wgt, zaf, k = Float64[], Float64[], Float64[], Float64[]
    for (sh, cxrs2) in splitbyshell(commonXrays(unk, std))
        zafU, zafS = unk.zafs[sh], std.zafs[sh]
        matU, matS = material(zafU), material(zafS)
        elm = element(sh)
        for cxr in cxrs2
            push!(unks, name(material(zafU)))
            push!(unkE0, beamEnergy(zafU))
            push!(stds, name(material(zafS)))
            push!(stdE0, beamEnergy(zafS))
            push!(xray, cxr)
            push!(wgt, weight(NormalizeByShell, cxr))
            push!(g, generation(zafU, zafS, inner(cxr)))
            push!(z, Z(zafU, zafS))
            push!(a, A(zafU, zafS, cxr, θunk, θstd))
            push!(f, F(zafU, zafS, cxr, θunk, θstd))
            push!(c, coating(zafU, zafS, cxr, θunk, θstd))
            tot = ZAFc(zafU, zafS, cxr, θunk, θstd)
            push!(zaf, tot)
            push!(k, tot * matU[elm] / matS[elm])
        end
    end
    return DataFrame(
        Unknown=unks,
        E0unk=0.001 * unkE0,
        Standard=stds,
        E0std=0.001 * stdE0,
        Xray=xray,
        Weight=wgt,
        Generation=g,
        Z=z,
        A=a,
        F=f,
        c=c,
        ZAF=zaf,
        k=k,
    )
end

"""
    detail(mzs::AbstractArray{Tuple{MultiZAF, MultiZAF}})::DataFrame

Tabulate the details of a matrix correction relative to the specified unknown and standard in a DataFrame.
"""
NeXLMatrixCorrection.detail(mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat) =
    mapreduce(tmm -> detail(tmm[1], tmm[2], θunk, θstd), append!, mzs)

"""
DataFrames.DataFrame(mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}}, θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to a specified Dict of unknowns and standards in a DataFrame.
"""
DataFrames.DataFrame(
    mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}},
    θunk::AbstractFloat,
    θstd::AbstractFloat,
) = mapreduce(tmm -> DataFrame(tmm[1], tmm[2], θunk, θstd), append!, mzs)


function DataFrames.DataFrame(
    unk::Material,
    std::Material,
    xrays::Vector{CharXRay},
    e0::Float64,
    toa::Float64;
    mc::Type{<:MatrixCorrection}=XPP,
    fc::Type{<:FluorescenceCorrection}=ReedFluorescence,
    coating::Type{<:CoatingCorrection}=Coating
)
    flines = filter(cxr -> energy(inner(cxr)) < e0, xrays)
    zafs = zafcorrection(mc, fc, coating, unk, std, flines, e0)
    DataFrames.DataFrame(zafs..., toa, toa)
end

"""
Represents the essential intermediary values for an XPP matrix correction of
characteristic X-rays from a particular atomic sub-shell in a particular material.
"""
function DataFrames.DataFrame(xpps::XPP...)
    return DataFrame(
        SubShell = [xpp.subshell for xpp in xpps],
        Material = [name(xpp.material) for xpp in xpps],
        BeamEnergy = [xpp.E0 for xpp in xpps],
        Phi0 = [xpp.ϕ0 for xpp in xpps],
        A = [xpp.A for xpp in xpps],
        a = [xpp.a for xpp in xpps],
        B = [xpp.B for xpp in xpps],
        b = [xpp.b for xpp in xpps],
        F = [xpp.F for xpp in xpps],
    )
end

"""
    DataFrame(unk::ZAFCorrection, std::ZAFCorrection, trans::AbstractVector{Transition},
    θunk::AbstractFloat, θstd::AbstractFloat)::DataFrame

Tabulate a matrix correction relative to the specified unknown and standard for the iterable of Transition, trans.
"""
function DataFrames.DataFrame(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    trans::NTuple{N,Transition},
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame where {N}
    @assert isequal(atomicsubshell(unk.za), atomicsubshell(std.za))
    "The atomic sub-shell for the standard and unknown don't match."
    cxrs = characteristic(
        element(atomicsubshell(unk.za)),
        trans,
        1.0e-9,
        0.999 * min(beamEnergy(unk.za), beamEnergy(std.za)),
    )
    stds, celm, stdE0, unks, unkE0, xray = String[], Float64[], Float64[], String[], Float64[], CharXRay[]
    z, a, f, c, zaf, k, unkToa, stdToa = #
        Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    for cxr in cxrs
        if isequal(inner(cxr), atomicsubshell(std.za))
            elm = element(cxr)
            push!(unks, name(material(unk.za)))
            push!(celm, material(unk.za)[elm])
            push!(unkE0, beamEnergy(unk.za))
            push!(unkToa, θunk)
            push!(stds, name(material(std.za)))
            push!(stdE0, beamEnergy(std.za))
            push!(stdToa, θstd)
            push!(xray, cxr)
            push!(z, Z(unk, std))
            push!(a, A(unk, std, cxr, θunk, θstd))
            push!(f, F(unk, std, cxr, θunk, θstd))
            push!(c, coating(unk, std, cxr, θunk, θstd))
            tot = ZAFc(unk, std, cxr, θunk, θstd)
            push!(zaf, tot)
            push!(k, tot * material(unk.za)[elm] / material(std.za)[elm])
        end
    end
    return DataFrame(
        Unknown = unks,
        E0unk = unkE0,
        TOAunk = unkToa,
        Standard = stds,
        E0std = stdE0,
        TOAstd = stdToa,
        Xray = xray,
        Z = z,
        A = a,
        F = f,
        c = c,
        ZAF = zaf,
        C = celm,
        k = k,
    )
end

DataFrames.DataFrame(
    unk::ZAFCorrection,
    std::ZAFCorrection,
    θunk::AbstractFloat,
    θstd::AbstractFloat
) = DataFrame(unk, std, alltransitions, θunk, θstd)

function DataFrames.DataFrame(
    unk::Material,
    stds::Dict{Element,Material},
    e0::Real,
    θunk::AbstractFloat,
    θstd::AbstractFloat;
    mctype::Type{<:MatrixCorrection} = XPP,
    fctype::Type{<:FluorescenceCorrection} = ReedFluorescence,
    cctype::Type{<:CoatingCorrection} = Coating,
)
    df = DataFrame()
    for (elm, std) in stds
        for ashell in atomicsubshells(elm, e0)
            append!(
                df,
                DataFrame(
                    zafcorrection(mctype, fctype, cctype, unk, std, ashell, e0)...,
                    alltransitions,
                    θunk,
                    θstd,
                ),
            )
        end
    end
    return df
end

function DataFrames.DataFrame(
    zafcorrs::Dict{ZAFCorrection, ZAFCorrection},
    θunk::AbstractFloat,
    θstd::AbstractFloat,
)::DataFrame
    df = DataFrame()
    for (unk, std) in zafcorrs
        append!(df, DataFrame(unk, std, θunk, θstd))
    end
    return df
end

# Depreciated
NeXLUncertainties.asa(::Type{DataFrame}, krs::AbstractVector{KRatio}, withComputedKs::Bool, mc::Type{<:MatrixCorrection} = XPP, fc::Type{<:FluorescenceCorrection}=ReedFluorescence) = DataFrame(krs, withComputedKs, mc, fc)
NeXLUncertainties.asa(::Type{DataFrame}, ir::IterationResult; kwargs...) = DataFrame(ir; kwargs...)
NeXLUncertainties.asa(::Type{DataFrame}, irs::AbstractVector{IterationResult}; kwargs...) = DataFrames.DataFrame(irs; kwargs...)
NeXLUncertainties.asa(::Type{DataFrame}, unk::MultiZAF, std::MultiZAF, θunk::AbstractFloat, θstd::AbstractFloat) = DataFrame(unk, std, θunk, θstd)
NeXLUncertainties.asa(::Type{DataFrame}, mzs::AbstractArray{Tuple{MultiZAF,MultiZAF}},θunk::AbstractFloat, θstd::AbstractFloat) = DataFrames.DataFrame(mzs ,θunk, θstd)
NeXLUncertainties.asa(::Type{DataFrame}, unk::Material, std::Material, xrays::Vector{CharXRay}, e0::Float64, toa::Float64; kwargs...) = DataFrames.DataFrame(unk, std, xrays, e0, toa; kwargs...)
NeXLUncertainties.asa(::Type{DataFrame}, xpps::XPP...) = DataFrames.DataFrame(xpps::XPP...)
NeXLUncertainties.asa(::Type{DataFrame}, unk::ZAFCorrection, std::ZAFCorrection, trans::NTuple{N,Transition}, θunk::AbstractFloat, θstd::AbstractFloat) where {N} = DataFrames.DataFrame(unk,std,trans,θunk,θstd)
NeXLUncertainties.asa(::Type{DataFrame}, unk::ZAFCorrection, std::ZAFCorrection, θunk::AbstractFloat, θstd::AbstractFloat) = DataFrames.DataFrame( unk, std, θunk, θstd)
NeXLUncertainties.asa(::Type{DataFrame}, zafs::Dict{ZAFCorrection,ZAFCorrection}, θunk::AbstractFloat, θstd::AbstractFloat) = DataFrames.DataFrame(zafs, θunk, θstd)
NeXLUncertainties.asa(::Type{DataFrame}, unk::Material, stds::Dict{Element,Material}, e0::Real, θunk::AbstractFloat, θstd::AbstractFloat; kwargs...) = DataFrames.DataFrame(unk, stds, e0, θunk, θstd; kwargs...)

end