module NeXLMatrixCorrectionGadflyExt

using NeXLMatrixCorrection
using Gadfly

using NeXLCore
using Colors

export plot2

function Gadfly.plot(tmc::Type{<:MatrixCorrection}, mat::Material, shells::Vector{AtomicSubShell}, beamEnergy::Float64)
    zz, sh, prz = Float64[], String[], Float64[]
    r = range(tmc, mat, beamEnergy, false)
    for shell in shells
        mc = matrixcorrection(tmc, mat, shell, beamEnergy)
        for z in range(0.0, stop = r, length = 100)
            push!(zz, z)
            push!(sh, repr(shell))
            push!(prz, ϕ(mc, z))
        end
    end
    df = DataFrame(ρz = zz, Shell = sh, ϕρz = prz)
    plot(
        df,
        x = :ρz,
        y = :ϕρz,
        color = :Shell,
        Geom.line,
        Coord.Cartesian(xmin = 0.0, xmax = r),
        Guide.xlabel("ρz [g/cm²]"),
        Guide.ylabel("ϕ(ρz)"),
    )
end

function Gadfly.plot(
    tmc::Type{<:MatrixCorrection},
    mat::Material,
    cxrs::Vector{CharXRay},
    beamEnergy::AbstractFloat,
    takeOffAngle,
)
    zz, sh, prz, lsty = Float64[], String[], Float64[], Int[]
    r = range(tmc, mat, beamEnergy, false)
    for cxr in cxrs
        shell = inner(cxr)
        mc = matrixcorrection(tmc, mat, shell, beamEnergy)
        for z in range(0.0, stop = r, length = 100)
            push!(zz, z)
            push!(sh, "$cxr")
            push!(prz, ϕ(mc, z))
            push!(lsty, 1)
            push!(zz, z)
            push!(sh, "$cxr")
            push!(prz, ϕabs(mc, z, cxr, takeOffAngle))
            push!(lsty, 2)
        end
    end
    df = DataFrame(ρz = zz, Line = sh, LineStyle = lsty, ϕρz = prz)
    plot(
        df,
        x = :ρz,
        y = :ϕρz,
        color = :Line,
        linestyle = :LineStyle,
        Geom.line,
        Scale.linestyle_discrete(),
        Coord.Cartesian(xmin = 0.0, xmax = r),
        Guide.xlabel("ρz [g/cm²]"),
        Guide.ylabel("ϕ(ρz)"),
        Guide.title("$(name(mat)) at $(0.001*beamEnergy) keV"),
    )
end

function Gadfly.plot(
    tmcs::AbstractVector{DataType},
    mat::Material,
    cxr::CharXRay,
    beamEnergy::AbstractFloat,
    takeOffAngle::AbstractFloat=deg2rad(40.0),
)
    zz, sh, prz, lsty = Float64[], String[], Float64[], Int[]
    r = 1.2*(range(Kanaya1972, mat, beamEnergy, false) - range(Kanaya1972, mat, energy(inner(cxr)), false))
    shell = inner(cxr)
    for z in range(0.0, stop = r, length = 100)
        for tmc in tmcs
            mc = matrixcorrection(tmc, mat, shell, beamEnergy)
            push!(zz, z)
            push!(sh, repr(tmc, context=:compact=>true)*":gen")
            push!(prz, ϕ(mc, z))
            push!(lsty, 1)
            push!(zz, z)
            push!(sh, repr(tmc, context=:compact=>true)*":abs")
            push!(prz, ϕabs(mc, z, cxr, takeOffAngle))
            push!(lsty, 2)
        end
    end
    df = DataFrame(ρz = zz, Line = sh, LineStyle = lsty, ϕρz = prz)
    plot(
        df,
        x = :ρz,
        y = :ϕρz,
        color = :Line,
        linestyle = :LineStyle,
        Geom.line,
        Scale.linestyle_discrete(),
        Coord.Cartesian(xmin = 0.0, xmax = r),
        Guide.xlabel("ρz [g/cm²]"),
        Guide.ylabel("ϕ(ρz)"),
        Guide.title("$(repr(cxr)) in $(name(mat)) at $(0.001*beamEnergy) keV"),
    )
end


function Gadfly.plot(krs::AbstractArray{KRatio}, mc=XPP, fc=ReedFluorescence, cc=Coating; palette=NeXLPalette)
    mfs, kok, dkok, color = String[], Float64[], Float64[], Color[]
    next = 1
    matcolors = Dict{String,RGB{Float64}}()
    for kr in filter(kr->haskey(kr.unkProps, :Composition), krs)
        unkComp = kr.unkProps[:Composition]
        if hasminrequired(mc, kr.unkProps) && hasminrequired(fc, kr.unkProps) && #
           hasminrequired(mc, kr.stdProps) && hasminrequired(fc, kr.stdProps) && #
           (value(unkComp[kr.element]) > 0.0) && (value(kr.standard[kr.element]) > 0.0)
             # Compute the k-ratio
             kc = gZAFc(kr, unkComp, mc, fc, cc) * (value(unkComp[kr.element]) / value(kr.standard[kr.element]))
             push!(mfs, name(shell(brightest(kr.xrays)))) # value(unkComp[kr.element]))
             push!(kok, value(kr.kratio) / kc)
             push!(dkok, σ(kr.kratio) / kc)
             matname = "C($(symbol(kr.element)), $(name(unkComp)), $(kr.unkProps[:BeamEnergy]/1000.0) keV, $(name(kr.stdProps[:Composition])))"
             if !haskey(matcolors, matname)
                 matcolors[matname] = palette[next]
                 next += 1
             end
             push!(color, matcolors[matname])
        end
    end
    plot(
        x = mfs,
        y = kok,
        ymin = kok .- dkok,
        ymax = kok .+ dkok,
        color = color,
        Geom.errorbar,
        Stat.x_jitter(range = 0.8),
        Guide.manual_color_key("Material", [keys(matcolors)...], [values(matcolors)...]), # Guide.yrug,
        Guide.xlabel("Shell"),
        Guide.ylabel("k[Measured]/k[Calculated]"),
    )
end

function Gadfly.plot(krs::AbstractArray{KRatio}, unkComp::Material, mc=XPP, fc=ReedFluorescence, cc=Coating; palette=NeXLPalette)
    mfs, kok, dkok, color = String[], Float64[], Float64[], Color[]
    next = 1
    matcolors = Dict{String,RGB{Float64}}()
    for kr in filter(kr->name(get(kr.unkProps, :Composition, unkComp)) == name(unkComp), krs)
        if hasminrequired(mc, kr.unkProps) &&
           hasminrequired(fc, kr.unkProps) && #
           hasminrequired(mc, kr.stdProps) &&
           hasminrequired(fc, kr.stdProps) && #
           (value(unkComp[kr.element]) > 0.0) &&
           (value(kr.standard[kr.element]) > 0.0)
            # Compute the k-ratio
            kc = gZAFc(kr, unkComp, mc, fc, cc) * (value(unkComp[kr.element]) / value(kr.standard[kr.element]))
            push!(mfs, name(shell(brightest(kr.xrays)))) # value(unkComp[kr.element]))
            push!(kok, value(kr.kratio) / kc)
            push!(dkok, σ(kr.kratio) / kc)
            matname = "$(name(unkComp)) $(kr.unkProps[:BeamEnergy]/1000.0) keV"
            if !haskey(matcolors, matname)
                matcolors[matname] = palette[next]
                next += 1
            end
            push!(color, matcolors[matname])
        end
    end
    plot(
        x = mfs,
        y = kok,
        ymin = kok .- dkok,
        ymax = kok .+ dkok,
        color = color,
        Geom.errorbar,
        Stat.x_jitter(range = 0.4),
        Guide.manual_color_key("Material", [keys(matcolors)...], [values(matcolors)...]), # Guide.yrug,
        Guide.xlabel("Shell"),
        Guide.ylabel("k[Measured]/k[Calculated]"),
    )
end

Gadfly.plot(irs::AbstractArray{IterationResult}; known::Union{Material,Missing} = missing, delta::Bool = false) =
    plot([ir.comp for ir in irs], known = known, delta = delta, label = "Measurement")

plot2(irs::AbstractVector{IterationResult}; known::Union{Material, Missing}=missing) =
    NeXLCore.plot2([ir.comp for ir in irs], known=known, label="Measurement")

end # module