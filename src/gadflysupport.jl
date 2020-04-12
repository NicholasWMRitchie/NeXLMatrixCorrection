using .Gadfly

using NeXLCore

function Gadfly.plot(rur::RecordingUpdateRule, measured::KRatio, celm=nothing)
    elm = measured.element
    ks, cs, lbls = [], [], []
    for i in eachindex(rur.comps)
        estkr = rur.comps[i][elm]*rur.zafs[i][elm]/rur.meas[elm].standard[elm]
        push!(ks,estkr)
        push!(cs,rur.comps[i][elm])
        push!(lbls,"$(i)")
    end
    if isnothing(celm)
        celm = rur.comps[end][elm]
    end
    l1 = layer(x=ks,y=cs,label=lbls, Geom.point, Geom.label(position=:above))
    l2 = layer(x=[value(measured.kratio)], y=[celm], Geom.point,  Theme(default_color="red"))
    plot(l2, l1, Guide.XLabel("k[$(elm.symbol)]"), Guide.YLabel("C[$(elm.symbol)]"))
end

function Gadfly.plot(rur::RecordingUpdateRule, measured::Vector{KRatio}, celm::Union{Nothing,Material}=nothing, width=3)
    plots =Union{Plot,Gadfly.Context}[plot(rur, kr, isnothing(celm) ? celm : value(celm[kr.element])) for kr in measured]
    height = convert(Int,floor((length(plots)+width-1)/width))
    for i in length(plots)+1:height*width
        push!(plots,Gadfly.context())
    end
    gridstack(reshape(plots, height, width))
end

function Gadfly.plot(tmc::Type{<:MatrixCorrection}, mat::Material, shells::Vector{AtomicSubShell}, beamEnergy::Float64)
    zz, sh, prz = Float64[], String[], Float64[]
    r = range(tmc, mat, beamEnergy)
    for shell in shells
        mc = matrixcorrection(tmc, mat, shell, beamEnergy)
        for z in range(0.0, stop = r, length = 100)
            push!(zz, z)
            push!(sh, repr(shell))
            push!(prz, ϕ(mc, z))
        end
    end
    df = DataFrame(ρz=zz, Shell=sh, ϕρz=prz)
    plot(df, x=:ρz, y=:ϕρz, color=:Shell, Geom.line,
        Coord.Cartesian(xmin=0.0, xmax=r), Guide.xlabel("ρz [g/cm²]"), Guide.ylabel("ϕ(ρz)"))
end

function Gadfly.plot(tmc::Type{<:MatrixCorrection}, mat::Material, cxrs::Vector{CharXRay}, beamEnergy::Float64, takeOffAngle)
    zz, sh, prz, lsty = Float64[], String[], Float64[], Int[]
    r = range(tmc, mat, beamEnergy)
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
    df = DataFrame(ρz=zz, Line=sh, LineStyle=lsty, ϕρz=prz)
    plot(df, x=:ρz, y=:ϕρz, color=:Line, linestyle=:LineStyle, Geom.line, Scale.linestyle_discrete(),
        Coord.Cartesian(xmin=0.0, xmax=r), Guide.xlabel("ρz [g/cm²]"), Guide.ylabel("ϕ(ρz)"),
        Guide.title("$(name(mat)) at $(0.001*beamEnergy) keV"))
end
