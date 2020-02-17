using .Gadfly

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
