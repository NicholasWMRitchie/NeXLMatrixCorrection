using Weave

let start_dir = pwd()
    cd(@__DIR__)
    outpath = normpath(joinpath(@__DIR__, "..", "docs", "build"))
    @show outpath
    if !isdirpath(outpath)
        mkpath(outpath)
    end

    weave("coatingthickness.jmd", out_path=joinpath(outpath,"coatingthickness.html"))
    weave("example.jmd", out_path=joinpath(outpath,"example.html"))

    cd(start_dir)
end
