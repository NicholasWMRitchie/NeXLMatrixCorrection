using Weave
using NeXLMatrixCorrection

let start_dir = pwd()
    cd(@__DIR__)
    outpath = normpath(joinpath(@__DIR__, "..", "docs", "src"))
    @show outpath
    if !isdirpath(outpath)
        mkpath(outpath)
    end

    weave("example.jmd", out_path=joinpath(outpath,"example.md"), doctype="github")
    weave("coatingthickness.jmd", out_path=joinpath(outpath,"coatingthickness.md"), doctype="github")
    weave("testagainstpap.jmd", out_path=joinpath(outpath,"testagainstpap.md"), doctype="github")
    weave("testagainstheinrich.jmd", out_path=joinpath(outpath,"testagainstheinrich.md"), doctype="github")

    cd(start_dir)
end
