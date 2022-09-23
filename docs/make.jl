using Documenter
using Weave
using NeXLMatrixCorrection

include("../weave/buildweave.jl")

pages = [
            "Example" => "example.md",
            "k-ratio Round Trip" => "roundtrip.md",
            "Advanced Iteration" => "advanced.md",
            "Coating Thickness" => "coatingthickness.md",
            "Compare To PAP" => "testagainstpap.md",
            "Compare to Heinrich" => "testagainstheinrich.md",
            "Emitted Intensities" => "emitted.md"
        ]

makedocs(modules = [NeXLMatrixCorrection], sitename = "NeXLMatrixCorrection.jl", pages = pages)

function addNISTHeaders(htmlfile::String)
    # read HTML
    html = transcode(String,read(htmlfile))
    # Find </head>
    i = findfirst(r"</[Hh][Ee][Aa][Dd]>", html)
    # Already added???
    j = findfirst("nist-header-footer", html)

    if isnothing(j) && (!isnothing(i))
        # Load header text
        header = open(joinpath(@__DIR__, "src", "assets", "nist_header.txt")) do io
            join(readlines(io),"\n")
        end
        # Insert nist-pages links right before </head>
        res = html[1:i.start-1]*"\n"*header*"\n"*html[i.start:end]
        write(htmlfile, res)
        println("Inserting NIST header/footer into $htmlfile")
    end
    return htmlfile
end

addNISTHeaders(joinpath(@__DIR__, "build","index.html"))
addNISTHeaders.(map(name->joinpath(@__DIR__, "build", splitext(name)[1], "index.html"), map(p->p.second, pages)))

# deploydocs(repo = "github.com/NeXLMatrixCorrection.jl.git")
