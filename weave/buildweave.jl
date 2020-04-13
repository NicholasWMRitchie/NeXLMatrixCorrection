let start_dir = pwd()
    cd(@__DIR__)

    weave("coatingthickness.jmd", out_path="../docs/src/coatingthickness.html")

    cd(start_dir)
end
