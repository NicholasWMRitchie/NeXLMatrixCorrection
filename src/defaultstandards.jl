using CSV
const __default_standards = Set{Material}()

function getstandards(elm::Element, cMin=0.01)
  if isempty(__default_standards)
    for row in CSV.File(joinpath(@__DIR__, "standard.csv"), header=1)
      try
          mat = parse(Material, row.Formula, name=row.Name, description=row.Comment)
          push!(__default_standards, mat)
      catch ex
          @info "getstandards(...)" exception=(ex, catch_backtrace())
          @error "Unable to parse $(row.Formula) in standard.csv as a Material."
      end
    end
  end
  filter(m -> value(m[elm]) > cMin, __default_standards)
end
