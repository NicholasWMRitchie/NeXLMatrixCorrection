const __default_standards = Set{Material}()

function getstandards(elm::Element, cMin=0.01)
  if isempty(__default_standards)
    for line in eachline(joinpath(@__DIR__, "standards.txt"))
      try
          mat = parse(Material, line)
          push!(__default_standards, mat)
      catch
          @error "Unable to parse $line in standards.txt as a Material."
      end
    end
  end
  filter(m -> m[elm] > cMin, __default_standards)
end
