struct Vertex
  position::Point
  value::Point
end
Base.isequal(v::Vertex, p::Point) = isequal(v.position, p)
function periodicindex(iterable, i::Int)
    l = Base.length(iterable)
    i < 1 && return l + i
    i > l && return i - l
    return i
end

