const default_initialnumvertices = 16
const default_minnumvertices = 4

mutable struct Helper
  numvertices::Int
  minnumvertices::Int
  timetaken::Float64
  function Helper(
      initialnumvertices::Int=default_initialnumvertices,
      minnumvertices::Int=default_minnumvertices,
      timetaken=0.0)
    @assert minnumvertices <= initialnumvertices
    return new(initialnumvertices, minnumvertices, timetaken)
  end
end

function indicesonhypercubesurface(lengths::AbstractVector{Int})
  indices = Set{Vector{Int}}()
  for i in CartesianIndices(Tuple(lengths))
    index = collect(Tuple(i))
    (any(index .== 1) || any(index .== lengths)) || continue
    push!(indices, index)
  end
  return indices
end

function sortgriddedpoints(A::Point, B::Point, n::Int)
  @assert n >= 2
  midpoint = middle(A, B)
  indices = indicesonhypercubesurface(n * ones(Int, 2))
  fractions = ((index .- 1) // (n - 1) for index ∈ indices)
  points = [Point(A + (B - A) * fraction) for fraction ∈ fractions]
  sort!(points, lt=(a, b)->isless(angle(a - midpoint), angle(b - midpoint)))
  return points
end

function distributevertices(f::T, A::Point, B::Point, helper::Helper,
                            vertices=Vector{Vertex}()) where {T<:Function}
  number_points_along_edge = 1 + Int(round(helper.numvertices / 4))
  output = Vector{Vertex}()
  for point ∈ sortgriddedpoints(A, B, number_points_along_edge)
    i = findfirst(v->isequal(v, point), vertices)
    isnothing(i) && push!(output, Vertex(point, Point(f(Vector(point)))))
    isnothing(i) || push!(output, vertices[i])
  end
  return output
end

function refine(A::Point, B::Point)
  M = (A + B) / 2
  return [(A, M), (Point(A.x, M.y), Point(M.x, B.y)),
          (Point(M.x, A.y), Point(B.x, M.y)), (M, B)]
end

function windingangle(vertices::AbstractVector{Vertex})
  θ = 0.0
  for (i, vertex) ∈ enumerate(vertices)
    j = periodicindex(vertices, i+1)
    θ += angle(vertex.value, vertices[j].value)
  end
  return θ
end

"""
The winding number of a complex function around a closed loop 
indicates the number of roots of the function encircled by the loop.

### Arguments
 * f: function that maps points in 2D space, ``R^2``,  to to a complex number
 * A: vector of real numbers that represents the lower left corner in ``R^2``
 * B: vector of real numbers that represents the upper right corner in ``R^2``

### Returns
The winding number of the function f: the number roots inside the square A->B
"""
function windingnumber(f::U, A::AbstractVector{T}, B::AbstractVector{T},
                       helper::Helper=Helper()) where {T<:Real, U<:Function}
  return windingnumber(distributevertices(f, Point(A), Point(B), helper))
end
function windingnumber(vertices::AbstractVector{Vertex})
  return windingangle(vertices)/(2π)
end


function boundingboxesofpoles!(solutions, f::T, A::Point, B::Point,
    vertices::AbstractVector{Vertex}, helper::Helper,
    xtol_rel, stopval, timelimit) where {T<:Function}
  t = @elapsed nfloat = round(windingnumber(vertices))
  n = isnan(nfloat) ? 0 : Int(nfloat)
  helper.timetaken += t
  ((helper.timetaken > timelimit) || iszero(n)) && return nothing
  isconverged = all(Vector(B - A) .< xtol_rel) ||
    (sum(abs.(f(Vector((A + B) / 2)))) <= stopval)
  if isconverged
    push!(solutions, (Vector(A), Vector(B)))
  else
    helper.numvertices = max(
      Int(round(helper.numvertices / 2)),
      Int(round(default_initialnumvertices / 2)))
    for q ∈ refine(A, B)
      newvertices = distributevertices(f, q..., helper, vertices)
      boundingboxesofpoles!(solutions, f, q..., newvertices, helper,
                            xtol_rel, stopval, timelimit)
    end
  end
  return nothing
end

function boundingboxesofpoles(f::T, A::AbstractVector{U}, B::AbstractVector{U};
                              kwargs...) where {T<:Function, U<:Real}
  @assert length(A) == length(B) == 2
  solutions = Vector{Tuple{Vector{Rational},Vector{Rational}}}()
  kwargs = Dict(kwargs)
  xtol_rel = get(kwargs, :xtol_rel, sqrt(eps()))
  stopval = get(kwargs, :stopval, sqrt(eps()))
  timelimit = get(kwargs, :timelimit, 600)
  initialnumvertices = get(kwargs, :initialnumvertices, default_initialnumvertices)
  minnumvertices = get(kwargs, :minnumvertices, min(initialnumvertices,
                                                    default_minnumvertices))
  helper = get(kwargs, :helper, Helper(initialnumvertices, minnumvertices))
  scale(fraction) = A .+ (B .- A) .* fraction
  f_local(fraction) = f(scale(fraction))
  lower = Point(zeros(Rational, 2))
  upper = Point(ones(Rational, 2))
  vertices = distributevertices(f_local, lower, upper, helper)
  boundingboxesofpoles!(solutions, f_local, lower, upper, vertices, helper,
                        xtol_rel, stopval, timelimit)
  return [(scale(solution[1]), scale(solution[2])) for solution ∈ solutions]
end
