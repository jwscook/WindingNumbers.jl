struct Point{T<:Real}
  x::T
  y::T
end
Point(x::T) where {T<:Complex} = Point(real(x), imag(x))
Point(p::Point{T}) where {T<:Real} = Point(p.x, p.y)
function Point(x::AbstractVector{T}) where {T}
  @assert length(x) == 2
  return Point(x[1], x[2])
end

Base.:+(a::Point, b::Point) = Point(a.x + b.x, a.y + b.y)
Base.:-(a::Point, b::Point) = Point(a.x - b.x, a.y - b.y)
Base.:*(a::Point, f::Vector{T}) where {T<:Number} = Point(a.x * f[1], a.y * f[2])
Base.:/(a::Point{}, f::T) where {T<:Number} = Point(a.x / f, a.y / f)

Vector(p::Point) = [p.x, p.y]

middle(a::Point, b::Point) = Point(0.5*(a.x + b.x), 0.5*(a.y + b.y))

angle(A::Point) = atan(A.y, A.x)

function angle(A::Point, B::Point)
  ϕ = angle(B) - angle(A)
  ϕ < -π && return ϕ + 2π
  ϕ > π && return ϕ - 2π
  return ϕ
end

