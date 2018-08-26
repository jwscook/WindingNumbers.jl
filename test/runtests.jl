using Random
using WindingNumbers
using WindingNumbers: angle, Point
Random.seed!(10)

using Test

@testset "Winding Numbers tests" begin

@testset "Basic angle tests" begin
  a = Point(1.0, 0.0)
  b = Point(0.0, 1.0)
  @test π/2 ≈ angle(b) - angle(a)

  a = Point(-1.0, 0.0)
  b = Point(1.0, 0.0)
  @test -π ≈ angle(b) - angle(a)

  a = Point(1.0, 0.0)
  b = Point(-1.0, 0.0)
  @test π ≈ angle(b) - angle(a)

  a = Point(1.0, 0.0)
  b = Point(1.0, 1.0)
  @test π/4 ≈ angle(b) - angle(a)
 end

function mock(x::Vector{Float64}, n, r=zeros(2))
  θ = atan(x[2] - r[2], x[1] - r[1])
  return [cos(n*θ), sin(n*θ)]
end

turnings(n::Int) = -n:n

@testset "Basic understanding - angles" begin
  a = Point(1.0, 1.0)
  b = Point(-1.0, 1.0)
  c = Point(-1.0, -1.0)
  d = Point(1.0, -1.0)
  @test angle(a, b) ≈ π/2
  @test angle(b, c) ≈ π/2
  @test angle(c, d) ≈ π/2
  @test angle(d, a) ≈ π/2
  @test angle(b, a) ≈ -π/2
  @test angle(c, b) ≈ -π/2
  @test angle(d, c) ≈ -π/2
  @test angle(a, d) ≈ -π/2
  ϕ = angle(Point(-1.0, 2*eps()), Point(-1.0, -2*eps()))
  @test ϕ ≈ 4*eps()
  @test ϕ > 0.0
  ϕ = angle(Point(-1.0, -2*eps()), Point(-1.0, 2*eps()))
  @test ϕ ≈ -4*eps()
  @test ϕ < 0.0
end

@testset "Basic understanding - square with pole" begin
  a = Point(1.0, 1.0)
  b = Point(-1.0, 1.0)
  c = Point(-1.0, -1.0)
  d = Point(1.0, -1.0)
  ϕ = angle(a, b) + angle(b, c) + angle(c, d) + angle(d, a)
  @test ϕ ≈ 2π
  ϕ = angle(b, a) + angle(c, b) + angle(d, c) + angle(a, d)
  @test ϕ ≈ -2π
end

@testset "Basic understanding - square without pole" begin
  a = Point(1.0, 1.0)
  b = Point(2.0, 1.0)
  c = Point(2.0, 2.0)
  d = Point(1.0, 2.0)
  ϕ = angle(a, b) + angle(b, c) + angle(c, d) + angle(d, a)
  @test ϕ ≈ 0.0
  ϕ = angle(b, a) + angle(c, b) + angle(d, c) + angle(a, d)
  @test ϕ ≈ 0.0
end

@testset "Pole at origin in grid" begin
  for i in turnings(5)
    n = WindingNumbers.windingnumber(x->mock(x, i), -ones(2), ones(2))
    @test Int(round(n)) == i
  end
end

@testset "Random location of pole inside grid" begin
  for i in turnings(2)
    r = rand(2) # solution
    a, b = r - rand(2), r + rand(2)
    @assert all(a .< r .<b)
    n = WindingNumbers.windingnumber(x->mock(x, i, r), a, b)
    @test Int(round(n)) == i
  end
end

@testset "Unphysical pole" begin
  for i in turnings(5)
    s = rand(2) # square
    @assert s[1] > 0.0
    @assert s[2] > 0.0
    n = WindingNumbers.windingnumber(x->[0.0, 0.0], s, 2*s)
    @test Int(round(n)) == 0
  end
end

@testset "Physical pole not in grid" begin
  for i in turnings(5)
    s = rand(2) # square
    @assert s[1] > 0.0
    @assert s[2] > 0.0
    n = WindingNumbers.windingnumber(x->mock(x, i), s, 2*s)
    @test Int(round(n)) == 0
  end
end

@testset "solve" begin
  for i ∈ 1:100
    rtol = 10^Float64(-rand(1:14))
    root = 2 .* rand(2) .- 1
    lower = root .- rand(2)
    upper = root .+ rand(2)
    allpointsinside = true
    function solvemock(x)
      allpointsinside &= all(lower .<= x .<= upper) 
      y = mock(x, 1, root)
      return y[1] + im * y[2]
    end
    solutions = WindingNumbers.boundingboxesofpoles(solvemock, lower, upper,
                                                    xtol_rel=rtol, stopval=eps())
    @test length(solutions) == 1
    @test allpointsinside
    length(solutions) == 1 || continue
    A, B = solutions[1]
    @test all(A .<= root .<= B)
    @test all((B .- A) ./ (upper .- lower) .< rtol)
  end
end

end
