# WindingNumbers.jl
![CI](https://github.com/jwscook/WindingNumbers.jl/workflows/CI/badge.svg)
[![Build Status](https://travis-ci.com/jwscook/WindingNumbers.jl.svg?branch=master)](https://travis-ci.com/jwscook/WindingNumbers.jl)
[![codecov.io](http://codecov.io/github/jwscook/WindingNumbers.jl/coverage.svg?branch=master)](http://codecov.io/github/jwscook/WindingNumbers.jl?branch=master)


The winding number of a the mapping of ``R^2`` to ``C^1`` around a closed rectangular loop
indicates the number of roots of the function encircled by the loop. The contour is rectangular for a bifurcating search to find the location of the roots.

`windingnumber` return the number of roots of a function ``f`` around the rectangular path with lower left corner `A` and uppert corner `B`.
