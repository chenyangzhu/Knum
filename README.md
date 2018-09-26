# Knum

Numerical Analysis for Julia

To use this you can use the following line,

```julia
Pkg > add https://github.com/Klaus271/Knum/src/Knum.jl
```

Available Algorithm

- Point
  - Bisection
  - Fixed Point
  - Newton's Method
  - Secant
  - False Position
  - Steffensen's Method
  - Aitken's Delta Method
- Polynomials
  - Neville's Method

# How to use it

This package is still under development.

```julia
julia> f(x) = (3.)^x
f (generic function with 1 method)

julia> Knum.Poly.Neville(f,1/2,[-2,-1,0,1,2])
5×5 Array{Float64,2}:
 0.111111  0.0       0.0      0.0      0.0
 0.333333  0.666667  0.0      0.0      0.0
 1.0       1.33333   1.5      0.0      0.0
 3.0       2.0       1.83333  1.77778  0.0
 9.0       0.0       1.5      1.66667  1.70833
```

# Acknowledgement

This project is using ForwardDiff.jl.

You can read their paper [Forward-Mode Automatic Differentiation in Julia on arXiv](https://arxiv.org/abs/1607.07892), or [Github Link](https://github.com/JuliaDiff/ForwardDiff.jl)
