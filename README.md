# Knum

Numerical Analysis for Julia

To use this you can use the following line,

```julia
Pkg > add https://github.com/chenyangzhu/Knum
```

Available Algorithm

- Point Approximation
  - Bisection
  - Fixed Point
  - Newton's Method
  - Secant's Method
  - False Position
  - Steffensen's Method
  - Aitken's Delta Method
- Polynomials
  - Neville's Method
- Differentiation
  - three and five Mid/End point Method
  - forward/backward differentiation
- Integration
  - Simpson's Method
  - Trapezoidal Method
- ODE's
  - Euler's Method
  - AB/AM Explicit Methods
- Matrix
  - LU Factorization

# How to use it

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
