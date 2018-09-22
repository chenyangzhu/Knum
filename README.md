# Knum

Numerical Analysis for Julia

To use this you can use the following line,

```julia
Pkg > add https://github.com/Klaus271/Knum/src/Knum.jl
```

Available Algorithm

- Bisection
- Fixed Point
- Newton's Method
- Secant
- False Position
- Steffensen's Method
- Aitken's Delta Method

# How to use it

This package is still under development.

```julia
julia> h(x::Vector) = 3^(-x[1])
h (generic function with 1 method)

julia> h(x) = 3^(-x)
h (generic function with 2 methods)

julia> Knum.fixed_point(h,0.,10)
0       0.0
1       1.0
2       0.3333333333333333
3       0.6933612743506348
4       0.46685562817398535
5       0.598760657911438
6       0.5179866461166681
7       0.5660536061282596
8       0.5369375702563746
9       0.5543903636000769
10      0.543861823031617
```

# Acknowledgement

This project is using ForwardDiff.jl.

You can read their paper [Forward-Mode Automatic Differentiation in Julia on arXiv](https://arxiv.org/abs/1607.07892), or [Github Link](https://github.com/JuliaDiff/ForwardDiff.jl)
