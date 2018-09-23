module Knum

using ForwardDiff

function bisection(f,a,b,n=100,tol=1e-7)
    p = 0.
    FA = f(a)
    for i in 1:n
        p = a+(b-a)/2
        FP = f(p)
        println(i,'\t',p)
        if FP == 0 || (b-a)/2 < tol
            return p
        end
        if sign(FA) * sign(FP) >0
            a = p
            FA = FP
        else
            b = p
        end
    end
    println("Method failed after iterations = ",n)
    return p
end

function fixed_point(g, p, N = 20, tol = 1e-7, show = false)
    #=
    Input:
    equa: The function you are going to compute
    p:    First point
    N:    Max # of iteration
    show: Whether or not to show calculation process

    Output:
    Approximation
    =#

    for i in 1:N+1
        if show
            println(i-1,'\t',p)
        end
        if abs(p-g(p))<tol
            return g(p)
        end
        p = g(p)
    end
    println("Method Failed After", n, "iterations, the result is returned.")
    return p
end

function newton(f, p, N = 20, tol = 1e-7, show = false)
    #=
    ---NOTICE---
    For Newton\'s method, your function need to enable 2 inputs,
    one with Float64, the other one with Vector.


    Input:
    equa: The function you are going to compute
    p:    First point
    N:    Max # of iteration
    tol:  Tolerance
    show: Whether or not to show calculation process

    Output:
    Approximation
    =#
    for i in 1:N+1
        println(i-1,'\t',p)
        dfp = ForwardDiff.gradient(equa,[p])[1]
        if abs(equa(p)/dfp) < tol
            if show
                println(i,'\t', p - equa(p)/dfp)
            end
            return p - f(p)/dfp
        else
            p = p - f(p)/dfp
        end
    end
    println("Method Failed After", n, "iterations, the result is returned.")
    return p
end

function secant(equa,p0,p1,N=20, tol = 1e-7 , show = false)
    #=
    Input:
    equa:  The function you are going to compute
    p0,p1: First  and second point
    N:     Max # of iteration
    tol:   Tolerance
    show:  Whether or not to show calculation process

    Output:
    Approximation

    =#
    if show
        println(0,'\t',p0)
    end
    for i in 1:N-1
        if show
            println(i,'\t',p1)
        end
        if abs(p1-p0) < tol
            return p1
        end
        dumy = p1
        p1 = p1 - (equa(p1)*(p1-p0))/(equa(p1)-equa(p0))
        p0 = dumy
    end
    if show
        println(N,'\t',p1)
    end
    return p1
end

function fpst(equa,p0,p1,N, tol = 1e-7, show = false)
    #=
    Input:
    equa:  The function you are going to compute
    p0,p1: First  and second point
    N:     Max # of iteration
    tol:   Tolerance
    show:  Whether or not to show calculation process

    Output:
    Approximation
    =#
    for i in 1:N+1
        if show
            println(i-1,'\t',p1)
        end
        dfp = equa(p1) * (p1 - p0) / (equa(p1) - equa(p0))
        if abs(dfp) < tol
            if show
                println(i,'\t', p1 - dfp)
            end
            return p1 - dfp
        else
            dumy = p1
            p1 = p1 - dfp
            p0 = dumy
        end
    end
    return p1
end

function Aitkens_delta(g, p, n)
    #=
    Input
    g - function
    p - First Point
    n - # of the sequence wanted

    Output
    aitken sequence
    =#

    # First compute the original sequence
    original_value = zeros(n+2)
    original_value[1] = g(p)
    for i in 1:n+1
        original_value[i+1] = g(original_value[i])
    end

    aitken_value = zeros(n)
    for i in 1:n
        aitken_value[i] = original_value[i] - ((original_value[i+1]-original_value[i])^2)/(original_value[i+2]-2*original_value[i+1]+original_value[i])
    end

    return aitken_value
end

function Steffensens(g, p0, n=10, tol=1e-7,show = false)
    #=
    Input
    g - function
    p - First Point
    n - # of the iteration
    tol - Tolerance

    Output
    Approximated Number
    =#
    for i in 1:n
        p1 = g(p0)
        p2 = g(p1)
        p = p0 - (p1-p0)^2/(p2-2p1+p0)
        if show
            println(i,'\t',p)
        end
        if abs(p - p0)<tol
            return p
        end
        p0 = p
    end
    println("Method Failed After", n, "iterations, the result is returned.")
    return p
end

function Muller(f, p0,p1,p2,n,tol, show = false)
    #=
    Input:
    p0,p1,p2 - Initial three guesses
    n - max # of iterations
    tol - tolerance

    Output:
    Roots
    =#

    # First make every input Complex
    p0 = Complex(p0)
    p1 = Complex(p1)
    p2 = Complex(p2)

    h1 = p1 - p0
    h2 = p2 - p1
    delta1 = (f(p1) - f(p0))/h1
    delta2 = (f(p2) - f(p1))/h2
    d = (delta2 - delta1)/(h1 + h2)

    # Create some temps
    b = Complex(0.)
    D = Complex(0.)
    E = Complex(0.)
    p = Complex(0.)

    # Start iteration
    for i in 3:n
        b = delta2 + h2 * d
        D = sqrt(b^2 - 4*f(p2)*d)
        if abs(b - D) < abs(b+D)
            E = b + D
        else
            E = b - D
        end
        h = -2* f(p2) / E
        p = p2 + h
        if abs(h)<tol
            return p
        end
        if show
            println(i,'\t',p)
        end

        # Prepare for next iteration
        p0 = p1
        p1 = p2
        p2 = p
        h1 = p1 - p0
        h2 = p2 - p1
        delta1 = (f(p1) - f(p0))/h1
        delta2 = (f(p2) - f(p1))/h2
        d = (delta2 - delta1)/(h1 + h2)
    end
    print("Method failed after N0 iterations, N0 =", n)
    return p
end

end # module
