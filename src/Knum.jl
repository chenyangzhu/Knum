module Knum

using ForwardDiff

function fixed_point(equa, p, N = 20, show = false)
    '''
    Input:
    equa: The function you are going to compute
    p:    First point
    N:    Max # of iteration
    show: Whether or not to show calculation process

    Output:
    Approximation
    '''
    for i in 1:N+1
        if show
            println(i-1,'\t',p)
        end
        p = equa(p)
    end
    return p
end

function newton(equa, p, N = 20, tol = 1e-7, show = false)
    '''
    Input:
    equa: The function you are going to compute
    p:    First point
    N:    Max # of iteration
    tol:  Tolerance
    show: Whether or not to show calculation process

    Output:
    Approximation
    '''
    for i in 1:N+1
        println(i-1,'\t',p)
        dfp = ForwardDiff.gradient(equa,[p])[1]
        if abs(equa(p)/dfp) < tol
            if show:
                println(i,'\t', p - equa(p)/dfp)
            end
            return p - equa(p)/dfp
        else
            p = p - equa(p)/dfp
        end
    end
    return p
end

function secant(equa,p0,p1,N=20, tol = 1e-7 , show = false)
    '''
    Input:
    equa:  The function you are going to compute
    p0,p1: First  and second point
    N:     Max # of iteration
    tol:   Tolerance
    show:  Whether or not to show calculation process

    Output:
    Approximation

    '''
    if show:
        println(0,'\t',p0)
    end
    for i in 1:N-1
        if show:
            println(i,'\t',p1)
        end
        if abs(p1-p0) < tol
            return p1
        end
        dumy = p1
        p1 = p1 - (equa(p1)*(p1-p0))/(equa(p1)-equa(p0))
        p0 = dumy
    end
    if show:
        println(N,'\t',p1)
    end
    return p1
end

function fpst(equa,p0,p1,N, tol = 1e-7, show = false)
    '''
    Input:
    equa:  The function you are going to compute
    p0,p1: First  and second point
    N:     Max # of iteration
    tol:   Tolerance
    show:  Whether or not to show calculation process

    Output:
    Approximation
    '''
    for i in 1:N+1
        if show:
            println(i-1,'\t',p1)
        end
        # compute derivative
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

end # module
