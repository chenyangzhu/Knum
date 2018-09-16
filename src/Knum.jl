module Knum

using ForwardDiff

function fixed_point(equa, p, N)
    for i in 1:N+1
        println(i-1,'\t',p)
        p = equa(p)
    end
end

function newton(equa,p,N,tol = 1e-7)
    for i in 1:N+1
        println(i-1,'\t',p)
        dfp = ForwardDiff.gradient(equa,[p])[1]
        if abs(equa(p)/dfp) < tol
            println(i,'\t', p - equa(p)/dfp)
            return p - equa(p)/dfp
        else
            p = p - equa(p)/dfp
        end
    end
    return p
end

function secant(equa,p0,p1,N, tol = 1e-5)
    println(0,'\t',p0)
    for i in 1:N-1
        println(i,'\t',p1)
        if abs(p1-p0) < tol
            return p1
        end
        dumy = p1
        p1 = p1 - (equa(p1)*(p1-p0))/(equa(p1)-equa(p0))
        p0 = dumy
    end
    println(N,'\t',p1)
    return p1
end

function fpst(equa,p0,p1,N, tol = 1e-7)
    for i in 1:N+1
        println(i-1,'\t',p1)
        # compute derivative
        dfp = equa(p1) * (p1 - p0) / (equa(p1) - equa(p0))
        if abs(dfp) < tol
            println(i,'\t', p1 - dfp)
            return p1 - dfp
        else
            dumy = p1
            p1 = p1 - dfp
            p0 = dumy
        end
    end
end

end # module
