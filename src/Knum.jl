module Knum

using ForwardDiff

module Point

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

end # end module Point


module Poly

function Neville(f,x,xi::Array)
    #=
    Input:
    f - function to value
    x - the initial point
    xi - Array of xs to evaluate

    Output:
    Q
    =#

    n = length(xi)-1
    Q = zeros(n+1,n+1)
    for i in 1:n+1
        Q[i,1] = f(xi[i])
    end
    for i in 1:n
        for j in 1:i
            Q[i+1,j+1] = ((x-xi[i-j+1])*Q[i+1,j] - (x-xi[i+1])*Q[i,j])/(xi[i+1]-xi[i-j+1])
        end
    end
    return Q
end

end # module Poly

module Diff

    function SimpleDiff(f,x,h)
        #=
        Input:
        f - function to evaluate
        x - point to evaluate
        h - step
        Output:
        Compute Foward/Backward Difference
        =#
        return (f(x+h)-f(x))/h
    end

    function ThreePointEnd(f,x,h)
        return (1/(2*h))*(-3*f(x)+4*f(x+h)-f(x+2h))
    end

    function ThreePointMid(f,x,h)
        return (1/(2*h))*(f(x+h)-f(x-h))
    end

    function FivePointMid(f,x,h)
        return (1/(12*h))*(f(x-2*h)-8*f(x-h)+8*(x+h)-f(x+2*h))
    end

    function FivePointEnd(f,x,h)
        return (1/(12*h))*(-25*f(x)+48*f(x-h)-36*f(x+2*h)+16f(x+3h)-3*f(x+4h))
    end

    function SecondMid(f,x,h)
        #=
        This function computes the second derivative using Midpoint method
        =#
        return 1/(h^2)*(f(x-h)-2*f(x)+f(x+h))
    end
end # Module Diff

module Integral
    function Legendre_coef(n)
        if n == 2
            return [0.5773502692, -0.5773502692],[1.,1.]
        elseif n == 3
            return [0.7745966692,0.,-0.7745966692],[0.5555555556, 0.8888888889, 0.5555555556]
        elseif n == 4
            return [0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116],[0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]
        elseif n == 5
            return [0.9061798459 0.5384693101 0.0000000000 -0.5384693101 -0.9061798459], [0.2369268850 0.4786286705 0.5688888889 0.4786286705 0.2369268850]
        end
    end

    function Trapezoidal(f,a,b)
        h = b-a
        return h/2*(f(a)+f(b))
    end

    function Simpsons(f,a,b)
        h=(b-a)/2
        return h/3*(f(a)+4*f(a+h)+f(b))
    end

    function CSimpsons(f,a,b,n)
        if n%2 == 0
            h = (b-a)/n
            odd = even = 0
            for j in 1:(n/2-1)
                odd += f(a+2*j*h)
                even += f(a+(2*j-1)*h)
            end
            even += f(a+(n-1)*h)
            return h/3*((f(a)+2*odd+4*even+f(b)))
        else
            Error("n should be even!")
        end
    end

    function CTrapezoidal(f,a,b,n)
        sum = 0
        h = (b-a)/n
        for j in 1:(n-1)
            x = a+j*h
            sum+=2*f(x)
        end
        return h/2*(sum+f(a)+f(b))
    end

    function CMidpoint(f,a,b,n)
        h = (b-a)/(n+2)
        sum = 0
        for j in 0:n/2
            sum += f(a+(2j+1)*h)
        end
        return sum*2*h
    end

    function Romberg(f,a,b,n)
        R = zeros(n,n)
        h = b-a
        R[1,1] = h/2*(f(a)+f(b))

        for row_number in 2:n
            sum = 0.
            for k in 1:2^(row_number-2)
                sum+=f(a+(k-0.5)*h)
            end
            R[row_number,1] = 0.5*(R[row_number-1,1]+h*sum)
            for j in 2:row_number
                R[row_number,j] = R[row_number,j-1]+(R[row_number,j-1]-R[row_number-1,j-1])/(4^(j-1)-1)
            end
            h = h/2
        end
        return R
    end

    function Adapt(f,a,b,tol,pts,n,return_detail=false)
        h = (b-a)/2
        x = [a:(b-a)/4:b;]
        fx = f.(x)
        n1 = n+5
        s = [1,0,4,0,1]' * fx * h/3
        s1 =[1,4,2,4,1]' * fx * h/6
        e1 = (s1-s)/15
        if abs(e1)<tol*2/3
            I = s1 + e1
            pts1 = [pts...,b]
        else
            I1,pts2,n2 = adapt(f,a,0.5*(a+b),0.5*tol,pts,n1,true)
            I2,pts1,n1 = adapt(f,0.5*(a+b),b,0.5*tol,pts2,n2,true)
            I = I1+I2
        end
        return_detail && return I,pts1,n1
        return I
    end

    function Gaussian(f,a,b,n)
        if n == 2
            tp = sqrt(3.)/3.
            tn = -sqrt(3.)/3.
            po = f((b+a)/2+tp*(b-a)/2)
            ne = f((b+a)/2+tn*(b-a)/2)
            return (po+ne)*(b-a)/2
        end
        if n == 5
            coe = [0.9061798459 0.5384693101 0.0000000000 -0.5384693101 -0.9061798459]
            c = [0.2369268850 0.4786286705 0.5688888889 0.4786286705 0.2369268850]
            trans_coe = ((b-a) .* coe .+ (b+a)) ./ 2
            ans = f.(trans_coe) * c' *(b-a)/2
            return ans[1]
        end
    end

    function SimpsonsDouble(f::Function,a::Float64,b::Float64,c::Function,d::Function,m::Int64,n::Int64)
        h = (b-a)/n
        J1 = 0
        J2 = 0
        J3 = 0
        for i in 0:n
            x = a+i*h
            HX = (d(x) - c(x))/m
            K1 = f(x,c(x)) + f(x,d(x))
            K2 = 0
            K3 = 0
            for j in 1:(m-1)
                y = c(x)+j*HX
                Q = f(x,y)
                if j%2 == 0
                    K2 = K2 + Q
                else
                    K3 = K3 + Q
                end
            end
            L = (K1 + 2*K2 + 4*K3)*HX/3
            if i == 0 || i == n
                J1 = J1 + L
            elseif i % 2 == 0
                J2 = J2 + L
            else
                J3 = J3 + L
            end # end if
        end # end for
        J = h*(J1 + 2*J2 + 4*J3)/3
        return J
    end  # end SimpsonsDouble

    function GaussianDouble(f::Function,a::Float64,b::Float64,c::Function,d::Function,m::Int64,n::Int64)
        rm,cm = Legendre_coef(m)
        rn,cn = Legendre_coef(n)
        h1 = (b-a)/2
        h2 = (b+a)/2
        J = 0
        for i in 1:m
            JX = 0
            x = h1 * rm[i] + h2
            d1 = d(x)
            c1 = c(x)
            k1 = (d1 - c1)/2
            k2 = (d1 + c1)/2
            for j in 1:n
                y = k1*rn[j]+k2
                Q = f(x,y)
                JX = JX + cn[j]*Q
            end
            J = J + cm[i]*k1*JX
        end # end for i
        J = h1*J
        return J
    end # end of GaussianDouble

end # module Integral

module ODE

    function Euler(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        @show (t,w)
        for i in 1:N
            w = w+h*f(t,w)
            t = t + h
            @show (t,w)
        end
        return w
    end # Euler

    function mEuler(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        for i in 1:N
            w = w+h/2*(f(t,w)+f(t+h,w+h*f(t,w)))
            t = t+h
        end
        return w
    end # mEuler

    function Taylor2(f,df,a::Float64,b::Float64,alpha::Float64,N)
        h = (b-a)/N
        t = a
        w = alpha
        @show (t,w)
        for j = 1:N
            w = w + h*f(t,w) + 0.5*h*h*df(t,w)
            t = t+h
            @show (t,w)
        end
        return w
    end # Taylor2

    function Taylor4(f,df,ddf,d3f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        h2 = h*h
        h3 = h*h2
        h4 = h2*h2
        @show t,w
        for j in 1:N
            w = w + h*f(t,w) + 0.5*h2*df(t,w)+h3/6*ddf(t,w)+h4/24*d3f(t,w)
            t = t+h
            @show t,w
        end
        return w
    end # End Taylor4

    function Midpoint(f,a,b,alpha,N)
        w = alpha
        h = (b-a)/N
        t = a
        @show (t,w)
        for i in 1:N
            w = w + h * f(t+h/2,w+h/2*f(t,w))
            t = t + h
            @show (t,w)
        end
        return w
    end # Taylor

    function mEuler(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        @show (t,w)
        for i in 1:N
            w = w+h/2*(f(t,w)+f(t+h,w+h*f(t,w)))
            t = t+h
            @show (t,w)
        end
        return w
    end # mEuler

    function Heun(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        @show (t,w)
        for i in 1:N
            w = w+h/4*(f(t,w)+3*f(t+2*h/3,w+2*h/3*f(t+h/3,w+h/3*f(t,w))))
            t = t+h
            @show (t,w)
        end
        return w
    end # Heun

    function RK4(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        @show (t,w)

        for i in 1:N
            k1 = h*f(t,w)
            k2 = h*f(t+h/2,w+1/2*k1)
            k3 = h*f(t+h/2,w+1/2*k2)
            k4 = h*f(t+h,w+k3)
            w = w+(k1+2*k2+2*k3+k4)/6
            t = t+h
            @show (t,w)
        end
        return w
    end # RK4

    function  AB3(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = zeros(N+1)
        w[1] = alpha
        for i in 2:3
            k1 = h*f(t,w[i-1])
            k2 = h*f(t+h/2,w[i-1]+1/2*k1)
            k3 = h*f(t+h/2,w[i-1]+1/2*k2)
            k4 = h*f(t+h,w[i-1]+k3)
            w[i] = w[i-1]+(k1+2*k2+2*k3+k4)/6
            t = t+h
        end
        for i in 4:N+1
            w[i] = w[i-1] + h/12 * (23*f(t,w[i-1])- 16*f(t-h,w[i-2])+ 5*f(t-2*h,w[i-3]))
            t = t + h
        end
        return w
    end # end AB3

    function AM3(f,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = zeros(N+1)
        w[1] = alpha
        for i in 2:3
            @show t,w
            k1 = h*f(t,w[i-1])
            k2 = h*f(t+h/2,w[i-1]+1/2*k1)
            k3 = h*f(t+h/2,w[i-1]+1/2*k2)
            k4 = h*f(t+h,w[i-1]+k3)
            w[i] = w[i-1]+(k1+2*k2+2*k3+k4)/6
            t = t + h
        end
        for i in 4:N+1
            @show t,w
            w[i] = w[i-1] + h/12 * (23*f(t,w[i-1]) - 16*f(t-h,w[i-2])+ 5*f(t-2*h,w[i-3]))
            w[i] = w[i-1] + h/24 * (9*f(t+h,w[i]) + 19*f(t,w[i-1]) - 5*f(t-h,w[i-2]) + f(t-2*h,w[i-3]))
            t = t + h
        end
        @show t,w
        return w
    end # AM3

    function rk4s(f,a,b,alpha,N)
        #=
        alpha is the initial values column vector
        f is the array of functions we use in column vector, m functions
        =#

        d = length(alpha)
        h = (b-a)/N
        t = a .+ [0:N;] .* h
        w = zeros(Float64,N+1,d)
        w[1,:] = alpha
        for j in 1:N
            k1 = f(t[j],w[j,:])
            k2 = f(t[j] .+ 0.5*h, w[j,:] .+ 0.5 * h .* k1)
            k3 = f(t[j] .+ 0.5*h, w[j,:] .+ 0.5 * h .* k2)
            k4 = f(t[j] .+     h, w[j,:] .+       h .* k3)
            w[j+1,:] = w[j,:] .+ h/6 * (k1 + 2* k2 + 2*k3 + k4)
            @show t[j], w[j+1,:]
        end
        return w,t
    end #End Rk4s

    function Trapezoidal(f,df,a,b,alpha,N,TOL,M = 1000)
        #=
        M is maximum number of iterations
        =#
        h = (b-a)/N
        t = a
        w = alpha

        for i in 1:N
            k1 = w+h/2*f(t,w)
            w0 = k1
            j = 1
            FLAG = 0
            while FLAG == 0
                w = w0- (w0 - h/2*f(t+h,w0)-k1)/(1-h/2*df(t+h,w0))
                if abs(w-w0) < TOL
                    FLAG = 1
                else
                    j = j+1
                    w0 = w
                    if j > M
                        println("The maximum number of iterations exceeded.")
                        return t,w
                    end # end if
                end # end if
            end # end while
            t = a+i*h
        end # end for
        return t,w
    end # end Trapezoidal

    function BkwdEulerStiff(lambda,a,b,alpha,N)
        h = (b-a)/N
        t = a
        w = alpha
        for i in 1:N
            w = w/(1-h*lambda)
        end
        return w
    end

end # End ODE

end # module Knum
