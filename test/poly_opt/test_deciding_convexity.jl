function getPOP_deciding_convexity(n,d,data)
    include(data*"/deciding_conv_var$(n)deg$(d).jl")
    
    println("***Problem setting***")

    println("Number of variables: n=",n)
    println("====================")

    @polyvar x[1:2*n]# variables

    # random quadratic objective function f
    v=reverse(monomials(x[1:n],2*d))
    obj=c'*v

    grad=differentiate.(obj,x[1:n])

    f=obj(x[1:n]=>x[n+1:2*n])-obj-grad'*(x[n+1:2*n]-x[1:n])


    # unit sphere constraint
    R=1.0
    h=[R-sum(x.^2)] #type of coefficients of each polynomial must be float

    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")
    return x,f,h,R
end

function test_test_deciding_convexity(n,d,data)
    
    x,f,h,R=getPOP_deciding_convexity(n,d,data)

    k=Int64(d)

    println("Relaxation order: k=",k)


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()



    g=Vector{Polynomial{true,Float64}}([])

    opt_val = SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    opt_val,opt_sol = CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Mix",tol=1e-5,scale=true) #Limited memory bundle method


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


function test_deciding_convexity(data)

    d=2
    N=[5;7;10;12]#;20;25;30]

    for n in N
        test_test_deciding_convexity(n,d,data)
    end
    test_test_deciding_convexity(5,3,data) 
    
end
