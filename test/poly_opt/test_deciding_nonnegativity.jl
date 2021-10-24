function getPOP_deciding_nonnegativity(n,d,data)
    include(data*"/deciding_nonneg_var$(n)deg$(d).jl")
    
    println("***Problem setting***")


    println("Number of variables: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    # random quadratic objective function f
    v=reverse(monomials(x,2*d))
    f=c'*v


    # unit sphere constraint
    R=1.0
    h=[R-sum(x.^2)] #type of coefficients of each polynomial must be float

    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")
    return x,f,h,R
end



function test_test_deciding_nonnegativity(n,d,data)
    
    x,f,h,R=getPOP_deciding_nonnegativity(n,d,data)

    
    

    k=Int64(d)

    println("Relaxation order: k=",k)


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()



    
    
    if n<=30
        g=Vector{Polynomial{true,Float64}}([])
        opt_val = SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek
    end

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


function test_deciding_nonnegativity(data)

    d=2
    N=[5;10;15;20;25;30;35;40;45;50]

    for n in N
        test_test_deciding_nonnegativity(n,d,data)
    end

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

    for j in 3:4
        test_test_deciding_nonnegativity(5,j,data)
    end
end
