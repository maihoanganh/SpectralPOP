function getPOP_deciding_copositivity(n,data)
    
    include(data*"/deciding_copos_size$(n).jl")
    
    println("***Problem setting***")

    println("Number of variables: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    # random quadratic objective function f
    f=sum(A[i,j]*x[i]^2*x[j]^2 for i=1:n for j=1:n)


    # unit sphere constraint
    R=1.0
    h=[R-sum(x.^2)] #type of coefficients of each polynomial must be float

    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")
    return x,f,h,R
end


function test_test_deciding_copositivity(n,data)

    x,f,h,R=getPOP_deciding_copositivity(n,data)
    

    k=Int64(2)

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


function test_deciding_copositivity(data)

N=[10;15;20;25]#;30]

for n in N
    test_test_deciding_copositivity(n,data)
end
end
