function getPOP_random_dense_quadratic_on_sphere(n::Int64,data)
    
    println("***Problem setting***")


    println("Number of variables: n=",n)
    println("====================")

    println("Number of inequality constraints: l_g=",0)
    println("====================")

    println("Number of equality constraints: l_h=",1)
    println("====================")

    
    include(data*"/densePOPsphere_deg2_var$(n)_nineq0_neq1.jl")

    x,f,g,h=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);
    return x,f,g,h,R
end    


function test_test_random_dense_quadratic_on_sphere(n::Int64,data)


    x,f,g,h,R=getPOP_random_dense_quadratic_on_sphere(n,data)
    
    k=1 # relaxed order

    println("Relaxation order: k=",k)


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=300
        g=Vector{Polynomial{true,Float64}}([])

        opt_val = SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Mix",tol=1e-5,scale=true) #Limited memory bundle method

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

    opt_val,opt_sol = CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-4,scale=true)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


function test_random_dense_quadratic_on_sphere(data)

N=[50;75;100;125;150;175;200;250;300;350;400;500;700;900;1200;1500]

for n in N
    test_test_random_dense_quadratic_on_sphere(n,data)
end
end
