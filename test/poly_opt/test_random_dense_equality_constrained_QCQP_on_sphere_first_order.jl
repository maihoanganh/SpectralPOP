function test_test_random_dense_equality_constrained_QCQP_on_sphere_first_order(n::Int64,data)


    println("***Problem setting***")

    l=ceil(Int32, n/4)+1

    println("Number of variable: n=",n)
    println("====================")

    println("Number of inequality constraints: l_g=",0)
    println("====================")

    println("Number of equality constraints: l_h=",l)
    println("====================")

    
    include(data*"/densePOPsphere_deg2_var$(n)_nineq0_neq$(l).jl")

    x,f,g,h=SpectralPOP.get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);
    
    
    k=Int64(1)

    println("Relaxed order: k=",k)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=300
        g=Vector{Polynomial{true,Float64}}([])

        opt_val = SpectralPOP.SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5) #Limited memory bundle method

    if n<=70
        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()

        opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-3)
    end
    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


function test_random_dense_equality_constrained_QCQP_on_sphere_first_order(data)


N=[50;60;70;80;100;120;150;200;300;400]

for n in N
    test_test_random_dense_equality_constrained_QCQP_on_sphere_first_order(n,data)
end
end