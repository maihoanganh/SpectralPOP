function test_random_dense_QCQP_unique_inequality_ball_constraint(data)

function test(n::Int64)

    println("***Problem setting***")


    l=ceil(Int32, n/4)

    println("Number of variable: n=",n)
    println("====================")

    println("Number of inequality constraints: l_g=",1)
    println("====================")

    println("Number of equality constraints: l_h=",l)
    println("====================")
    
    include(data*"/densePOPsphere_deg2_var$(n)_nineq1_neq$(l).jl")

    x,f,g,h=SpectralPOP.get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);

    k=2

    println("Relaxed order: k=",k)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=30

        opt_val = SpectralPOP.SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = SpectralPOP.CTP_POP_UniqueBallIneq(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5)

    if n<=20
        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()

        opt_val,opt_sol = SpectralPOP.BTP_POP(x,f,g,h,k,R)
    end
    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


N=[5;10;15;20;25;30;35;40;45;50]

for n in N
    test(n)
end
end