function getPOP_random_dense_QCQP_unique_inequality_ball_constraint(n::Int64,data)
    println("***Problem setting***")


    l=ceil(Int32, n/4)

    println("Number of variables: n=",n)
    println("====================")

    println("Number of inequality constraints: l_g=",1)
    println("====================")

    println("Number of equality constraints: l_h=",l)
    println("====================")
    
    include(data*"/densePOPsphere_deg2_var$(n)_nineq1_neq$(l).jl")

    x,f,g,h=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);
    return x,f,g,h,R
end



function test_test_random_dense_QCQP_unique_inequality_ball_constraint(n::Int64,data)
    
    x,f,g,h,R=getPOP_random_dense_QCQP_unique_inequality_ball_constraint(n,data)

    k=2

    println("Relaxation order: k=",k)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=30

        opt_val = SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = CTP_POP_UniqueBallIneq(x,f,h,k,R;method="LMBM",EigAlg="Mix",tol=1e-5,scale=true)

    if n<=25
        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
        try
            opt_val,opt_sol = BTP_POP(x,f,g,h,k,R,scale=true)
        catch
            println("Error with BTP_POP")
        end
    end
    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


function test_random_dense_QCQP_unique_inequality_ball_constraint(data)


N=[5;10;15;20;25;30;35;40;45;50]

for n in N
    test_test_random_dense_QCQP_unique_inequality_ball_constraint(n,data)
end
end
