function test_random_dense_QCQP_on_ball(data)


function test(n::Int64)

    println("***Problem setting***")

    m=ceil(Int32, n/8)
    l=ceil(Int32, n/8)

    println("Number of variable: n=",n)
    println("====================")

    println("Number of inequality constraints: l_g=",m)
    println("====================")

    println("Number of equality constraints: l_h=",l)
    println("====================")

    include(data*"/densePOPsphere_deg2_var$(n)_nineq$(m)_neq$(l).jl")

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

    opt_val,opt_sol = SpectralPOP.CTP_POP_on_Ball(x,f,g,h,k,R)

    if n<=20
        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()

        opt_val,opt_sol = SpectralPOP.BTP_POP(x,f,g,h,k,R,EigAlg="Mix",tol=1e-5,showNormGrad=false)
    end
    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


N=[10;15;20;25;30;35]

for n in N
    test(n)
end
end