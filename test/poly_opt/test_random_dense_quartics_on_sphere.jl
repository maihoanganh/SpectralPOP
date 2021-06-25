function test_test_random_dense_quartics_on_sphere(n::Int64,data)


    println("***Problem setting***")


    println("Number of variable: n=",n)
    println("====================")

    println("Number of inequality constraints: l_g=",0)
    println("====================")

    println("Number of equality constraints: l_h=",1)
    println("====================")


    
    include(data*"/densePOPsphere_deg4_var$(n)_nineq0_neq1.jl")

    x,f,g,h=SpectralPOP.get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);

    
        
    k=Int64(2)

    println("Relaxed order: k=",k)


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()



    g=Vector{Polynomial{true,Float64}}([])

    opt_val = SpectralPOP.SumofSquares_POP2(x,f,g,h,k) # SumOfSquares.jl + Mosek

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5) #Limited memory bundle method

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-3)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


function test_random_dense_quartics_on_sphere(data)


N=[5;10;15;20;25;30]

for n in N
    test_test_random_dense_quartics_on_sphere(n,data)
end
end