function getPOP_Evaluation_comparisons(n::Int64,data)
    println("***Problem setting***")

    l=ceil(Int32, n/4)+1

    println("Number of variables: n=",n)
    println("====================")
        
    println("Number of inequality constraints: l_g=",0)
    println("====================")
        
    println("Number of equality constraints: l_h=",l)
    println("====================")
    
    
    include(data*"/densePOPsphere_deg2_var$(n)_nineq0_neq$(l).jl")

    x,f,g,h=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);

    return x,f,g,h,R
end



function test_Evaluation_comparisons(n::Int64,data)

    x,f,g,h,R=getPOP_Evaluation_comparisons(n,data)
    
    k=Int64(2)

    println("Relaxation order: k=",k)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()



    opt_val,opt_sol = CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Mix",tol=1e-5,showEvaluation=true,scale=true) #Limited memory bundle method


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

    opt_val,opt_sol = CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-3,showEvaluation=true,scale=true)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


function Evaluation_comparisons(data)


N=[5;10;15;20;25]

for n in N
    test_Evaluation_comparisons(n,data)
end
end
