function Norm_Subgrad(data)

println("***POP with single inequality (ball) constraint***")

n=10
m=1
l=ceil(Int32, n/4)
    



println("Number of variables: n=",n)
println("====================")

println("Number of inequality constraints: l_g=",1)
println("====================")

println("Number of equality constraints: l_h=",l)
println("====================")


include(data*"/densePOPsphere_deg2_var$(n)_nineq1_neq$(l).jl")

x,f,g,h=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);

k=2

println("Relaxation order: k=",k)



opt_val,opt_sol = CTP_POP_UniqueBallIneq(x,f,h,k,R;method="LMBM",EigAlg="Mix",showNormGrad=true,scale=true)


println()
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println()



println("***POP on a ball***")



m=ceil(Int32, n/8)
l=ceil(Int32, n/8)


println("Number of variables: n=",n)
println("====================")

println("Number of inequality constraints: l_g=",m)
println("====================")

println("Number of equality constraints: l_h=",l)
println("====================")


include(data*"/densePOPsphere_deg2_var$(n)_nineq$(m)_neq$(l).jl")

x,f,g,h=get_POP(n,m,l,lmon_g,supp_g,coe_g,lmon_h,supp_h,coe_h,lmon_f,supp_f,coe_f);

k=2

println("Relaxation order: k=",k)

println()
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println()


opt_val,opt_sol = CTP_POP_on_Ball(x,f,g,h,k,R,EigAlg="Mix",showNormGrad=true,scale=true)
end
