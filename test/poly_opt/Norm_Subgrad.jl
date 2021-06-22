function Norm_Subgrad(data)

println("***POP with single inequality (ball) constraint***")

n=10
m=1
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



opt_val,opt_sol = SpectralPOP.CTP_POP_UniqueBallIneq(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5,showNormGrad=true)


println()
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println()



println("***POP on a ball***")



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


opt_val,opt_sol = SpectralPOP.CTP_POP_on_Ball(x,f,g,h,k,R,showNormGrad=true)
end