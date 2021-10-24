function test_examples()

println("****Test polynomial optimization****")

@polyvar x[1:2] # variables

f=x[1]^2+0.5*x[1]*x[2]-0.25*x[2]^2+0.75*x[1]-0.3*x[2] # objective function

R=1.0 # squared radius of a sphere constraint
h=[R-sum(x.^2);(x[1]-1.0)*x[2]] # equality constraints (including the sphere constraint)

k=2 # relaxed order


g=Vector{Polynomial{true,Float64}}([])
opt_val = SumofSquares_POP_WithExtraction(x,f,g,h,k) # SumOfSquares.jl + Mosek
println()
println(".................................")
println()


# get approximations of the optimal value and an optimal solution
opt_val,opt_sol = CTP_POP(x,f,h,k,R,method="LMBM",EigAlg="Mix",tol=1e-5)# Limited Memory Bundle Method
println()
println(".................................")
println()
opt_val,opt_sol = CTP_POP(x,f,h,k,R,method="SketchyCGAL",EigAlg="Normal",tol=1e-3)# Limited Memory Bundle Method

# println()
# println()
# println()
# println()
# println()

# println("****Test polynomial system****")


# @polyvar x[1:2] # variables

# # mickey equations
# h=[x[1]^2 + 4*x[2]^2 - 4;
#         2*x[2]^2 - x[1]]

# L=10 # squared radius of a ball centered at origin containing at least one real root
# k=1 # relaxed order

# # get a real root
# root = ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-5)
end
