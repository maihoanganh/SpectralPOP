module SpectralPOP


using SumOfSquares, DynamicPolynomials, LinearAlgebra, MosekTools, JuMP, Arpack, RowEchelon, Libdl, Printf, Compat, OSQP, SparseArrays


export CTP_POP, ASC_PolySys


#solvers
include("../solvers/LMBM/build.jl")
include("../solvers/LMBM/deps.jl")
include("../solvers/LMBM/LMBM.jl")

include("../solvers/ProximalBundleMethod/ProximalMethod.jl")

include("../solvers/SketchyCGAL/NystromSketch.jl")
include("../solvers/SketchyCGAL/CGAL.jl")

# src
include("basicfuncs.jl")
include("POP_CTP.jl")
include("POP_BTP.jl")
include("extrac_optimizers.jl")
include("POP_SumOfSquares.jl")
include("SpectralASC_PolynomialSystems.jl")




#Run tests
include("../test/test_examples.jl")
#Polynomial optimization
include("../test/poly_opt/test_random_dense_quadratic_on_sphere.jl")
include("../test/poly_opt/test_random_dense_equality_constrained_QCQP_on_sphere_first_order.jl")
include("../test/poly_opt/test_random_dense_equality_constrained_QCQP_on_sphere_second_order.jl")
include("../test/poly_opt/test_random_dense_QCQP_unique_inequality_(ball)_constraint.jl")
include("../test/poly_opt/test_random_dense_QCQP_on_ball.jl")
include("../test/poly_opt/test_random_dense_quartics_on_sphere.jl")
include("../test/poly_opt/Evaluation_comparisons.jl")
include("../test/poly_opt/Norm_Subgrad.jl")
include("../test/poly_opt/test_deciding_convexity.jl")
include("../test/poly_opt/test_deciding_nonnegativity.jl")
include("../test/poly_opt/test_deciding_copositivity.jl")

#Polynomial systems
include("../test/poly_sys/chemkin.jl")
include("../test/poly_sys/d1.jl")
include("../test/poly_sys/des22_24.jl")
include("../test/poly_sys/i1.jl")
include("../test/poly_sys/katsura10.jl")
include("../test/poly_sys/katsura9.jl")
include("../test/poly_sys/katsura8.jl")
include("../test/poly_sys/katsura7.jl")
include("../test/poly_sys/katsura6.jl")
include("../test/poly_sys/kin1.jl")
include("../test/poly_sys/ku10.jl")
include("../test/poly_sys/pole27sys.jl")
include("../test/poly_sys/pole28sys.jl")
include("../test/poly_sys/stewgou.jl")


end
