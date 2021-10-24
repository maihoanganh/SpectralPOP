"""
solving systems of polynomial equations
"""


function test_katsura9()

@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 # variables
x=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10]



println("katsura9")

h=[-x1+2*x10^2+2*x9^2+2*x8^2+2*x7^2+2*x6^2+2*x5^2+2*x4^2+2*x3^2+2*x2^2+x1^2;
 -x2+2*x10*x9+2*x9*x8+2*x8*x7+2*x7*x6+2*x6*x5+2*x5*x4+2*x4*x3+2*x3*x2+2*x2*x1;
 -x3+2*x10*x8+2*x9*x7+2*x8*x6+2*x7*x5+2*x6*x4+2*x5*x3+2*x4*x2+2*x3*x1+x2^2;
 -x4+2*x10*x7+2*x9*x6+2*x8*x5+2*x7*x4+2*x6*x3+2*x5*x2+2*x4*x1+2*x3*x2;
 -x5+2*x10*x6+2*x9*x5+2*x8*x4+2*x7*x3+2*x6*x2+2*x5*x1+2*x4*x2+x3^2;
 -x6+2*x10*x5+2*x9*x4+2*x8*x3+2*x7*x2+2*x6*x1+2*x5*x2+2*x4*x3;
 -x7+2*x10*x4+2*x9*x3+2*x8*x2+2*x7*x1+2*x6*x2+2*x5*x3+x4^2;
 -x8+2*x10*x3+2*x9*x2+2*x8*x1+2*x7*x2+2*x6*x3+2*x5*x4;
 -x9+2*x10*x2+2*x9*x1+2*x8*x2+2*x7*x3+2*x6*x4+x5^2;
 -1+2*x10+2*x9+2*x8+2*x7+2*x6+2*x5+2*x4+2*x3+2*x2+x1]

L=1e4 # Squared radius of a ball containing at least one real root
k=2 # relaxed order

sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end
