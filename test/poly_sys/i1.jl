"""
solving systems of polynomial equations
"""

function test_i1()


@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 # variables
x=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10]


println("i1")

h=[x1 - 0.25428722 - 0.18324757*x4*x3*x9;
 x2 - 0.37842197 - 0.16275449*x1*x10*x6;
 x3 - 0.27162577 - 0.16955071*x1*x2*x10;
 x4 - 0.19807914 - 0.15585316*x7*x1*x6;
 x5 - 0.44166728 - 0.19950920*x7*x6*x3;
 x6 - 0.14654113 - 0.18922793*x8*x5*x10;
 x7 - 0.42937161 - 0.21180484*x2*x5*x8;
 x8 - 0.07056438 - 0.17081208*x1*x7*x6;
 x9 - 0.34504906 - 0.19612740*x10*x6*x8;
 x10 - 0.42651102 - 0.21466544*x4*x8*x1]

L=1e4 # Squared radius of a ball containing at least one real root
k=3 # relaxed order


sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end
