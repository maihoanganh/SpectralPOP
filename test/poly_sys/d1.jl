"""
solving systems of polynomial equations
"""

function test_d1()


@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 # variables
x=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12]

println("d1")

h=[ x1^2  + x2^2 - 1;
 x3^2  + x4^2 - 1;
 x5^2  + x6^2 - 1;
 x7^2  + x8^2 - 1;
 x9^2  + x10^2 - 1;
 x11^2 + x12^2 - 1;
 3*x3 + 2*x5 + x7 - 3.9701;
 3*x1*x4 + 2*x1*x6 + x1*x8 - 1.7172;
 3*x2*x4 + 2*x2*x6 + x2*x8 - 4.0616;
 x3*x9 + x5*x9 + x7*x9 - 1.9791;
 x2*x4*x9 + x2*x6*x9 + x2*x8*x9 + x1*x10 - 1.9115;
 - x3*x10*x11 - x5*x10*x11 - x7*x10*x11 + x4*x12 + x6*x12 + x8*x12 - 0.4077]

L=1e4 # Squared radius of a ball containing at least one real root
k=2 # relaxed order

sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end
