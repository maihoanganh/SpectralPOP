"""
solving systems of polynomial equations
"""

function test_chemkin()

@polyvar x3 x4 y2 y3 y4 y5 z2 z3 z4 z5# variables
x=[x3; x4; y2; y3; y4; y5; z2; z3; z4; z5]


println("chemkin")

h=[9*y2^2 - 5.656854249492381*y2 + z2;
x3^2 + y3^2 + z3^2 - 1;
x4^2 + y4^2 + z4^2 - 1;
       y5^2 + z5^2 - 0.888888888888889;
x3 - 2.828427124746190*y2*x3 + y2*y3 + z2*z3 - 1/3;
x3*x4 + y3*y4 + z3*z4 - 1/3;
1/3*x4 + y4*y5 + z4*z5 - 1/3;
8/3 - 2.828427124746190*y2 + x3 + x4;
y2 + y3 + y4 + y5 + 0.8888888888888889;
z2 + z3 + z4 + z5]

L=1e4 # Squared radius of a ball containing at least one real root
k=2 # relaxed order

sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end
