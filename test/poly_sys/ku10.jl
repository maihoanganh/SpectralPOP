"""
solving systems of polynomial equations
"""

function test_ku10()


@polyvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 # variables
x=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10]


println("ku10")

h=[5*x1*x2+ 5*x1+ 3*x2+ 55;
 7*x2*x3+ 9*x2+ 9*x3+ 19;
 3*x3*x4+ 6*x3+ 5*x4-4;
 6*x4*x5+ 6*x4+ 7*x5+ 118;
x5*x6+ 3*x5+ 9*x6+ 27;
 6*x6*x7+ 7*x6+x7+ 72;
 9*x7*x8+ 7*x7+x8+ 35;
 4*x8*x9+ 4*x8+ 6*x9+ 16;
 8*x9*x10+ 4*x9+ 3*x10-51;
 3*x1*x10-6*x1+x10+ 5]

L=1e4 # Squared radius of a ball containing at least one real root
k=1 # relaxed order



sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end
