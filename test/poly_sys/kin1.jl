"""
solving systems of polynomial equations
"""


function test_kin1()

@polyvar s1 s2 s3 s4 s5 s6 c1 c2 c3 c4 c5 c6 # variables
x=[s1;s2;s3;s4;s5;s6;c1;c2;c3;c4;c5;c6]


println("kin1")

h=[s1^2 + c1^2 - 1;
 s2^2 + c2^2 - 1;
 s3^2 + c3^2 - 1;
 s4^2 + c4^2 - 1;
 s5^2 + c5^2 - 1;
 s6^2 + c6^2 - 1;
 s2*c5*s6 - s3*c5*s6 - s4*c5*s6 + c2*c6 + c3*c6 + c4*c6 - 0.4077;
 c1*c2*s5 + c1*c3*s5 + c1*c4*s5 + s1*c5 - 1.9115;
 s2*s5 + s3*s5 + s4*s5 - 1.9791;
 c1*c2 + c1*c3 + c1*c4 + c1*c2 + c1*c3 + c1*c2 - 4.0616;
 s1*c2 + s1*c3 + s1*c4 + s1*c2 + s1*c3 + s1*c2 - 1.7172;
 s2 + s3 + s4 + s2 + s3 + s2 - 3.9701]

L=1e4 # Squared radius of a ball containing at least one real root
k=2 # relaxed order

sol=SpectralPOP.ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end