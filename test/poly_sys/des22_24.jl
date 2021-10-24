"""
solving systems of polynomial equations
"""

function test_des22_24()


@polyvar a20 a21 a22 a23 a30 a31 a32 a33 a34 a35 u# variables
x=[a20;a21;a22;a23;a30;a31;a32;a33;a34;a35]


println("des22_24")

h=[16*a20*a32 + 18*a21*a31 + 20*a22*a30;
-80*a23 + 180*a34 + 855*a35;
 7*a20*a31 + 8*a21*a30;
 210*a35 - 210;
 40*a20*a34 + 44*a21*a33 + 48*a22*a32 + 52*a23*a31 + 280*a30;
 27*a20*a33 + 30*a21*a32 + 33*a22*a31 + 36*a23*a30;
 55*a20*a35 + 60*a21*a34 + 65*a22*a33 + 70*a23*a32 + 80*a30 + 375*a31;
 78*a21*a35 + 84*a22*a34 + 90*a23*a33 - 170*a20 + 102*a31 + 480*a32;
136*a23*a35 - 114*a22 + 152*a33 + 720*a34;
105*a22*a35 + 112*a23*a34 - 144*a21 + 126*a32 + 595*a33]

L=1e4 # Squared radius of a ball containing at least one real root
k=3 # relaxed order


sol=ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end
