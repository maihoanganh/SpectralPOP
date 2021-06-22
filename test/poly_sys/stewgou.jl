function test_stewgou()



@polyvar x1 x2 x3 y1 y2 y3 z1 z2 z3# variables
x=[x1;x2;x3;y1;y2;y3;z1;z2;z3]

println("stewgou")
h=[x1^2+y1^2+z1^2-31.0;
x2^2+y2^2+z2^2-39;
x3^2+y3^2+z3^2-29;
x1*x2+y1*y2+z1*z2+6*x1-6*x2-51;
x1*x3+y1*y3+z1*z3+7*x1-2*y1-7*x3+2*y3-50;
x2*x3+y2*y3+z2*z3+x2-2*y2-x3+2*y3-34;
-12*x1+15*y1-10*x2-25*y2+18*x3+18*y3+32;
-14*x1+35*y1-36*x2-45*y2+30*x3+18*y3-8;
2*x1+2*y1-14*x2-2*y2+8*x3-y3-20]


L=1e4 # Squared radius of a ball containing at least one real root
k=2 # relaxed order

sol=SpectralPOP.ASC_PolySys(x,h,k,L,method="LMBM",EigAlg="Arpack",tol=1e-3);
end