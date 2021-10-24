#=if isfile(joinpath(dirname(@__FILE__),"deps.jl"))
    include("./deps.jl")
else
    error("LMBM is not properly installed. Please run Pkg.build(\"LMBM\")")
end=#



"""
This function is for Limited memory bundle method with the following parameter:

        N             Number of variables.
        NA            Maximum bundle dimension, NA >= 2.
        MCU           Upper limit for maximum number of stored
                         corrections, MCU >= 3.

Real parameters in the Cdouble array of length 8
           RPAR(1)       Tolerance for change of function values.
           RPAR(2)       Second Tolerance for change of function values.
           RPAR(3)       Tolerance for the function value.
           RPAR(4)       Tolerance for the first termination criterion.
           RPAR(5)       Tolerance for the second termination criterion.
           RPAR(6)       Distance measure parameter, 0 <= RPAR(6).
           RPAR(7)       Line search parameter, 0 < RPAR(7) < 0.25.
           RPAR(8)       Maximum stepsize, 1 < RPAR(8).
                           If RPAR(I) <= 0 for I=1,3,4,5,7, and 8 the
                           default value of the parameter will be used.
                           If RPAR(2) < 0 the the parameter and the
                           corresponding termination criterion will be
                           ignored. If RPAR(2) = 0 default value will
                           be used. If RPAR(6) < 0 the default value
                           will be used.

Integer parameters in the int array of length 7
          IPAR(1)       Exponent for distance measure.
          IPAR(2)       Maximum number of iterations.
          IPAR(3)       Maximum number of function evaluations.
          IPAR(4)       Maximum number of iterations with changes of
                         function values smaller than RPAR(1).
          IPAR(5)       Printout specification:
                             -1  - No printout.
                              0  - Only the error messages.
                              1  - The final values of the objective
                                   function.
                              2  - The final values of the objective
                                   function and the most serious
                                   warning messages.
                              3  - The whole final solution.
                              4  - At each iteration values of the
                                   objective function.
                              5  - At each iteration the whole
                                   solution
          IPAR(6)       Selection of the method:
                              0  - Limited memory bundle method.
                              1  - L-BFGS bundle method.
          IPAR(7)       Selection of the scaling:
                              0  - Scaling at every iteration with STU/UTU.
                              1  - Scaling at every iteration with STS/STU.
                              2  - Interval scaling with STU/UTU.
                              3  - Interval scaling with STS/STU.
                              4  - Preliminary scaling with STU/UTU.
                              5  - Preliminary scaling with STS/STU.
                              6  - No scaling.
                           If IPAR(I) <= 0 the default value of the
                           parameter will be used.
"""

function lmbm(optFun::Function,x::Array{Float64};na::Int=2,mcu::Int=5,printinfo::Bool=false,maxtime::Float64=300000.0,tol::Float64=1e-7)


  global fOptPtr = @cfunction($optFun,Cdouble,(Cint,Ptr{Cdouble}, Ptr{Cdouble}))
  ccall((:set_obj_func,liblmbm),Cvoid,(Ptr{Cvoid},),fOptPtr)

  ipar=convert(Array{Cint},[0,5000000,5000000,0,-1,0,0])
 
  rpar=convert(Array{Cdouble},[0,0,0,tol,tol,0.5,0.24,2])

  maxtime = convert(Cfloat,maxtime);
  rtim=zeros(Cfloat,2)
  iout=zeros(Cint,3)
  n=convert(Cint,length(x));
  na =convert(Cint,na);
  mcu =convert(Cint,mcu);
  mc =convert(Cint,7);
  nw =convert(Cint,1 + 9*n + 2*n*na + 3*na + 2*n*(mcu+1) + 3*(mcu+2)*(mcu+1)/2 + 9*(mcu+1));
  w=zeros(Cdouble,nw);
  fVal=zeros(Cdouble,1);

  if printinfo
    @printf("---------------\n");
    @printf("| Parameters: |\n");
    @printf("---------------\n");
    @printf("%-8s %d\n", "n:", n);
    @printf("%-8s %f\n", "maxtime:", maxtime);
    @printf("%-8s %d\n", "na:", na);
    @printf("%-8s %d\n", "mcu:", mcu);
    @printf("%-8s %d\n", "mc:", mc);
    @printf("rpar: ");
    # printDoubleVector(rpar, 8);
    @printf("\n");
    @printf("ipar: ");
    # printIntVector(ipar, 7);
    @printf("\n");
  end

  ccall((:lmbmu_,liblmbm),Cvoid,(
      Ptr{Cint}, # n --- number of variables
      Ptr{Cint}, # na --- Maximum bundle dimension
      Ptr{Cint}, # mcu --- Upper limit for maximum number of stored
      Ptr{Cint}, # mc --- Upper limit for maximum number of stored
      Ptr{Cint}, # nw --- length of the vector n
      Ptr{Cdouble},# x --- initial solution
      Ptr{Cdouble}, # f --- return optimized value, which we cannot read
      Ptr{Cdouble}, # rpar --- real value parameters,
      Ptr{Cint}, # ipar --- integer value parameters
      Ptr{Cint}, # iout --- integer output parameters
      Ptr{Cfloat}, # maxtime --- maximum execution time
      Ptr{Cfloat}, # rtime --- returned timings
      Ptr{Cdouble}), # w --- working vector timings
      Ref(n), Ref(na), Ref(mcu), Ref(mc), Ref(nw), x,fVal, rpar, ipar,
         iout, Ref(maxtime), rtim, w)

  if printinfo
    @printf("-----------\n");
    @printf("| Output: |\n");
    @printf("-----------\n");
    @printf("%-16s %d\n", "Termination:", iout[3]);
    @printf("%-16s %d\n", "N. iter.:", iout[1]);
    @printf("%-16s %d\n", "N. func. eval.:", iout[2]);
    @printf("%-16s %f\n", "Final value:", fVal[1]);
    @printf("%-16s %f\n", "Execution time:", rtim[1]);
  end
   return fVal[1], x
end
