# SpectralSOS
- SpectralSOS is a Julia package of solving equality constrained polynomial optimization problems (POPs) on an Euclidean sphere:

**inf_x { f(x) : x in R^n, hj(x) = 0, j = 1,...,l } with h1 := R - |x|^2.**

- The main idea of SpectralSOS is to solve an SDP (moment) relaxation of the form:

**v = sup_X { <C,X> : X is psd, AX = b },**

which has constant trace property:

**AX = b => trace(X) = a,**

by using spectral (the largest eigenvalue) minimization:

**v = inf_z { aλ1(C - A'z) + b'z }**

with Limited memory bundle method instead of costly Interior-point methods.

- Compared to SumOfSquares (Mosek) and SketchyCGAL, SpectralSOS is cheaper, faster, but maintains the same accuracy as SumOfSquares on a tested sample of random dense equality constrained QCQPs on the unit sphere.

# Required softwares
SpectralSOS has been implemented on a desktop computer with the following softwares:
- Ubuntu 18.04.4
- Julia 1.3.1
- Fortran 2018

The following sofwares are used for comparison purposes:
- [SumOfSquares.jl](https://github.com/JuliaOpt/SumOfSquares.jl)
- [Mosek 9.1](https://www.mosek.com)
- [SketchyCGAL](https://github.com/alpyurtsever/SketchyCGAL)

# Remark
- Limited memory bundle method is only supported on Ubuntu.

# Installation
- To use SpectralSOS in Julia, run
```ruby
Pkg> add https://github.com/maihoanganh/SpectralSOS.git
```

# Usage
The following examples briefly guide to use SpectralSOS (see more examples in files test/test_....ipynb):

## Polynomial optimization
Consider an equality constrained POP on the unit sphere as follows:
```ruby
using DynamicPolynomials

@polyvar x[1:2] # variables

f=x[1]^2+0.5*x[1]*x[2]-0.25*x[2]^2+0.75*x[1]-0.3*x[2] # objective function to minimize

R=1.0 # squared radius of a sphere constraint
h=[R-sum(x.^2);(x[1]-1.0)*x[2]] # equality constraints (including the sphere constraint)

k=2 # relaxation order

using SpectralSOS

# get the optimal value and an optimal solution
opt_val,opt_sol = CTP_POP(x,f,h,k,R,method="LMBM",EigAlg="Arpack",tol=1e-5,scale=true)
```
Here:

- ```method="LMBM"```: [Limited memory bundle method](https://github.com/maihoanganh/LMBMinterface) (solver of spectral minimization). You can also choose ```method="PB"``` ([Proximal bundle method](https://github.com/maihoanganh/ProximalBundleMethod)) or ```method="SketchyCGAL"```.

- ```EigAlg="Arpack"```: [Arpack package](https://github.com/JuliaLinearAlgebra/Arpack.jl) (to compute the largest eigenvalue). You can also choose ```EigAlg="Normal"``` to use essential function of computing eigenvalues in Julia or ```EigAlg="Mix"``` to use the combination of the two methods.

- ```tol=1e-5```: the precision of the solver for spectral minimization.

- ```scale=true```: to scale the input of the SDP relaxation.

See other options in the [link](https://github.com/maihoanganh/SpectralSOS/blob/master/examples/test_random_dense_quadratic_on_sphere.ipynb).

The output from ```CTP_POP``` and other functions is explained in the [link](https://github.com/maihoanganh/SpectralSOS/blob/master/examples/test_random_dense_quadratic_on_sphere.ipynb).


To solve other types of POPs, see the links below:
- [Constrained POPs with single inequality (ball) constraint](https://github.com/maihoanganh/SpectralSOS/blob/master/examples/test_random_dense_QCQP_unique_inequality_(ball)_constraint.ipynb);
- [Constrained POPs on a ball](https://github.com/maihoanganh/SpectralSOS/blob/master/examples/test_random_dense_QCQP_on_ball.ipynb).



# References
For more details, please refer to:

**N. H. A. Mai, J.-B. Lasserre, and V. Magron. A hierarchy of spectral relaxations for polynomial optimization, 2020. Submitted.**

To get the paper's benchmarks, download the zip file in this [link](https://drive.google.com/file/d/11RqaDaXAngPAKSh-6RWPvVJILs6c95_Y/view?usp=sharing) and unzip the file.

The following commands allow one to reproduce the paper's benchmarks:
```ruby
using SpectralSOS

SpectralSOS.test_examples()

data="/home/hoanganh/Desktop/math-topics/SpectralPOP/codes/dataPOP" # path of data 

#The path needs to be changed on the user's computer

# Polynomial optimization
SpectralSOS.test_random_dense_quadratic_on_sphere(data) # Table 2
SpectralSOS.test_random_dense_equality_constrained_QCQP_on_sphere_first_order(data) # Table 3
SpectralSOS.test_random_dense_equality_constrained_QCQP_on_sphere_second_order(data) # Table 4 and 5
SpectralSOS.Evaluation_comparisons(data) # Table 6
SpectralSOS.test_random_dense_QCQP_unique_inequality_ball_constraint(data) # Table 7
SpectralSOS.test_random_dense_QCQP_on_ball(data) # Table 8
SpectralSOS.Norm_Subgrad(data) # Table 9
SpectralSOS.test_random_dense_quartics_on_sphere(data) # Table 10
SpectralSOS.test_deciding_nonnegativity(data) # Table 11
SpectralSOS.test_deciding_convexity(data) #Table 12
SpectralSOS.test_deciding_copositivity(data) #Table 13
```

On the display, two solvers are separated by a line: 

```~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```

and two problems are separated by three lines:

```~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```
```~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```
```~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```

# SpectralSOS
# SpectralPOP
