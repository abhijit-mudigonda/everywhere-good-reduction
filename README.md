This repository contains code associated to the paper _Quadratic Fields Admitting Elliptic Curves with Rational j-Invariant and Good Reduction Everywhere_, by [Benjamin Matschke](https://math.bu.edu/people/matschke/) and [Abhijit Mudigonda](https://abhijit-mudigonda.github.io/math/). 

##Quickstart

In order to obtain the numerical bounds in Corollary D and Corollary E of the paper, run the following steps after cloning this repository. 
1. Make sure you have [SageMath](https://www.sagemath.org/) installed and using a Python 3 kernel, as well as the Python modules [argh](https://pythonhosted.org/argh/) and [tqdm](https://tqdm.github.io/).
2. To compute lower and upper bounds on the constants from the paper, run 
`sage compute_const.sage computeConstant`
3. To run unit tests, run `sage compute_const.sage runTests`

##Other stuff

To compute upper and lower bounds on the constant in Theorem 7.3, run 
`sage kappa.sage`

To compute





