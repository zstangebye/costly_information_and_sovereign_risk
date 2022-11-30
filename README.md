# costly_information_and_sovereign_risk
Solution code for "Costly information and sovereign risk" with Grace Weishi Gu

Found in solution_code subfolder. Run the jobscript ier_ml.script or compile with MPI in the following order 

asa047.f (Legacy, open-source code for Nelder-Mead multivariate optimization)
global_vars.f90 (contains all global variables)
aux_functions.f90 (contains all relevant function calls)
main.f90 (main program, where the value-function iteration is located)

The code is designed with a modular structure to facilitate its translation to other dynamic non-linear models. The solution uses the mode of Gaussian Processes to approximate endogenous functions, such as value, price, and policy functions. To create a new function in this vein (let's call it function G), proceed as follows:

NEW VARIABLES (done in global_vars.f90)
0) Set lower and upper bounds for the function G: GL, GH
1) Create a 2-D global vector housing the current kernel coefficients for function G: G_kernCoeff
2) Create a global vector with length equal to the number of grid points to house the coefficients for the GPR approximation: G_gprCoeff
3) Create a global vector with length equal to the number of grid points to house the (logit-transformed) values of X at the grid points: tset_Gz
4) Create a pair of scalars to house the sample mean and standard deviation of tset_Gz: G_mean, G_std. These are used to clean, i.e., de-mean and standardize the logit-transformed output values of G, such that the mean-zero prior has no significant effect in interpolation
5) Create a global vector with length equal to the number of grid points to house the "cleaned" values of tset_Gz: tset_G. This is what the GPR will operate on

NEW FUNCTIONS (done in aux_functions.f90)
0) Create your new function to be called G(). Follow the template for vRepay or vDefault, i.e., translate input variables into the unit hypercube on which the GPR operates and evaluate the GPR at the current (a) kernel parameters, (b) GPR coefficients, (c) mean and standard-deviation scalars. Note GPR_approx takes the dimension of the problem as well as the hypercube-translated gridpoints as arguments 
1) Create the function/subroutine that delivers new values for G given other current equilibrium objects, e.g., Bellman equation optimization, prices based on expectations, etc. This will be called in the VFI to provide data to run the GPR on
2) Create a log-likelihood function tailored to your function G: negLL_G(). Follow the template for negLL_VR() or negLL_VD(). These eliminate extra arguments from the general negLL() function so as to be suitable to be passed to the Nelder-Mead minimization in the main VFI.

CHANGES IN VFI (done in main.f90)
0) Define locally a covariance matrix: K_G. This will be used to generate the (global) GPR coefficients
1) Normalize global mean and standard deviation parameters to zero and one respectively. Set G_gprCoeff to initial conditions (bearing in mind the logit-transform)
2) Add the update step for your function G in the VFI loop wherever it belongs. Follow the examples in the code for either new_q or whatever is nearest your application. This portion is automatically parallelized, with logit-transformed results going into the val1 vector, which recollected once all processors have executed.
3) Recollect val1 vector into tset_Gz
4) Clean tset_Gz by de-meaning/standardizing it using G_mean and G_std to create tset_G. Again, follow the template already there
5) Optimize the likelihood (negLL_G) of tset_G over the log-10 kernel parameters. This is done in a loop over all starting points with finite likelihoods. Again, follow the template for some other function, e.g., negLL_q
6) Use the optimal kernel parameters (G_kernCoeff) to create the covariance matrix K_G
7) Use the matrix K_G to create the GPR coefficients (G_gprCoeff). Once this is done, the function G() will be automatically updated and may be called from anywhere in the code, especially in the next loop.
8) That's it!
