! This module defines the (relatively few) global variables;
! The majority of variables are local in scope to facilitate efficiency

module global_vars  

    implicit none
    ! IF RUNNING IN PARALLEL, uncomment " include 'mpif.h' ", set mpi_on = .TRUE.,
    ! and comment the definition of the four MPI variables (which are also defined in mpif.h)

    include 'mpif.h' 
    
    logical, parameter :: mpi_on = .TRUE. 
    ! integer :: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_LOGICAL, MPI_INTEGER ! These variables are not used 
    ! if mpi is not on but have to be defined, nevertheless

    ! Preference parameters
    double precision, parameter :: ggamma = 2.0 ! Country risk-aversion 
    double precision, parameter :: CRRA_LEND = 0.97 ! Lender risk-aversion
    double precision, parameter :: bbeta = 0.963 ! Country impatience

    ! Output process parameters
    double precision, parameter :: rrho_z = 0.9212, ssigma_z = .0226 ! Country output process params (y = exp(z))
    double precision, parameter :: w = 1.0 ! Lender wealth

    ! Debt parameters 
    double precision, parameter :: r = .01 ! Quarterly risk-free rate
    double precision, parameter :: coup = r ! debt coupon (paid even when debt matures)
    double precision, parameter :: llambda = 0.066 - coup ! Debt maturity (DbtSrv-Coupon)

    ! minimum issuance price to prevent dilution before default (from HMS (2016))
    double precision, parameter :: min_q = .7 ! 70% of risk-free price

    ! Information cost 
    double precision, parameter :: kkappa = 1.5e-5 ! Information cost

    ! Default cost parameters
    double precision, parameter :: ythresh = .9 ! break below which default is free
    double precision, parameter :: MP = .225 ! maximum proportional cost (applied at yH due to standard curvature!)
    double precision, parameter :: pphi = 0.083 ! Probability of re-entry deal


    ! Recovery parameters
    double precision, parameter :: norm_mu_est = 0.79287 ! Asonuma (2012) # 0.66945 C&T (2013)
    double precision, parameter :: norm_sigma_est = 1.15429 ! Asonuma (2012) # 1.13600 C&T (2013)
    ! Asonuma implies mean of 0.651 and volatility of .214 in logit-normal
    ! C&T implies mean of 0.63 and volatility of .216 in logit-normal

    ! Information flow upper bound 
    double precision, parameter :: rho_cap = .99 ! Highest correlation
    double precision, parameter ::  info_cap = log( 1.0 / sqrt(1.0 - rho_cap**2.0)) ! Highest implied information

    ! "Size of model" parameters
    integer :: ccounter ! Iteration counter
    integer, parameter :: nDims = 3 ! The dimensionality of the problem
    integer, parameter :: nInputs = nDims*20 ! Number of training inputs (grid points)
    integer, parameter :: quad_size = 15 ! Size of quadrature used in expectations
    integer, parameter :: N_validate = 10000 ! Number of points to check in each iteration
    integer, parameter :: Nsimul = 500000 ! Length of simulation in each iteration
    integer, parameter :: max_iters = 600 ! Maximum number of iterations  
    double precision, parameter :: conv_tol = 1.0e-6 ! Convergence tolerance
    double precision, parameter :: opt_tol = 1.0e-9 ! Tolerance level for optimization problem (golden-section)
    double precision, parameter :: typ_thresh = 1.0e-5 ! Level of typicality required in grid point draws
    integer, parameter :: nMoments = 8 ! Number of simulated moments to report
    integer, parameter :: nKernCoeffsDraws = 96 ! Number of random kernel parameter coefficients to draw

    double precision, dimension(quad_size) :: GCweights, GCnodes ! Gaussian weights and nodes
    double precision, dimension(quad_size) :: m_GCnodes, x_GCnodes ! Weights and nodes

    double precision, parameter :: bL = 0.0, bH = 0.8 ! Bounds on debt levels 
    double precision, parameter :: VL = -100.0, VH = -1.0e-6 ! Bounds on value functions

    ! Debt grid measures
    ! Threshold for debt grid: bint_thresh uniform points on [bL,b_thresh] 
    ! and nInputs-bint_thresh uniform points on [bthresh,bH]
    double precision, parameter :: b_thresh = 0.6 
    integer, parameter :: bint_thresh = nDims*10

    ! Relevant constants
    double precision, parameter :: ppi = 3.14159265358979323846
    double precision, parameter :: e = 2.71828182845904523536 
    double precision, parameter :: stds_out = 3.0 ! How many stdevs out to set bounds for shocks
    double precision, parameter :: ssigma_uncond = ssigma_z/sqrt(1.0-rrho_z**2.0)
    double precision, parameter :: x_entropy =  .5*log(2.0*ppi*e*norm_sigma_est**2.0), & 
        z_entropy =  .5*log(2.0*ppi*e*ssigma_uncond**2.0)
    ! Bounds on output levels (determined by AR-process and stds_out)
    double precision, parameter :: yL=exp(-stds_out*ssigma_uncond) , yH = exp(stds_out*ssigma_uncond) 
    double precision, parameter :: mL = 0.0, mH = 1.0, xL = norm_mu_est - stds_out*norm_sigma_est,  & 
        xH = norm_mu_est + stds_out*norm_sigma_est
    ! Divide by this to ensure that normal distributions integrate to one 
    ! (Gauss-Chebyshev mass within 2*stds_out = 2*3 stdev range)
    double precision, parameter :: pdf_adjust = 0.997327028223715 


    ! Default cost parameters (as related to boundaries)
    double precision, parameter :: ylev = -MP*ythresh/(yH-ythresh) ! Default cost level parameter
    double precision, parameter :: ycurv = MP/(yH-ythresh) ! Default cost curvature parameter

    integer, parameter :: bN_opt = 50 ! This is the grid size over which we conduct global searches

    ! GPR Parameters
    double precision, parameter :: sn_star = 1.0e-4 ! Assumed measurement error (when not updating hyper-parameters)
    double precision, dimension(nDims,nInputs) :: gridPointsX, gridPointsM, scaledGridPointsX, scaledGridPointsM ! Input grids
    double precision, dimension(2) :: VD_kernCoeff, q_kernCoeff, qd_kernCoeff, &
        qRN_kernCoeff, qdRN_kernCoeff, def_kernCoeff ! GPR kernel parameters 
    double precision, dimension(2) :: VR_kernCoeff, A_kernCoeff, I_kernCoeff ! more GPR kernel parameters
    double precision, dimension(nInputs) :: VR_gprCoeff, VD_gprCoeff, A_gprCoeff, I_gprCoeff, &
        q_gprCoeff, qd_gprCoeff, qRN_gprCoeff, qdRN_gprCoeff, def_gprCoeff ! GPR coefficients 
    ! Output values: tset_...z will represent the non-cleaned value and tset_... will be the cleaned 
    ! value that the GPR acts on
    double precision, dimension(nInputs) :: tset_VR, tset_VD, tset_A, tset_I, tset_q, tset_qd, &
        tset_qRN, tset_qdRN, tset_def, tset_VRz, tset_VDz, tset_Az, tset_Iz, &
        tset_qz, tset_qdz, tset_qRNz, tset_qdRNz, tset_defz 
    ! These are used to clean the output values before applying the GPR
    double precision :: VR_mean, VD_mean, A_mean, I_mean, q_mean, qd_mean, qRN_mean, qdRN_mean, def_mean, & 
        VR_std, VD_std, A_std, I_std, q_std, qd_std, qRN_std, qdRN_std, def_std
    ! Hyperparameter bounds (in log-terms)
    double precision, parameter :: kernSSLogLB = -1.0, kernSSLogUB = 9.0, kernLogLB = -1.0, kernLogUB = 9.0
        
    end module global_vars