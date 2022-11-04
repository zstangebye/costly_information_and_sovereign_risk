!  Machine learning code to solve sovereign debt models
!  This specific code solves it for:
!  'Costly information in sovereign debt markets'

! It implements Gaussian-Process Dynamic Programming, similar to the algorithm 
! described in Scheidigger and Bilionis (2019) but with noteworthy differences 
! Adjustments to the Scheidigger and Bilionis (2019) algorithm can be found in
! Appendix A of the Gu and Stangebye (2022)

! The code here is parallelized and designed to be run in Fortran MPI

!  Authors:
!  Zach Stangebye and Grace Gu
!  Nov 2022
 
program main

    ! Import auxiliary functions for use in main program
    use aux_functions

    implicit none

    double precision, dimension(nDims-1,N_validate) :: convergence_points ! Inputs for checking convergence
    double precision, dimension(N_validate) :: v_old_points, v_new_points ! Old and new value functions for convergence

    double precision, dimension(bN_opt) :: b_opt_grid, v_opt_grid ! Grid for optimization in global step
    double precision :: best_v, bL_local, bH_local, minimized_LL ! place-holder for optimization in global step
    integer :: max_i ! place-holder for optimization in global step

    double precision :: ddist, time1, time2  ! Current sup-norm distance between value functions and times

    integer :: i, j, i1 ! Just indices for loops

    ! Temporary variables used in solution
    double precision :: ytod, bissued, btod, xtod, mtod, optInfoSoln, optInfoSolnLevel, optBorrowingSoln, bellmanValue
    double precision, dimension(nInputs,nInputs) :: K_VR, K_A, K_I, K_VD, K_q, K_qd, K_qRN, K_qdRN, K_def, eyeN, ttempK
    double precision, dimension(nKernCoeffsDraws*2) :: random_param_line, opt_param_line
    double precision, dimension(nKernCoeffsDraws,2) :: random_param_block, opt_param_block
    double precision, dimension(nKernCoeffsDraws) :: param_LL
    integer, dimension(1) :: maxLL_ind

    double precision, dimension(nDims) :: test_point

    ! Nelder-Mead minimization variables
    integer ::  icount, ifault, kcount, konvge, numres
    double precision :: reqmin
    double precision, dimension(2) :: step, temp_coeff_vec, temp_coeff_vec_soln, rand_param_draw
    logical :: found_hyper_params

    ! Parallelization variables 
    ! ierr reads error signals; id is processor number; nproc is total number of processors 
    ! nii is how many grid points are handed to each processor; itop and iend are gridpoints assigned 
    ! to current processor
    integer:: ierr, id, nproc, nii, itop, iend, nii2, itop2, iend2, nii3, itop3, iend3
    double precision, dimension(:), allocatable :: val1, val2, val3, val4, val5, & 
        val1agg, val2agg, val3agg, val4agg, val5agg

    ! Simulation variables
    double precision, dimension(nMoments) :: moment_vec
    double precision, dimension(1) :: seed
    integer :: status(MPI_STATUS_SIZE) 


    ! Begin by assigning different processors the indices in the group arrays that each will compute
    ! Start by initializing parallelization variables
    if(mpi_on) then
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,id,ierr)
        call mpi_comm_size(mpi_comm_world,nproc,ierr)
    else
        id = 0
        nproc = 1
    end if

    nii  = int(DBLE(nInputs-1)/DBLE(nproc))+1 ! Number of inputs handled by each processor (economic optimization)
    itop = id*nii+1 ! This processor's first index
    iend = min((id+1)*nii,nInputs) ! This processor's last index

    ! How many inputs each processor works with depends on how many processors we're using
    allocate(val1(nii), val1agg(nii*nproc), val2(nii), val2agg(nii*nproc), & 
        val3(nii), val3agg(nii*nproc))

    nii2  = int(DBLE(nKernCoeffsDraws-1)/DBLE(nproc))+1 ! Number of inputs handled by each processor (hyperparam optimization 1)
    itop2 = id*nii2+1 ! This processor's first index
    iend2 = min((id+1)*nii2,nKernCoeffsDraws) ! This processor's last index

    allocate(val4(nii2), val4agg(nii2*nproc))

    nii3  = int(DBLE(nKernCoeffsDraws*2-1)/DBLE(nproc))+1 ! Number of inputs handled by each processor (hyperparam optimization 2)
    itop3 = id*nii3+1 ! This processor's first index
    iend3 = min((id+1)*nii3,nKernCoeffsDraws*2) ! This processor's last index

    allocate(val5(nii3), val5agg(nii3*nproc))
    ! Now the parallel variables are initiated!

    ! First, we'll set the technical Nelder-Mead variables
    reqmin = 1.0D-06

    step(1) = 1.0D+00
    step(2) = 1.0D+00

    konvge = 10
    kcount = 500

    ! Initial kernel parameters meant only to keep NaNs away
    ! since GPR coefficients will be zero
    VR_kernCoeff = 1.0 
    A_kernCoeff = 1.0
    I_kernCoeff = 1.0
    VD_kernCoeff = 1.0
    q_kernCoeff = 1.0 
    qd_kernCoeff = 1.0
    qRN_kernCoeff = 1.0 
    qdRN_kernCoeff = 1.0

    ! We  initialize the relevant GPR coefficients. The only ones used before optimization 
    ! are VR and VD. We set large values to correspond to zero (the upper bound) 
    ! after passing through the logit-transform 

    VR_gprCoeff = 500.0
    A_gprCoeff = 0.0
    I_gprCoeff = 0.0
    VD_gprCoeff = 500.0
    q_gprCoeff = 0.0 
    qd_gprCoeff = 0.0
    qRN_gprCoeff = 0.0 
    qdRN_gprCoeff = 0.0

    ! Data cleaning parameters. Initialized at standard normal and updated 
    ! every GPR
    VR_mean = 0.0
    VD_mean = 0.0
    A_mean = 0.0
    I_mean = 0.0
    q_mean = 0.0
    qd_mean = 0.0
    qRN_mean = 0.0
    qdRN_mean = 0.0

    VR_std = 1.0
    VD_std = 1.0
    A_std = 1.0
    I_std = 1.0
    q_std = 1.0
    qd_std = 1.0
    qRN_std = 1.0
    qdRN_std = 1.0

    ! Initialize the quadrature and draw the gridpoints
    call GaussChebyshevQuad(GCweights,GCnodes)

    ! Notice that we do not separately seed the random processes. Thus they will give every processor the same
    ! grid points
    call drawInputs(gridPointsX,gridPointsM)
    scaledGridPointsX(1,:) = (gridPointsX(1,:) - yL)/(yH-yL)
    scaledGridPointsM(1,:) = (gridPointsM(1,:) - yL)/(yH-yL)
    scaledGridPointsX(2,:) = (gridPointsX(2,:) - bL)/(bH-bL)
    scaledGridPointsM(2,:) = (gridPointsM(2,:) - bL)/(bH-bL)
    scaledGridPointsX(3,:) = (gridPointsX(3,:) - xL)/(xH-xL)
    scaledGridPointsM(3,:) = (gridPointsM(3,:) - mL)/(mH-mL)

    ! Printed grid points are SCALED into [0,1] so that matlab code can directly create the GP functions 
    ! for plotting and analysis
    if(id .EQ. 0) then 
        open(200,file="grid_points_X.txt")
        
            do i=1,nInputs
                write(200,"(3F9.6)") scaledGridPointsX(:,i)
            end do 

        close(200)

        open(200,file="grid_points_M.txt")
        
            do i=1,nInputs
                write(200,"(3F9.6)") scaledGridPointsM(:,i)
            end do 

        close(200)

        open(200,file="bounds.txt") 
        
            write(200,*) bL, bH 
            write(200,*) yL, yH 
            write(200,*) xL, xH 
            write(200,*) mL, mH
            write(200,*) VL, VH

        close(200)
    end if 

    ! Give each processor a different seed now to ensure randomness when needed
    call random_number(seed) 
    seed = int(seed*1000)*(id+1) 
    call random_seed(put=int(seed))

    ! Now we set the initial values of the hyperparameters to be optimized over for every function: One set of initial values for 
    ! every processor. These will be the same throughout. As such, we select initial values that imply finite likelihoods. This ensures
    ! stability when paired with a Simple/Nelder-Mead optimization algorithm for the likelihood maximization

    tset_I = 1.0 ! We pick an arbitrary vector of outputs to set to unity for this routine

    ! Setting these to one implies that calling their likelihood just delivers log(det(Ssigma)) + 1'*Ssigma^(-1)*1, which depends only 
    ! on the hyperparameters we're evaluating here.

    ! This whole procedure is done by one processor and then shared with the other processors

    if( id .EQ. 0) then 

        call random_seed(put= (/ 1284 /) )

        do i = 1, nKernCoeffsDraws 

            found_hyper_params = .FALSE.
            
            do while( .NOT. found_hyper_params)

                call random_number(rand_param_draw)
                temp_coeff_vec(1) = kernSSLogLB + rand_param_draw(1)*(kernSSLogUB-kernSSLogLB)
                temp_coeff_vec(2) = kernLogLB + rand_param_draw(2)*(kernLogUB-kernLogLB) 

                if( negLL_info(temp_coeff_vec) < 9.9e9 ) then 
                    found_hyper_params = .TRUE.
                end if 

            end do 
            val4agg( i ) = negLL_info(temp_coeff_vec)
            val5agg( (i-1)*2+1 : i*2 ) = temp_coeff_vec

        end do 

        do i=1,nproc-1
            call MPI_SEND( val4agg, nKernCoeffsDraws, MPI_DOUBLE_PRECISION, i, &
                        i, MPI_COMM_WORLD, ierr )
            
            call MPI_SEND( val5agg, nKernCoeffsDraws*2, MPI_DOUBLE_PRECISION, i, &
                        i*10, MPI_COMM_WORLD, ierr )
        end do 
        
    else 

        call MPI_RECV(val4agg, nKernCoeffsDraws, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )  
        
        call MPI_RECV(val5agg, nKernCoeffsDraws*2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )  

    end if 

    param_LL = val4agg(1:nKernCoeffsDraws)
    random_param_line = val5agg(1:nKernCoeffsDraws*2)
  

    ! Now reformulate the random draws into blocks
    i1 = 1
    do i=1,nKernCoeffsDraws 
        do j=1,2 
            random_param_block(i,j) = random_param_line(i1)
            i1 = i1 + 1
        end do 
    end do 

    ! Here, we create a coarse grid for optimization; we perform a continuous
    ! optimization on the region highlighted as optimal by this grid search
    do i=1,bN_opt 
        b_opt_grid(i) = bL + (bH-bL)/(DBLE(bN_opt)-1.0)*(DBLE(i)-1.0)
    end do

    ! Here, we print the likelihoods of these various kernels given unity vectors 
    ! to check that our routine worked
    if(id .EQ. 0) then 
        open(3,file="init_kern.txt") 

            write(3,*) param_LL

        close(3)
        
        open(3,file="param_draws_and_parallel_architecture.txt")

            write(3,*) nii 
            write(3,*) nii2 
            write(3,*) nii3
            write(3,*) itop, iend 
            write(3,*) itop2, iend2  
            write(3,*) itop3, iend3

            do i=1,nKernCoeffsDraws
                write(3,*) random_param_block(i,:)
            end do 

            write(3,*) random_param_line

        close(3) 

        open(3,file="b_opt_grid.txt")

            write(3,*) b_opt_grid

        close(3)

    end if 

    ! Double-check params are same across processors
    if(id .EQ. 2) then 

        open(3,file="init_kern2.txt") 

            write(3,*) random_param_block

        close(3)

    end if 

    ! Now every processor has the parameter block! 
              
    ! The Gaussian quadrature nodes scaled to the region of m and x
    m_GCnodes =  mL + ( 1.0 + GCnodes)/2.0 *(mH-mL)
    x_GCnodes =  xL + ( 1.0 + GCnodes)/2.0 *(xH-xL)

    ! Make a quick identity matrix to be used in updating the GPR coefficients
    eyeN = 0.0
    do i=1,nInputs 
        eyeN(i,i) = 1.0
    end do 

    ! Here begins the actual iterative procedure
    ddist = 1.0 
    ccounter = 0
    do while((ddist > conv_tol) .AND. (ccounter < max_iters))

        call cpu_time(time1) ! Time the iterations

        ! Start by drawing test points initially for convergence check
        if(id .EQ. 0) then 
            call drawConvergenceInputs(convergence_points)
            ! convergence_points = gridPointsX(1:2,:)
            do i=1,N_validate 
                v_old_points(i) = vRepay(convergence_points(:,i))
            end do 

            open(3,file="convergence_points.txt")

                do i=1,N_validate
                    write(3,*) convergence_points(:,i)
                end do 
                    
            close(3)
        end if 

        ccounter = ccounter + 1

        ! We begin the backward iterations by updating the information acquisition function;
        ! Notice we do not need to do this (or the pricing function) in the final period

            ! First, we solve the forecaster's problem at the training inputs
            do i=itop,iend 

                ! This problem only operates from penultimate period on
                if(ccounter > 1) then 
                    call localMaximizeFun(DBLE(1.0e-6),rho_cap ,optInfoSoln , optInfoSolnLevel, gridpointsX(:,i), 1) 
                    ! Last "1" indicates that we're solving forecaster's problem

                else 
                    optInfoSoln = 1.0e-6
                end if 

                ! Recall we need to do a logit-transform on the GP since this is a bounded function
                val1(i-itop+1) = log(optInfoSoln/(1.0-optInfoSoln))  

            end do 

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val1(1),nii,mpi_double_precision, &
                    val1agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                tset_Iz = val1agg(1:nInputs)

            else
                tset_Iz = val1
            end if

            ! Next, we "clean" the point estimates (de-mean and standardize)
            I_mean = sum(tset_Iz)/DBLE(nInputs)
            I_std = sqrt( sum((tset_Iz - I_mean)**2.0)/DBLE(nInputs-1) )
            if(I_std < 1.0e-6) I_std = 1.0 
            tset_I = (tset_Iz - I_mean)/I_std

            ! Next we update the Gaussian process (this automatically changes the info policy function itself)

            ! This bit optimizes the hyper-parameters
            ! by doing a parallelized search over hyper-parameters
            do i=itop2,iend2 
                temp_coeff_vec = random_param_block(i,:)

                call nelmin(negLL_info,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                    konvge, kcount, icount, numres, ifault ) 

                ! Store the likelihood 
                val4(i-itop2+1) = minimized_LL

                ! And store the parameters
                val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln
            end do

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                    val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
                param_LL = val4agg(1:nKernCoeffsDraws)

                call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                    val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

                opt_param_line = val5agg(1:nKernCoeffsDraws*2)

            else
                param_LL = val4
                opt_param_line = val5
            end if


            ! Now reformulate the optimal hyperparameters into blocks
            do i=1,nKernCoeffsDraws 
                opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
            end do 

            ! Select the hyperparameters with the biggest likelihood
            maxLL_ind = minloc(param_LL)
            temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

            temp_coeff_vec = random_param_block(maxLL_ind(1),:)

            I_kernCoeff = 10.0**temp_coeff_vec_soln

            ! Update the gpr coefficients, which automatically updates the policy itself
            call createK(K_I, I_kernCoeff,1)
            I_gprCoeff = matmul(inv(K_I + sn_star**2.0*eyeN),tset_I)

            ! The info policy function is now updated! So print the results

            if(id .EQ. 0) then 

                open(3,file="tset_I.txt")

                    write(3,*) 1.0/(1.0+exp(-tset_Iz))

                close(3)

                open(3,file="I_kern_info.txt")

                    write(3,*) param_LL
                    write(3,*) maxLL_ind(1)
                    write(3,*) temp_coeff_vec
                    write(3,*) temp_coeff_vec_soln

                close(3)

            end if 

            ! Next, compute expected default probabilites tomorrow. This is not needed for the solution 
            ! but we examine it in the analysis
            do i=itop,iend 

                ! First the repayment price (Y_t, B_{t+1}, x_t)
                ytod = gridPointsX(1,i)
                bissued = gridPointsX(2,i)
                xtod = gridPointsX(3,i)

                if( ccounter > 1) then 
                    val1(i-itop+1) = -log(1.0/max(defProbTomorrow(ytod,bissued,xtod), 1.0e-6) - 1.0)
                else 
                    val1(i-itop+1) = -log(1.0/(1.0e-6) - 1.0)
                end if 

            end do 

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val1(1),nii,mpi_double_precision, &
                    val1agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                tset_defz = val1agg(1:nInputs)

            else
                tset_defz = val1
            end if
            ! Next, we "clean" the point estimates (de-mean and standardize)
            def_mean = sum(tset_defz)/DBLE(nInputs)
            def_std = sqrt( sum((tset_defz - def_mean)**2.0)/DBLE(nInputs-1) )
            if(def_std < 1.0e-6) def_std = 1.0 
            tset_def = (tset_defz - def_mean)/def_std

            
            do i=itop2,iend2 
                temp_coeff_vec = random_param_block(i,:)

                call nelmin(negLL_def,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                    konvge, kcount, icount, numres, ifault ) 

                val4(i-itop2+1) = minimized_LL

                val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln
            end do

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                    val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
                param_LL = val4agg(1:nKernCoeffsDraws)

                call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                    val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

                opt_param_line = val5agg(1:nKernCoeffsDraws*2)

            else
                param_LL = val4
                opt_param_line = val5
            end if

            ! Select the hyperparameters with the biggest likelihood
            maxLL_ind = minloc(param_LL)

            do i=1,nKernCoeffsDraws 
                opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
            end do 

            temp_coeff_vec = random_param_block(maxLL_ind(1),:)
            temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

            def_kernCoeff = 10.0**temp_coeff_vec_soln

            ! Update the gpr coefficients, which automatically updates the policy itself
            call createK(K_def, def_kernCoeff,2)
            def_gprCoeff = matmul(inv(K_def + sn_star**2.0*eyeN),tset_def)

            if(id .EQ. 0) then  

                open(3,file="tset_def.txt")

                    ! write(3,*) tset_q
                    write(3,*) 1.0/(1.0 + exp(-tset_defz))

                close(3)

                open(3,file="def_kern_info.txt")

                    write(3,*) param_LL
                    write(3,*) maxLL_ind(1)
                    write(3,*) temp_coeff_vec
                    write(3,*) temp_coeff_vec_soln
                    write(3,*) negLL_def(random_param_block(maxLL_ind(1),:))
                    write(3,*) negLL_def(opt_param_block(maxLL_ind(1),:))

                close(3)

            end if 


            ! Next, we solve the pricing equations (default and repayment)
            do i=itop,iend 

                ! First the repayment price (Y_t, B_{t+1}, x_t)
                ytod = gridPointsX(1,i)
                bissued = gridPointsX(2,i)
                xtod = gridPointsX(3,i)

                if( ccounter > 1) then 
                    val1(i-itop+1) = -log(1.0/max(new_q(ytod,bissued,xtod), 1.0e-6) - 1.0)
                else 
                    val1(i-itop+1) = -log(1.0/(1.0e-6) - 1.0)
                end if 

                ! Now the default price (Y_t, B_t, m_t)
                btod = bissued
                mtod = gridPointsM(3,i)

                if( ccounter > 1) then 
                    val2(i-itop+1) = -log(1.0/max(new_qd(ytod,btod,mtod),1.0e-6) - 1.0)
                else 
                    val2(i-itop+1) = -log(1.0/(1.0e-6) - 1.0)
                end if 

            end do 

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val1(1),nii,mpi_double_precision, &
                    val1agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                tset_qz = val1agg(1:nInputs)

                call mpi_allgather(val2(1),nii,mpi_double_precision, &
                    val2agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                tset_qdz = val2agg(1:nInputs)

            else
                tset_qz = val1
                tset_qdz = val2
            end if

            ! Next, we "clean" the point estimates (de-mean and standardize)
            q_mean = sum(tset_qz)/DBLE(nInputs)
            q_std = sqrt( sum((tset_qz - q_mean)**2.0)/DBLE(nInputs-1) )
            if(q_std < 1.0e-6) q_std = 1.0 
            tset_q = (tset_qz - q_mean)/q_std

            qd_mean = sum(tset_qdz)/DBLE(nInputs)
            qd_std = sqrt( sum((tset_qdz - qd_mean)**2.0)/DBLE(nInputs-1) )
            if(qd_std < 1.0e-6) qd_std = 1.0 
            tset_qd = (tset_qdz - qd_mean)/qd_std

            ! Next we update the Gaussian process (this automatically changes the info policy function itself)
            ! This bit optimizes the hyper-parameters
            ! by doing a parallelized search over hyper-parameters

            do i=itop2,iend2 
                temp_coeff_vec = random_param_block(i,:)

                call nelmin(negLL_q,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                    konvge, kcount, icount, numres, ifault ) 

                val4(i-itop2+1) = minimized_LL

                val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln
            end do

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                    val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
                param_LL = val4agg(1:nKernCoeffsDraws)

                call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                    val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

                opt_param_line = val5agg(1:nKernCoeffsDraws*2)

            else
                param_LL = val4
                opt_param_line = val5
            end if

            ! Select the hyperparameters with the biggest likelihood
            maxLL_ind = minloc(param_LL)

            do i=1,nKernCoeffsDraws 
                opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
            end do 

            temp_coeff_vec = random_param_block(maxLL_ind(1),:)
            temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

            q_kernCoeff = 10.0**temp_coeff_vec_soln

            ! Update the gpr coefficients, which automatically updates the policy itself
            call createK(K_q, q_kernCoeff,2)
            q_gprCoeff = matmul(inv(K_q + sn_star**2.0*eyeN),tset_q)

            if(id .EQ. 0) then  

                open(3,file="tset_q.txt")

                    ! write(3,*) tset_q
                    write(3,*) 1.0/(1.0 + exp(-tset_qz))

                close(3)

                open(3,file="q_kern_info.txt")

                    write(3,*) param_LL
                    write(3,*) maxLL_ind(1)
                    write(3,*) temp_coeff_vec
                    write(3,*) temp_coeff_vec_soln
                    write(3,*) negLL_q(random_param_block(maxLL_ind(1),:))
                    write(3,*) negLL_q(opt_param_block(maxLL_ind(1),:))

                close(3)

            end if 
        
            do i=itop2,iend2 

                temp_coeff_vec = random_param_block(i,:)

                call nelmin(negLL_qd,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                    konvge, kcount, icount, numres, ifault ) 

                val4(i-itop2+1) = minimized_LL
                val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln

            end do

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                    val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
                param_LL = val4agg(1:nKernCoeffsDraws)

                call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                    val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

                opt_param_line = val5agg(1:nKernCoeffsDraws*2)

            else
                param_LL = val4
                opt_param_line = val5
            end if

            ! Select the hyperparameters with the biggest likelihood
            maxLL_ind = minloc(param_LL)

            do i=1,nKernCoeffsDraws 
                opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
            end do 

            temp_coeff_vec = random_param_block(maxLL_ind(1),:)
            temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

            qd_kernCoeff = 10.0**temp_coeff_vec_soln

            ! Update the gpr coefficients, which automatically updates the policy itself
            call createK(K_qd, qd_kernCoeff,3)
            qd_gprCoeff = matmul(inv(K_qd + sn_star**2.0*eyeN),tset_qd)

            if(id .EQ. 0) then 

                open(3,file="qd_kern_info.txt")
                    write(3,*) param_LL
                    write(3,*) maxLL_ind(1)
                    write(3,*) temp_coeff_vec
                    write(3,*) temp_coeff_vec_soln

                close(3)

                open(3,file="tset_qd.txt")

                    write(3,*) 1.0/(1.0+exp(-tset_qdz))
                    ! write(3,*) tset_qd

                close(3)

            end if 

            ! Next, we solve the RN pricing equations (default and repayment)
            do i=itop,iend 

                ! First the repayment price (Y_t, B_{t+1}, x_t)
                ytod = gridPointsX(1,i)
                bissued = gridPointsX(2,i)
                xtod = gridPointsX(3,i)

                if(ccounter > 1) then 
                    val1(i-itop+1) = -log(1.0/max(new_qRN(ytod,bissued,xtod),1.0e-6) - 1.0)
                else 
                    val1(i-itop+1) = -log(1.0/(1.0e-6) - 1.0)
                end if 

                ! Now the default price (Y_t, B_t, m_t)
                btod = bissued
                mtod = gridPointsM(3,i)

                if(ccounter > 1) then 
                    val2(i-itop+1) = -log(1.0/max(new_qdRN(ytod,btod,mtod),1.0e-6) - 1.0)
                else 
                    val2(i-itop+1) = -log(1.0/(1.0e-6) - 1.0)
                end if 

            end do 

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val1(1),nii,mpi_double_precision, &
                    val1agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                tset_qRNz = val1agg(1:nInputs)

                call mpi_allgather(val2(1),nii,mpi_double_precision, &
                    val2agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
                tset_qdRNz = val2agg(1:nInputs)

            else
                tset_qRNz = val1
                tset_qdRNz = val2
            end if

            ! Next, we "clean" the point estimates (de-mean and standardize)
            qRN_mean = sum(tset_qRNz)/DBLE(nInputs)
            qRN_std = sqrt( sum((tset_qRNz - qRN_mean)**2.0)/DBLE(nInputs-1) )
            if(qRN_std < 1.0e-6) qRN_std = 1.0 
            tset_qRN = (tset_qRNz - qRN_mean)/qRN_std

            qdRN_mean = sum(tset_qdRNz)/DBLE(nInputs)
            qdRN_std = sqrt( sum((tset_qdRNz - qdRN_mean)**2.0)/DBLE(nInputs-1) )
            if(qdRN_std < 1.0e-6) qdRN_std = 1.0 
            tset_qdRN = (tset_qdRNz - qdRN_mean)/qdRN_std

            ! Next we update the Gaussian process (this automatically changes the info policy function itself)
            ! This bit optimizes the hyper-parameters
            ! by doing a parallelized search over hyper-parameters

            do i=itop2,iend2 
                temp_coeff_vec = random_param_block(i,:)

                call nelmin(negLL_qRN,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                    konvge, kcount, icount, numres, ifault ) 

                val4(i-itop2+1) = minimized_LL

                val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln
            end do

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                    val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
                param_LL = val4agg(1:nKernCoeffsDraws)

                call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                    val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

                opt_param_line = val5agg(1:nKernCoeffsDraws*2)

            else
                param_LL = val4
                opt_param_line = val5
            end if

            ! Select the hyperparameters with the biggest likelihood
            maxLL_ind = minloc(param_LL)

            do i=1,nKernCoeffsDraws 
                opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
            end do 

            temp_coeff_vec = random_param_block(maxLL_ind(1),:)
            temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

            qRN_kernCoeff = 10.0**temp_coeff_vec_soln

            ! Update the gpr coefficients, which automatically updates the policy itself
            call createK(K_qRN, qRN_kernCoeff,2)
            qRN_gprCoeff = matmul(inv(K_qRN + sn_star**2.0*eyeN),tset_qRN)

            if(id .EQ. 0) then 

                open(3,file="tset_qRN.txt")

                    write(3,*) 1.0/(1.0 + exp(-tset_qRNz))

                close(3)

                open(3,file="qRN_kern_info.txt")

                    write(3,*) param_LL
                    write(3,*) maxLL_ind(1)
                    write(3,*) temp_coeff_vec
                    write(3,*) temp_coeff_vec_soln
                    write(3,*) negLL_qRN(random_param_block(maxLL_ind(1),:))
                    write(3,*) negLL_qRN(opt_param_block(maxLL_ind(1),:))

                close(3)

            end if 
        
            do i=itop2,iend2 

                temp_coeff_vec = random_param_block(i,:)

                call nelmin(negLL_qdRN,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                    konvge, kcount, icount, numres, ifault ) 

                val4(i-itop2+1) = minimized_LL
                val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln

            end do

            ! Make sure all processes stop to wait here in parallelization; then collect all results
            if(mpi_on) then
                call MPI_BARRIER(mpi_comm_world, ierr)
            end if

            ! Collect all processors' activity together
            if(mpi_on) then

                call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                    val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
                param_LL = val4agg(1:nKernCoeffsDraws)

                call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                    val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

                opt_param_line = val5agg(1:nKernCoeffsDraws*2)

            else
                param_LL = val4
                opt_param_line = val5
            end if


            ! Select the hyperparameters with the biggest likelihood
            maxLL_ind = minloc(param_LL)

            do i=1,nKernCoeffsDraws 
                opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
            end do 

            temp_coeff_vec = random_param_block(maxLL_ind(1),:)
            temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

            qdRN_kernCoeff = 10.0**temp_coeff_vec_soln

            ! Update the gpr coefficients, which automatically updates the policy itself
            call createK(K_qdRN, qdRN_kernCoeff,3)
            qdRN_gprCoeff = matmul(inv(K_qdRN + sn_star**2.0*eyeN),tset_qdRN)

            if(id .EQ. 0) then 

                open(3,file="qdRN_kern_info.txt")
                    write(3,*) param_LL
                    write(3,*) maxLL_ind(1)
                    write(3,*) temp_coeff_vec
                    write(3,*) temp_coeff_vec_soln

                close(3)

                open(3,file="tset_qdRN.txt")

                    write(3,*) 1.0/(1.0+exp(-tset_qdRNz))
                    ! write(3,*) tset_qdRN

                close(3)

            end if 


        ! Finally, we solve the sovereign Bellman equations equations (default and repayment)
        do i=itop,iend 

            ! First the repayment value (Y_t, B_{t})

            ! First, conduct a coarse global grid search before a fine, non-grid, local one
            best_v = -1.0e10 
            max_i = 1
            do i1=1,bN_opt 
                v_opt_grid(i1) = repayBellman(b_opt_grid(i1),gridpointsX(1,i),gridPointsX(2,i))

                if(v_opt_grid(i1) .GE. best_v) then 
                    best_v = v_opt_grid(i1) 
                    max_i = i1 
                end if 

            end do 

            if(max_i .EQ. 1) then 
                bL_local = bL
                bH_local = b_opt_grid(2)
            elseif(max_i .EQ. bN_opt) then 
                bL_local = b_opt_grid(bN_opt-1)
                bH_local = b_opt_grid(bN_opt)
            else 
                bL_local = b_opt_grid(max_i-1)
                bH_local = b_opt_grid(max_i+1)
            end if 

            ! Here, we check to ensure our windows never box us into an area that is outside 
            ! of the constraints for all search directions. If it is, restrain it to the area of 
            ! allowable issuance
            test_point = gridpointsX(:,i)
            test_point(2) = (1.0-llambda)*gridpointsX(2,i)
            if( (bH_local > test_point(2)) .AND. (q(test_point) < min_q) ) then 
                bH_local = test_point(2) 
            end if 

            call localMaximizeFun(bL_local,bH_local,optBorrowingSoln, bellmanValue, gridpointsX(:,i), 2) 
            ! Last "2" indicates that we're solving Bellman problem

            ! ! Need to check zero points with price constraint binding to ensure policy continuity
            test_point = gridpointsX(:,i)
            test_point(2) = optBorrowingSoln
            if( (q(test_point) < min_q) .AND. (gridpointsX(2,i) < 1.0e-6 ) ) then 
                optBorrowingSoln = bL + 1.0e-6
            end if

            ! ! Here, we've violated the issuance at low prices constraint
            if( (bellmanValue .ge. VH) .OR. (bellmanValue .LE. VL) ) then 

                optBorrowingSoln = (1.0-llambda)*gridPointsX(2,i)-1.0e-6
                bellmanValue = repayBellman(optBorrowingSoln,gridpointsX(1,i),gridPointsX(2,i))
                 
            end if 

            ! Recall we need to do a logit-transform on the policy function GP since this is a bounded function
            val1(i-itop+1) = -log( (bH - bL)/min(max(optBorrowingSoln - bL,1.0e-6),bH-bL-1.0e-6) - 1.0)

            ! logit-transformation not needed for value functions, but we'll do it for consistency (with wide bounds)
            val2(i-itop+1) = -log( (VH-VL)/( bellmanValue - VL) - 1.0)

            ! Now the default value (Y_t, B_t, m_t)
            ytod = gridPointsM(1,i)
            btod = gridPointsM(2,i)
            mtod = gridPointsM(3,i)

            val3(i-itop+1) = -log( (VH-VL)/(defaultBellman(ytod,btod,mtod) - VL) - 1.0)

        end do 

        ! Make sure all processes stop to wait here in parallelization; then collect all results
        if(mpi_on) then
            call MPI_BARRIER(mpi_comm_world, ierr)
        end if

        ! Collect all processors' activity together
        if(mpi_on) then

            call mpi_allgather(val1(1),nii,mpi_double_precision, &
                val1agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
            tset_Az = val1agg(1:nInputs)

            call mpi_allgather(val2(1),nii,mpi_double_precision, &
                val2agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
            tset_VRz = val2agg(1:nInputs)

            call mpi_allgather(val3(1),nii,mpi_double_precision, &
                val3agg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
            tset_VDz = val3agg(1:nInputs)

        else
            tset_A = val1
            tset_VR = val2
            tset_VD = val3
        end if

        ! Next, we "clean" the point estimates (de-mean and standardize)
        A_mean = sum(tset_Az)/DBLE(nInputs)
        A_std = sqrt( sum((tset_Az - A_mean)**2.0)/DBLE(nInputs-1) )
        if(A_std < 1.0e-6) A_std = 1.0 
        tset_A = (tset_Az - A_mean)/A_std

        VR_mean = sum(tset_VRz)/DBLE(nInputs)
        VR_std = sqrt( sum((tset_VRz - VR_mean)**2.0)/DBLE(nInputs-1) )
        if(VR_std < 1.0e-6) VR_std = 1.0 
        tset_VR = (tset_VRz - VR_mean)/VR_std

        VD_mean = sum(tset_VDz)/DBLE(nInputs)
        VD_std = sqrt( sum((tset_VDz - VD_mean)**2.0)/DBLE(nInputs-1) )
        if(VD_std < 1.0e-6) VD_std = 1.0 
        tset_VD = (tset_VDz - VD_mean)/VD_std

        do i=itop2,iend2 

            temp_coeff_vec = random_param_block(i,:)

            call nelmin(negLL_VD,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                  konvge, kcount, icount, numres, ifault ) 

            val4(i-itop2+1) = minimized_LL
            val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln
        end do

        ! Make sure all processes stop to wait here in parallelization; then collect all results
        if(mpi_on) then
            call MPI_BARRIER(mpi_comm_world, ierr)
        end if

        ! Collect all processors' activity together
        if(mpi_on) then

            call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
            param_LL = val4agg(1:nKernCoeffsDraws)

            call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

            opt_param_line = val5agg(1:nKernCoeffsDraws*2)

        else
            param_LL = val4
            opt_param_line = val5
        end if

        ! Select the hyperparameters with the biggest likelihood
        maxLL_ind = minloc(param_LL)

        do i=1,nKernCoeffsDraws 
            opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
        end do 

        temp_coeff_vec = random_param_block(maxLL_ind(1),:)
        temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)

        VD_kernCoeff = 10.0**temp_coeff_vec_soln

        ! Update the gpr coefficients, which automatically updates the policy itself
        call createK(K_VD, VD_kernCoeff,3)
        VD_gprCoeff = matmul(inv(K_VD + sn_star**2.0*eyeN),tset_VD)

        if(id .EQ. 0) then
            
            open(3,file="VD_kern_info.txt")
                write(3,*) param_LL
                write(3,*) maxLL_ind(1)
                write(3,*) temp_coeff_vec
                write(3,*) temp_coeff_vec_soln
                write(3,*) negLL_VD(random_param_block(maxLL_ind(1),:))
                write(3,*) negLL_VD(opt_param_block(maxLL_ind(1),:))

                do i=1,nInputs
                    test_point(1) = 1.0
                    test_point(2) = 0.2
                    test_point(3) = DBLE(i/nInputs)
                    write(3,*) vDefault( test_point )
                end do 

            close(3)

        end if 

        do i=itop2,iend2 
            temp_coeff_vec = random_param_block(i,:)

            call nelmin(negLL_VR,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                  konvge, kcount, icount, numres, ifault ) 

            val4(i-itop2+1) = minimized_LL
            val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln 
        end do

        ! Make sure all processes stop to wait here in parallelization; then collect all results
        if(mpi_on) then
            call MPI_BARRIER(mpi_comm_world, ierr)
        end if

        ! Collect all processors' activity together
        if(mpi_on) then

            call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
            param_LL = val4agg(1:nKernCoeffsDraws)

            call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

            opt_param_line = val5agg(1:nKernCoeffsDraws*2)

        else
            param_LL = val4
            opt_param_line = val5
        end if

        do i=1,nKernCoeffsDraws 
            opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
        end do 

        ! Select the hyperparameters with the biggest likelihood
        
        maxLL_ind = minloc(param_LL)
        temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)
        temp_coeff_vec = random_param_block(maxLL_ind(1),:)

        VR_kernCoeff = 10.0**temp_coeff_vec_soln

        ! Update the gpr coefficients, which automatically updates the policy itself
        call createK(K_VR, VR_kernCoeff,1)
        VR_gprCoeff = matmul(inv(K_VR + sn_star**2.0*eyeN),tset_VR)

        if(id .EQ. 0) then 

            open(3,file="VR_kern_info.txt")
                write(3,*) param_LL
                write(3,*) maxLL_ind(1)
                write(3,*) temp_coeff_vec
                write(3,*) temp_coeff_vec_soln

                write(3,*) negLL_VR(temp_coeff_vec)
                write(3,*) negLL_VR(temp_coeff_vec_soln)

                do i=1,nInputs 
                    test_point(1) = gridPointsX(1,1)
                    test_point(2) = gridPointsX(2,i)
                    write(3,*) vRepay( test_point(1:2) )
                end do 

            close(3)

        end if 

        do i=itop2,iend2 
            temp_coeff_vec = random_param_block(i,:)

            call nelmin(negLL_A,2,temp_coeff_vec, temp_coeff_vec_soln, minimized_LL, reqmin, step, &
                  konvge, kcount, icount, numres, ifault ) 

            val4(i-itop2+1) = minimized_LL
            val5( (i-itop2)*2+1 : (i-itop2 + 1)*2 ) = temp_coeff_vec_soln
        end do

        ! Make sure all processes stop to wait here in parallelization; then collect all results
        if(mpi_on) then
            call MPI_BARRIER(mpi_comm_world, ierr)
        end if

        ! Collect all processors' activity together
        if(mpi_on) then

            call mpi_allgather(val4(1),nii2,mpi_double_precision, &
                val4agg(1),nii2,mpi_double_precision,mpi_comm_world,ierr)
            param_LL = val4agg(1:nKernCoeffsDraws)

            call mpi_allgather(val5(1),nii3,mpi_double_precision, &
                val5agg(1),nii3,mpi_double_precision,mpi_comm_world,ierr)

            opt_param_line = val5agg(1:nKernCoeffsDraws*2)

        else
            param_LL = val4
            opt_param_line = val5
        end if

        do i=1,nKernCoeffsDraws 
            opt_param_block(i,:) = opt_param_line( (i-1)*2 + 1: i*2 )
        end do 

        ! Select the hyperparameters with the biggest likelihood
        maxLL_ind = minloc(param_LL)
        temp_coeff_vec_soln = opt_param_block(maxLL_ind(1),:)
                  
        A_kernCoeff = 10.0**temp_coeff_vec_soln

        ! Update the gpr coefficients, which automatically updates the policy itself
        call createK(K_A, A_kernCoeff,1)
        A_gprCoeff = matmul(inv(K_A + sn_star**2.0*eyeN),tset_A)

        if(id .EQ. 0) then 

            open(3,file="A_kern_info.txt")
                write(3,*) param_LL
                write(3,*) maxLL_ind(1)
                write(3,*) temp_coeff_vec
                write(3,*) temp_coeff_vec_soln

            close(3)

            open(3,file="tset_VR.txt")

                do i=1,nInputs
                    write(3,*) VL + (VH-VL)/(1.0 + exp(-tset_VRz(i)))
                end do 

            close(3)

            open(3,file="tset_VD.txt")

                do i=1,nInputs
                    write(3,*) VL + (VH-VL)/(1.0 + exp(-tset_VDz(i)))
                end do 

            close(3)

            open(3,file="tset_A.txt")

                ! write(3,*) gridPointsX(2,:)*(1.0-llambda) + tset_A
                ! write(3,*) tset_A
                do i=1,nInputs
                    write(3,*) bL + (bH-bL)/(1.0 + exp(-tset_Az(i)))
                end do 

            close(3)

        end if 

        ! Now check convergence by evaluating value function again on same points 
        ! as drawn at beginning
        if( id .EQ. 0) then 
            do i=1,N_validate 
                v_new_points(i) = vRepay(convergence_points(:,i))
            end do 

            ! Sup-norm convergence
            ddist = maxval(abs(v_new_points-v_old_points))

            do i=1,nproc-1
                call MPI_SEND( ddist, 1, MPI_DOUBLE_PRECISION, i, &
                    i, MPI_COMM_WORLD, ierr )
            end do 
        
        else 

            call MPI_RECV(ddist, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr ) 

        end if 


        ! Have one processor print out the key results so they don't overwrite each other
        if(id .EQ. 0) then 

            call cpu_time(time2)

            open(3,file="iters.txt")

                write(3,*) "Iteration Number ", ccounter, ddist, time2-time1

            close(3)

            open(200,file="VD_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*) VD_gprCoeff(i)
                end do

            close(200)

            open(200,file="q_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*)  q_gprCoeff(i)
                end do

            close(200)

            open(200,file="qd_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*) qd_gprCoeff(i)
                end do

            close(200)

            open(200,file="qRN_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*)  qRN_gprCoeff(i)
                end do

            close(200)

            open(200,file="qdRN_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*) qdRN_gprCoeff(i)
                end do

            close(200)

            open(200,file="VR_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*) VR_gprCoeff(i)
                end do

            close(200)

            open(200,file="A_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*)  A_gprCoeff(i)
                end do

            close(200)

            open(200,file="I_gpr_coeffs.txt")
        
                do i=1,nInputs
                    write(200,*) I_gprCoeff(i)
                end do

            close(200)


            open(200,file="VD_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  VD_kernCoeff(i)
                end do

            close(200)

            open(200,file="q_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  q_kernCoeff(i)
                end do

            close(200)

            open(200,file="qd_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  qd_kernCoeff(i)
                end do

            close(200)

            open(200,file="qRN_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  qRN_kernCoeff(i)
                end do

            close(200)

            open(200,file="qdRN_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  qdRN_kernCoeff(i)
                end do

            close(200)

            open(200,file="VR_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  VR_kernCoeff(i)
                end do

            close(200)

            open(200,file="A_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  A_kernCoeff(i)
                end do

            close(200)

            open(200,file="I_kern_coeffs.txt")
        
                do i=1,2
                    write(200,*)  I_kernCoeff(i)
                end do

            close(200)

            open(200,file="means.txt") 
        
                write(200,*) VR_mean 
                write(200,*) VD_mean 
                write(200,*) A_mean 
                write(200,*) I_mean
                write(200,*) q_mean 
                write(200,*) qd_mean 
                write(200,*) qRN_mean 
                write(200,*) qdRN_mean

            close(200)

            open(200,file="stds.txt") 
        
                write(200,*) VR_std 
                write(200,*) VD_std
                write(200,*) A_std 
                write(200,*) I_std
                write(200,*) q_std
                write(200,*) qd_std 
                write(200,*) qRN_std 
                write(200,*) qdRN_std

            close(200)

        end if

    end do


    ! With the VFI having converged, we now simulate the model to get the relevant moments
    ! Note this needs only be done on one processor
    if(id .EQ. 0) then 

        call simulated_moments(moment_vec)

        open(21,file="sim_moments.txt") 

            write(21,*) "Iterations: ", ccounter
            write(21,*) "Value function Sup-norm: ", ddist
            write(21,*) "Average Debt/GDP: ", moment_vec(1)
            write(21,*) "Ann. Def. Freq: ", moment_vec(2)
            write(21,*) "Mean Spread: ", moment_vec(3)
            write(21,*) "Spread Vol: ", moment_vec(4)
            write(21,*) "Attn Frac: ", moment_vec(5)
            write(21,*) "Spread Cyclicality: ", moment_vec(6)
            write(21,*) "Trade Share Cyclicality: ", moment_vec(7)
            write(21,*) "RN Spread Share: ", moment_vec(8)
            write(21,*) " "
            write(21,*) "Discount Factor: ", bbeta 
            write(21,*) "Lender CRRA: ", CRRA_LEND
            write(21,*) "ythresh: ", ythresh
            write(21,*) "MP: ", MP
            write(21,*) "Implied y0: ", ylev
            write(21,*) "Implied y1: ", ycurv 
            write(21,*) "Info Cost: ", kkappa 
            write(21,*) "Maturity: ", llambda 
            write(21,*) "sigma-z: ", ssigma_z 
            write(21,*) "rho-z: ", rrho_z 
            write(21,*) "SovCRRA: ", ggamma 
            write(21,*) " "
            write(21,*) "Grid Points: ", nInputs 
            write(21,*) "Quadrature Size: ", quad_size 

        close(21)

    end if 

    if(mpi_on) then
        call mpi_finalize(ierr)
    end if

end program main
