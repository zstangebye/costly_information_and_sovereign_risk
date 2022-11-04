! This file contains auxiliary functions used to solve the
! sovereign debt model with machine learning methods

module aux_functions

    use global_vars

    implicit none

    contains 

    !!!! ECONOMIC/INTUITIVE FUNCTIONS

    ! The country's flow utility
    function sov_utility(c) result(sov_utility1)
        implicit none

        double precision :: c
        double precision :: sov_utility1 
        
        if(ggamma > 1.01 .OR. ggamma < .99) then 
            if(c > 0) then

                sov_utility1 = c**(1.0-ggamma)/(1.0-ggamma)

            else

                sov_utility1 = -1.0e10

            end if
        else 
            if(c > 0) then

                sov_utility1 = log(c)

            else

                sov_utility1 = -1.0e10

            end if
        end if 

    end function sov_utility

    ! Lenders' marginal utility
    function lender_marginal_u(c) result(lender_marginal_u1)
        implicit none 

        double precision, intent(in) :: c
        double precision :: lender_marginal_u1

        if(c <= 0.0) then
            lender_marginal_u1 = 1.0e8
        else
            lender_marginal_u1 = c**(-CRRA_LEND)
        end if

    end function lender_marginal_u

    ! Output in default, i.e., output - cost
    function ydef(y) result(ydef1)
        implicit none

        double precision, intent(in) :: y
        double precision :: ydef1

        ydef1 = y - max(ylev*y + ycurv*y**2.0,0.0)

    end function ydef


    ! Now we define value/pricing/policy functions
    function vRepay(xIn) result(vRepay1)! Repayment value (Y_t, B_t, l_t)
        implicit none 

        double precision, intent(in), dimension(nDims-1) :: xIn
        double precision, dimension(nDims-1) :: xTransformed
        double precision :: vRepay1

            xTransformed(1) = (xIn(1) - yL)/(yH - yL)
            xTransformed(2) = (xIn(2) - bL)/(bH - bL)


            vRepay1 = VL + (VH - VL)/(1.0 + & 
                exp(-VR_std*GPR_approx(xTransformed,VR_gprCoeff,VR_kernCoeff,scaledGridPointsX(1:2,:),2)-VR_mean))

    end function vRepay  

    function vDefault(xIn) result(vDefault1) ! Value of default (Y_t, B_{t},m_t)
        implicit none 

        double precision, intent(in), dimension(nDims) :: xIn
        double precision, dimension(nDims) :: xTransformed
        double precision :: vDefault1

            xTransformed(1) = (xIn(1) - yL)/(yH - yL)
            xTransformed(2) = (xIn(2) - bL)/(bH - bL)
            xTransformed(3) = (xIn(3) - mL)/(mH - mL)

            vDefault1 = VL + (VH - VL)/(1.0 + & 
                exp(-VD_std*GPR_approx(xTransformed,VD_gprCoeff,VD_kernCoeff,scaledGridPointsM,3)-VD_mean))


    end function vDefault

    ! Borrowing policy (Y_t, B_t)
    ! We will require this to be bounded in [bL,bH] by using a logit-structure
    function APol(xIn) result(APol1) ! Borrowing policy (Y_t, B_t)
        implicit none 

        double precision, intent(in), dimension(nDims-1) :: xIn
        double precision, dimension(nDims-1) :: xTransformed
        double precision :: APol1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)

        APol1 = bL + (bH - bL)/(1.0 + & 
            exp(-A_std*GPR_approx(xTransformed,A_gprCoeff,A_kernCoeff,scaledGridPointsX(1:2,:),2)-A_mean))

    end function APol 

    ! Information acquisition policy (Y_t, B_{t+1})
    ! We will require this to be bounded in [0,1] by using a logit-structure
    function infoPol(xIn) result(infoPol1)! Information acquisition policy (Y_t, B_{t+1})
        implicit none 

        double precision, intent(in), dimension(nDims-1) :: xIn
        double precision, dimension(nDims-1) :: xTransformed
        double precision :: infoPol1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)

        infoPol1 = 1.0/(1.0 + exp(-I_std*GPR_approx(xTransformed,I_gprCoeff,I_kernCoeff,scaledGridPointsX(1:2,:),2)-I_mean) )

    end function infoPol 

    ! This is a binary default function
    function defFun(todayY,todayB,todayM) result(defFun1) ! Default indicator (Y_t, B_t, m_t)
        implicit none 

        double precision, intent(in) :: todayY, todayB, todayM
        double precision, dimension(nDims) :: xIn 
        double precision :: defFun1

        xIn(1) = todayY
        xIn(2) = todayB 
        xIn(3) = todayM 

        if(vDefault(xIn) > vRepay(xIn(1:2))) then 
            defFun1 = 1.0
        else
            defFun1 = 0.0 
        end if 

    end function defFun

    ! Equilibrium repayment/issuance price (Y_t, B_{t+1}, x_t)
    ! Required to be bounded in [0,1] by using a logit-structure 
    ! This is valid because coupons are set to risk-free rate -> risk-free price = 1
    function defProb(xIn) result(def1) ! Equilibrium repayment/issuance price (Y_t, B_{t+1}, x_t)
        implicit none 

        double precision, intent(in), dimension(nDims) :: xIn
        double precision, dimension(nDims) :: xTransformed
        double precision :: def1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)
        xTransformed(3) = (xIn(3) - xL)/(xH - xL)

        def1 = 1.0/(1.0 + exp(-def_std*GPR_approx(xTransformed,def_gprCoeff,def_kernCoeff,scaledGridPointsX,3)-def_mean))

    end function defProb 

    ! Equilibrium repayment/issuance price (Y_t, B_{t+1}, x_t)
    ! Required to be bounded in [0,1] by using a logit-structure 
    ! This is valid because coupons are set to risk-free rate -> risk-free price = 1
    function q(xIn) result(q1) ! Equilibrium repayment/issuance price (Y_t, B_{t+1}, x_t)
        implicit none 

        double precision, intent(in), dimension(nDims) :: xIn
        double precision, dimension(nDims) :: xTransformed
        double precision :: q1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)
        xTransformed(3) = (xIn(3) - xL)/(xH - xL)

        q1 = 1.0/(1.0 + exp(-q_std*GPR_approx(xTransformed,q_gprCoeff,q_kernCoeff,scaledGridPointsX,3)-q_mean))

    end function q 

    ! Equilibrium repayment/issuance price (Y_t, B_{t+1}, x_t)
    ! Required to be bounded in [0,1] by using a logit-structure 
    function qd(xIn) result(qd1) ! Equilibrium default price (Y_t, B_t, m_t)
        implicit none 

        double precision, intent(in), dimension(nDims) :: xIn
        double precision, dimension(nDims) :: xTransformed
        double precision :: qd1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)
        xTransformed(3) = (xIn(3) - mL)/(mH - mL)

        qd1 = 1.0/(1.0 + exp(-qd_std*GPR_approx(xTransformed,qd_gprCoeff,qd_kernCoeff,scaledGridPointsM,3)-qd_mean))

    end function qd

    ! Off-equilibrium RN repayment/issuance price (Y_t, B_{t+1}, x_t)
    ! Required to be bounded in [0,1] by using a logit-structure 
    ! This is valid because coupons are set to risk-free rate -> risk-free price = 1
    function qRN(xIn) result(q1) ! Equilibrium repayment/issuance RN price (Y_t, B_{t+1}, x_t)
        implicit none 

        double precision, intent(in), dimension(nDims) :: xIn
        double precision, dimension(nDims) :: xTransformed
        double precision :: q1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)
        xTransformed(3) = (xIn(3) - xL)/(xH - xL)

        q1 = 1.0/(1.0 + exp(-qRN_std*GPR_approx(xTransformed,qRN_gprCoeff,qRN_kernCoeff,scaledGridPointsX,3)-qRN_mean))

    end function qRN  

    ! Off-equilibrium RN repayment/issuance price (Y_t, B_{t+1}, x_t)
    ! Required to be bounded in [0,1] by using a logit-structure 
    function qdRN(xIn) result(qd1) ! Equilibrium default RN price (Y_t, B_t, m_t)
        implicit none 

        double precision, intent(in), dimension(nDims) :: xIn
        double precision, dimension(nDims) :: xTransformed
        double precision :: qd1

        xTransformed(1) = (xIn(1) - yL)/(yH - yL)
        xTransformed(2) = (xIn(2) - bL)/(bH - bL)
        xTransformed(3) = (xIn(3) - mL)/(mH - mL)

        qd1 = 1.0/(1.0 + exp(-qdRN_std*GPR_approx(xTransformed,qdRN_gprCoeff,qdRN_kernCoeff,scaledGridPointsM,3)-qdRN_mean))

    end function qdRN  

    !! TECHNICAL/SOLUTION-BASED FUNCTIONS

    ! The general, 1-D optimization routine using golden-section search, which is gradient-free and more reliable. 
    ! It applies to both the forecaster's problem and the sovereign Bellman
    subroutine localMaximizeFun(xLL, xHH, xSoln, optVal, curr_state, whichProb)
        implicit none

        double precision, intent(in) :: xLL, xHH ! The bounds over which we're maximizing
        double precision, intent(inout) :: xSoln, optVal ! The results of the optimization

        integer, intent(in) :: whichProb ! Specifices which problem we're maximizing
        ! whichProb = 1 -> Forecaster's problem
        ! whichProb = 2 -> Sovereign's Bellman problem
        double precision, intent(in), dimension(nDims) :: curr_state

        double precision, parameter :: aalpha1 = (3.0-sqrt(5.0))/2.0
        double precision, parameter :: aalpha2 = (sqrt(5.0)-1.0)/2.0
        double precision :: x1, x2, v1, v2, d ! Used in golden-section algorithm
        integer :: gs_counter ! Keeps track of iterations

        x1 = xLL + aalpha1*(xHH-xLL)
        x2 = xLL + aalpha2*(xHH-xLL)

        if(whichProb .EQ. 1) then 
            v1 = -defMSE(x1,curr_state(1),curr_state(2))
            v2 = -defMSE(x2,curr_state(1),curr_state(2))
        else 
            v1 = repayBellman(x1,curr_state(1),curr_state(2))
            v2 = repayBellman(x2,curr_state(1),curr_state(2))
        end if 

        d = aalpha1*aalpha2*(xHH-xLL)
        gs_counter = 2

        ! Here's the loop
        do while(opt_tol < d)
            gs_counter = gs_counter + 1
            d = aalpha2*d
            if(v2 > v1) then

                x1 = x2
                x2 = x2+d
                v1 = v2
                if(whichProb .EQ. 1) then 
                    v2 = -defMSE(x2,curr_state(1),curr_state(2))
                 else 
                    v2 = repayBellman(x2,curr_state(1),curr_state(2))
                end if 
                
            else

                x2 = x1
                x1 = x1-d
                v2 = v1
                if(whichProb .EQ. 1) then 
                    v1 = -defMSE(x1,curr_state(1),curr_state(2))
                 else 
                    v1 = repayBellman(x1,curr_state(1),curr_state(2))
                end if 
                
            end if
        end do 

        ! Post-convergence, go with whatever interior option is better
        if(v2 > v1) then
            xSoln = x2
            optVal = v2
        else
            xSoln = x1
            optVal = v1
        end if

    end subroutine localMaximizeFun

    ! This function delivers the sovereign Bellman in repayment given a state (Y,B) and a 
    ! possible debt issuance level B'. It is maximized over B' in practice
    function repayBellman(issued_b,current_y,current_b) result(repayBellman1)
        implicit none 

        double precision, intent(in) :: issued_b,current_y,current_b
        double precision :: zL1, zH1, ytempObj, mtempObj, curr_u, curr_c, valTomorrow, RvalTomorrow
        double precision :: repayBellman1 

        double precision, dimension(quad_size) :: z_gc_points, yvec, mvec, xvec, mvec2
        double precision, dimension(nDims-1) :: tom_repay_state
        double precision, dimension(nDims) :: tom_default_state, curr_q_state, curr_def_state

        integer :: iy, im, ix ! Loop counter 

        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        tom_repay_state(2) = issued_b 
        tom_default_state(2) = issued_b 

        curr_q_state(1) = current_y 
        curr_q_state(2) = issued_b 
        

        ! mtom and xtod can be treated as independent for purpose of unconditional expectation since they are additively 
        ! separable. This simplifies computation substantially and is valid even if they are not independent

        ! Start with the continuation value
        do iy=1,quad_size 
            tom_repay_state(1) = exp(z_gc_points(iy))
            tom_default_state(1) = tom_repay_state(1)

            RvalTomorrow = vRepay(tom_repay_state)

            do im=1,quad_size
                tom_default_state(3) = m_GCnodes(im)

                mtempObj = max(RvalTomorrow,vDefault(tom_default_state))

                mvec(im) = GCweights(im)*mtempObj*sqrt(1.0-GCnodes(im)**2.0)*logit_pdf_est(m_GCnodes(im))

                mvec2(im) = vDefault(tom_default_state)

            end do 
            ytempObj = (mH-mL)/2.0*sum(mvec)

            yvec(iy) = GCweights(iy)*ytempObj*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust
            
        end do 
        valTomorrow = (zH1-zL1)/2.0*sum(yvec)

        ! Now do the flow utility; first check if this is a significant dilution state 
        ! and if so rule it out as in HMS(2016)

        curr_q_state(3) = x_GCnodes(1) ! Assume the worst possible current price 
        
        ! Need to check issuance in zero-debt state separately; if price is not high enough, 
        ! then set policy to zero issuance

        if( (q(curr_q_state) < min_q) .AND. (current_b < opt_tol ) ) then 

            curr_u = sov_utility(current_y)
            
        elseif( (q(curr_q_state) < min_q) .AND. (issued_b > (1.0-llambda)*current_b ) ) then 
        
            ! Ensure that this could never be optimal, i.e., impose constraint
            curr_u = -5.0e10 

        else 

            if( infoPol(curr_q_state(1:2)) > .001) then 

                do ix=1,quad_size
                    curr_q_state(3) = x_GCnodes(ix)
                    curr_c = current_y - (llambda + coup)*current_b + q(curr_q_state)*(issued_b - (1.0-llambda)*current_b)

                    xvec(ix) = GCweights(ix)*sov_utility(curr_c)*sqrt(1.0-GCnodes(ix)**2.0)* &
                        normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 
                
                end do 
                curr_u = (xH-xL)/2.0*sum(xvec)

            else 

                curr_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)
                curr_c = current_y - (llambda + coup)*current_b + q(curr_q_state)*(issued_b - (1.0-llambda)*current_b)

                curr_u = sov_utility(curr_c)


            end if 
        
        end if 

        repayBellman1 = curr_u + bbeta*valTomorrow

    end function repayBellman

    ! This function computes the mean-squared-error plus the information cost for the forecaster's problem. It is optimized 
    ! over current_rho in the solution given a state (Y,B')
    function defMSE(current_rho,current_y,issued_b) result(defMSE1)
        implicit none 

        double precision, intent(in) :: current_rho, current_y, issued_b ! current state
        double precision :: zL1, zH1, ytempObj, mtempObj, xtempObj, chosen_rho_tom, & 
            mean_norm_m, vol_norm_m, defInd, tom_qd, marg_return
        double precision :: defMSE1, return_bar, ex_post_return_dev, MSE, MI

        integer :: ix, ixtom, im, iy, il ! Loop counters

        double precision, dimension(quad_size) :: z_gc_points, mvec, yvec, xvec, xvec_tom
        double precision, dimension(nDims-1) ::  tom_info_state, tom_borr_state
        double precision, dimension(nDims) ::  tom_q_state, tom_qd_state



        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        ! If the country repays even in the worst state, then the recovery rate cannot possibly
        ! affect the returns and so we skip the expensive maximization here in place of an objective minimized 
        ! at the lower bound

        if( defFun(exp(z_gc_points(1)),issued_b,m_GCnodes(1)) < 0.5 ) then  

            defMSE1 = 5.0 + current_rho 

        else 

            MI = -.5*log(1-current_rho**2.0)

            tom_borr_state(2) = issued_b
            tom_qd_state(2) = issued_b ! Tomorrow, if the country defaults it is nominally on the level it issued today

            do iy=1,quad_size 
            ! This expectation is written to minimize the number of equilibrium function calls, since they are GPR calls and 
            ! can be more expensive than storing scalars
                tom_borr_state(1) = exp(z_gc_points(iy)) 
                tom_info_state(1) = exp(z_gc_points(iy))
                tom_q_state(1) = exp(z_gc_points(iy))
                tom_qd_state(1) = exp(z_gc_points(iy))

                tom_info_state(2) = APol(tom_borr_state)
                chosen_rho_tom = infoPol(tom_info_state)
                tom_q_state(2) = tom_info_state(2)

                do im=1,quad_size 
                    tom_qd_state(3) = m_GCnodes(im)
                    tom_qd = qd(tom_qd_state)

                        defInd = defFun(exp(z_gc_points(iy)),issued_b,m_GCnodes(im))

                        if( defInd > 0.5) then 

                            mtempObj = tom_qd 

                        else 

                            if(chosen_rho_tom > 0.01) then 
                            ! If this is not the case, we can save a dimension since we know it will be constant in x

                                do ix=1,quad_size ! Here, it will not be constant in x, so we take an expectation over x

                                    tom_q_state(3) = x_GCnodes(ix)

                                    marg_return = llambda + coup + (1.0-llambda)*q(tom_q_state) 

                                    xvec(ix) = GCweights(ix)*marg_return*sqrt(1.0-GCnodes(ix)**2.0)* &
                                        normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 
                                end do 
                                mtempObj = (xH-xL)/2.0*sum(xvec)

                            else ! Here, it is constant in x, so we evaluate it at roughly the midpoint of x;
                            ! If quad_size is odd, it will be exactly the midpoint
                            
                                tom_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)

                                mtempObj = llambda + coup + (1.0-llambda)*q(tom_q_state) 

                            end if 

                        end if 

                    mvec(im) = GCweights(im)*mtempObj*sqrt(1.0-GCnodes(im)**2.0)* & 
                        logit_pdf_est(m_GCnodes(im))

                end do 
                ytempObj = (mH-mL)/2.0*sum(mvec)

                yvec(iy) = GCweights(iy)*ytempObj*sqrt(1.0-GCnodes(iy)**2.0)* & 
                    normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

            end do 

            return_bar = (zH1-zL1)/2.0*sum(yvec)

            ! Now, compute the variance
            do iy=1,quad_size 
            ! This expectation is written to minimize the number of equilibrium function calls, since they are GPR calls and 
            ! can be more expensive than storing scalars
                tom_borr_state(1) = exp(z_gc_points(iy)) 
                tom_info_state(1) = exp(z_gc_points(iy))
                tom_q_state(1) = exp(z_gc_points(iy))
                tom_qd_state(1) = exp(z_gc_points(iy))

                tom_info_state(2) = APol(tom_borr_state)
                chosen_rho_tom = infoPol(tom_info_state)
                tom_q_state(2) = tom_info_state(2)

                do ix = 1,quad_size

                    mean_norm_m = norm_mu_est + current_rho*(x_GCnodes(ix) - norm_mu_est) ! Mean and volatility of x and m are the same
                    vol_norm_m = sqrt(1.0-current_rho**2.0)*norm_sigma_est

                    do im=1,quad_size 
                        tom_qd_state(3) = m_GCnodes(im)
                        tom_qd = qd(tom_qd_state)


                            defInd = defFun(exp(z_gc_points(iy)),issued_b,m_GCnodes(im))

                            if( defInd > 0.5) then 

                                mtempObj = (tom_qd - return_bar)**2.0 

                            else 

                                if(chosen_rho_tom > 0.01) then 
                                ! If this is not the case, we can save a dimension since we know it will be constant in x

                                    do ixtom=1,quad_size ! Here, it will not be constant in x, so we take an expectation over x

                                        tom_q_state(3) = x_GCnodes(ixtom)

                                        ex_post_return_dev = (llambda + coup + (1.0-llambda)*q(tom_q_state) - return_bar)**2.0

                                        xvec_tom(ixtom) = GCweights(ixtom)*ex_post_return_dev*sqrt(1.0-GCnodes(ixtom)**2.0)* &
                                            normPDF( (x_GCnodes(ixtom) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 
                                    end do 
                                    mtempObj = (xH-xL)/2.0*sum(xvec_tom)

                                else ! Here, it is constant in x, so we evaluate it at roughly the midpoint of x;
                                ! If quad_size is odd, it will be exactly the midpoint
                                
                                    tom_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)

                                    mtempObj = (llambda + coup + (1.0-llambda)*q(tom_q_state) - return_bar)**2.0

                                end if 

                            end if 

                        mvec(im) = GCweights(im)*mtempObj*sqrt(1.0-GCnodes(im)**2.0)* & 
                            logit_pdf(m_GCnodes(im),vol_norm_m,mean_norm_m)

                    end do 
                    xtempObj = (mH-mL)/2.0*sum(mvec)

                    xvec(ix) = GCweights(ix)*xtempObj*sqrt(1.0-GCnodes(ix)**2.0)* &
                        normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                end do 

                ytempObj = (xH-xL)/2.0*sum(xvec)

                yvec(iy) = GCweights(iy)*ytempObj*sqrt(1.0-GCnodes(iy)**2.0)* & 
                    normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

            end do
            MSE = (zH1-zL1)/2.0*sum(yvec)

            defMSE1 = sqrt(MSE) + kkappa*MI

        end if 

    end function defMSE 

    ! This function delivers the Bellman equation in default. It requires no maximization
    function defaultBellman(current_y,current_b,current_m) result(defaultBellman1)
        implicit none 

        double precision, intent(in) :: current_y,current_b, current_m
        double precision :: zL1, zH1, ytempObj, curr_u, valTomorrow
        double precision :: defaultBellman1 

        double precision, dimension(quad_size) :: z_gc_points, yvec
        double precision, dimension(nDims-1) :: tom_repay_stateR
        double precision, dimension(nDims) :: tom_default_stateD, tom_default_stateR

        integer :: iy, im, ix ! Loop counter 

        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z
 
        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        tom_repay_stateR(2) = current_b*current_m 
        tom_default_stateR(2) = current_b*current_m 
        tom_default_stateD(2) = current_b 

        tom_default_stateR(3) = current_m
        tom_default_stateD(3) = current_m


        ! Start with the continuation value
        do iy=1,quad_size 
            tom_repay_stateR(1) = exp(z_gc_points(iy))
            tom_default_stateR(1) = tom_repay_stateR(1)
            tom_default_stateD(1) = tom_repay_stateR(1)

            ytempObj = pphi*max(vRepay(tom_repay_stateR),vDefault(tom_default_stateR)) + (1.0-pphi)*vDefault(tom_default_stateD)

            yvec(iy) = GCweights(iy)*ytempObj*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 
        end do 
        valTomorrow = (zH1-zL1)/2.0*sum(yvec)

        ! Now do the flow utility
        if(ydef(current_y) > 0.0) then 
            curr_u = sov_utility(ydef(current_y))
        else  
            curr_u = -1.0e10 
        end if

        defaultBellman1 = curr_u + bbeta*valTomorrow

    end function defaultBellman

    ! Use interval bisection to find the market clearing q for a given debt issuance
    function new_q(current_y,issued_b,current_x) result(new_q1) 
        implicit none 

        double precision, intent(in) :: current_y, issued_b, current_x ! current state
        double precision :: qL, qH, q_guess, curr_gap, new_q1 

        qL = 0.0 
        qH = 1.0 
        curr_gap = 1.0

        do while( qH - qL > opt_tol )

            q_guess = .5*(qL + qH)
            curr_gap = q_gap( q_guess, current_y,issued_b,current_x)

            ! Note the q-gap function is a decreasing one and we're searching 
            ! for its root
            if(curr_gap > 0.0) then
                qL = q_guess 
            else
                qH = q_guess 
            end if 

            ! Add a check if we stumble across zero early
            if(abs(curr_gap) .LE. opt_tol) then 
                qL = q_guess 
                qH = q_guess 
            end if 

        end do 

        new_q1 = .5*(qL + qH)

    end function new_q

    ! This function computes the gap between a conjectured q and the one the lenders will take 
    ! given the equilibrium risk profile
    function q_gap(curr_q,current_y,issued_b,current_x) result(q_gap1)
        implicit none

        double precision, intent(in) :: curr_q, current_y, issued_b, current_x ! current state
        double precision :: zL1, zH1, ytempObj_top, ytempObj_bottom, mtempObj_top, mtempObj_bottom, chosen_rho_tom, chosen_rho, & 
            mean_norm_m, vol_norm_m, MU_num, MU_denom, lmu_top, lmu_bottom, defInd, tom_qd, marg_return
        double precision :: q_gap1

        integer :: ix, im, iy ! Loop counters

        double precision, dimension(quad_size) :: z_gc_points, mvec_top, mvec_bottom, yvec_top, yvec_bottom, xvec_top, & 
            xvec_bottom
        double precision, dimension(nDims-1) :: curr_info_state, tom_borr_state, tom_info_state
        double precision, dimension(nDims) ::  tom_q_state, tom_qd_state 

        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        curr_info_state(1) = current_y 
        curr_info_state(2) = issued_b

        ! Check what the information policy is for tomorrow; update expectations accordingly
        chosen_rho = infoPol(curr_info_state)

        mean_norm_m = norm_mu_est + chosen_rho*(current_x - norm_mu_est) ! Mean and volatility of x and m are the same
        vol_norm_m = sqrt(1.0-chosen_rho**2.0)*norm_sigma_est

        do iy=1,quad_size 
        ! This expectation is written to minimize the number of equilibrium function calls, since they are GPR calls and 
        ! can be more expensive than storing scalars
            tom_borr_state(1) = exp(z_gc_points(iy))
            tom_borr_state(2) = issued_b 

            tom_info_state(1) = exp(z_gc_points(iy))
            tom_info_state(2) = APol(tom_borr_state)

            chosen_rho_tom = infoPol(tom_info_state)

            tom_q_state(1) = exp(z_gc_points(iy))
            tom_q_state(2) = tom_info_state(2)

            tom_qd_state(1) = exp(z_gc_points(iy))
            tom_qd_state(2) = issued_b ! Tomorrow, if the country defaults it is nominally on the level it issued today

            do im=1,quad_size 
                defInd = defFun(exp(z_gc_points(iy)),issued_b,m_GCnodes(im))

                tom_qd_state(3) = m_GCnodes(im)

                tom_qd = qd(tom_qd_state)

                if(chosen_rho_tom > 0.01) then 
                ! If this is not the case, we can save a dimension since we know it will be constant in x

                    do ix=1,quad_size ! Here, it will not be constant in x, so we take an expectation over x

                        tom_q_state(3) = x_GCnodes(ix)

                        marg_return = (1.0-defInd)*(llambda + coup + (1.0-llambda)*q(tom_q_state) ) + defInd*tom_qd 

                        lmu_bottom = lender_marginal_u( (w-issued_b*curr_q)*(1+r) + issued_b*marg_return )

                        lmu_top = lmu_bottom*marg_return

                        xvec_top(ix) = GCweights(ix)*lmu_top*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                        xvec_bottom(ix) = GCweights(ix)*lmu_bottom*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                    end do 
                    mtempObj_top = (xH-xL)/2.0*sum(xvec_top)
                    mtempObj_bottom = (xH-xL)/2.0*sum(xvec_bottom)

                else ! Here, it is constant in x, so we evaluate it at roughly the midpoint of x;
                ! If quad_size is odd, it will be exactly the midpoint
                
                    tom_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)

                    marg_return = (1.0-defInd)*(llambda + coup + (1.0-llambda)*q(tom_q_state) ) + defInd*tom_qd 

                    mtempObj_bottom = lender_marginal_u( (w-issued_b*curr_q)*(1+r) + issued_b*marg_return )
                    mtempObj_top = mtempObj_bottom*marg_return

                end if 

                mvec_top(im) = GCweights(im)*mtempObj_top*sqrt(1.0-GCnodes(im)**2.0)* & 
                    logit_pdf(m_GCnodes(im),vol_norm_m,mean_norm_m)
                mvec_bottom(im) = GCweights(im)*mtempObj_bottom*sqrt(1.0-GCnodes(im)**2.0)* & 
                    logit_pdf(m_GCnodes(im),vol_norm_m,mean_norm_m)

            end do 
            ytempObj_top = (mH-mL)/2.0*sum(mvec_top)
            ytempObj_bottom = (mH-mL)/2.0*sum(mvec_bottom)

            yvec_top(iy) = GCweights(iy)*ytempObj_top*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 
            yvec_bottom(iy) = GCweights(iy)*ytempObj_bottom*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

        end do 
        MU_num = (zH1-zL1)/2.0*sum(yvec_top)
        MU_denom = (zH1-zL1)/2.0*sum(yvec_bottom)

        q_gap1 = MU_num/((1.0+r)*MU_denom) - curr_q 
 
    end function q_gap

    ! Use interval bisection to find the market clearing q for a given debt issuance
    function new_qd(current_y,issued_b,current_m) result(new_qd1) 
        implicit none 

        double precision, intent(in) :: current_y, issued_b, current_m ! current state
        double precision :: qL, qH, q_guess, curr_gap, new_qd1 

        qL = 0.0 
        qH = 1.0 
        curr_gap = 1.0

        do while( qH - qL > opt_tol )

            q_guess = .5*(qL + qH)
            curr_gap = qd_gap( q_guess, current_y,issued_b,current_m)

            ! Note the q-gap function is a decreasing one and we're searching 
            ! for its root
            if(curr_gap > 0.0) then
                qL = q_guess 
            else
                qH = q_guess 
            end if 

            ! Add a check if we stumble across zero early
            if(abs(curr_gap) .LE. opt_tol) then 
                qL = q_guess 
                qH = q_guess 
            end if 

        end do 

        new_qd1 = .5*(qL + qH)

    end function new_qd

    ! This updates the pricing function in default
    function qd_gap(curr_qd,current_y,current_b,current_m) result(qd_gap1)
        implicit none 

        double precision, intent(in) :: current_y, current_b, current_m, curr_qd ! Current state 
        double precision :: zL1, zH1, ytempObj_top, ytempObj_bottom, chosen_rho_tom, MU_num, & 
            MU_denom, lmu_top, lmu_bottom, defInd, tom_q, tom_qdR, tom_qdD, marg_return_R, marg_return_D
        double precision :: qd_gap1

        integer :: ix, iy ! Loop counters

        double precision, dimension(quad_size) :: z_gc_points, yvec_top, yvec_bottom, xvec_top, xvec_bottom
        double precision, dimension(nDims-1) :: future_info_state, future_borrow_state
        double precision, dimension(nDims) :: tom_qd_state, tom_q_state
        
        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        do iy=1,quad_size 

            future_borrow_state(1) = exp(z_gc_points(iy)) 
            future_borrow_state(2) = current_m*current_b

            future_info_state(1) = future_borrow_state(1)
            future_info_state(2) = APol(future_borrow_state)

            chosen_rho_tom = infoPol(future_info_state)

            tom_q_state(1) = future_borrow_state(1)
            tom_q_state(2) = future_info_state(2)

            defInd = defFun(exp(z_gc_points(iy)),current_b*current_m,current_m)

            tom_qd_state(1) = exp(z_gc_points(iy))
            tom_qd_state(2) = current_b*current_m
            tom_qd_state(3) = current_m

            tom_qdR = qd(tom_qd_state)

            tom_qd_state(2) = current_b

            tom_qdD = qd(tom_qd_state)


            if(chosen_rho_tom > 0.01) then 

                do ix=1,quad_size 

                    tom_q_state(3) = x_GCnodes(ix)

                    marg_return_R = current_m*((1.0-defInd)*(llambda + coup + (1.0-llambda)*q(tom_q_state) ) + defInd*tom_qdR)
                    marg_return_D = tom_qdD

                    lmu_bottom = pphi*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_R) + & 
                        (1.0-pphi)*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_D)
                    
                    lmu_top = pphi*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_R)* &
                        marg_return_R + (1.0-pphi)*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_D)* &
                        marg_return_D

                    xvec_top(ix) = GCweights(ix)*lmu_top*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                    xvec_bottom(ix) = GCweights(ix)*lmu_bottom*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                end do 

                ytempObj_top = (xH-xL)/2.0*sum(xvec_top)
                ytempObj_bottom = (xH-xL)/2.0*sum(xvec_bottom)

            else 

                tom_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)

                marg_return_R = current_m*((1.0-defInd)*(llambda + coup + (1.0-llambda)*q(tom_q_state) ) + defInd*tom_qdR )               
                marg_return_D = tom_qdD

                ytempObj_bottom = (pphi*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_R) + & 
                    (1.0-pphi)*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_D))
                    
                ytempObj_top = (pphi*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_R)* &
                    marg_return_R + (1.0-pphi)*lender_marginal_u( (w-curr_qd*current_b)*(1.0+r) + current_b*marg_return_D)*marg_return_D)

            end if

            yvec_top(iy) = GCweights(iy)*ytempObj_top*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 
            yvec_bottom(iy) = GCweights(iy)*ytempObj_bottom*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

        end do 

        MU_num = (zH1-zL1)/2.0*sum(yvec_top)
        MU_denom = (zH1-zL1)/2.0*sum(yvec_bottom)


        qd_gap1 = MU_num/((1.0+r)*MU_denom) - curr_qd

    end function qd_gap

    ! This function computes the expected default probability tomorrow given the current state today, i.e.,
    ! output level, issuance, and signal. It will give a sense of how `surprising' default tends to be but is not 
    ! used in the solution of the model
    function defProbTomorrow(current_y,issued_b,current_x) result(new_def1)
        implicit none

        double precision, intent(in) :: current_y, issued_b, current_x ! current state
        double precision :: zL1, zH1, ytempObj, mtempObj, chosen_rho_tom, chosen_rho, & 
            mean_norm_m, vol_norm_m
        double precision :: new_def1

        integer :: ix, im, iy ! Loop counters

        double precision, dimension(quad_size) :: z_gc_points, mvec, yvec
        double precision, dimension(nDims-1) :: curr_info_state

        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        curr_info_state(1) = current_y 
        curr_info_state(2) = issued_b 

        ! Check what the information policy is for tomorrow; update expectations accordingly
        chosen_rho = infoPol(curr_info_state)

        mean_norm_m = norm_mu_est + chosen_rho*(current_x - norm_mu_est) ! Mean and volatility of x and m are the same
        vol_norm_m = sqrt(1.0-chosen_rho**2.0)*norm_sigma_est

        do iy=1,quad_size 
            do im=1,quad_size 

                mtempObj = defFun(exp(z_gc_points(iy)),issued_b,m_GCnodes(im))

                mvec(im) = GCweights(im)*mtempObj*sqrt(1.0-GCnodes(im)**2.0)* & 
                    logit_pdf(m_GCnodes(im),vol_norm_m,mean_norm_m)

            end do 
            ytempObj = (mH-mL)/2.0*sum(mvec)

            yvec(iy) = GCweights(iy)*ytempObj*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

        end do 
        new_def1 = (zH1-zL1)/2.0*sum(yvec)

    end function defProbTomorrow

    ! This function updates the pricing function in repayment as priced by a risk-neutral agent 
    ! It does not influence equilibrium outcomes, but allows us to measure the equilibrium risk premium
    ! Notice that we can define RN prices without reference to the current price, so 
    ! no rootfinding is required
    function new_qRN(current_y,issued_b,current_x) result(new_q1)
        implicit none

        double precision, intent(in) :: current_y, issued_b, current_x ! current state
        double precision :: zL1, zH1, ytempObj_top, ytempObj_bottom, mtempObj_top, mtempObj_bottom, chosen_rho_tom, chosen_rho, & 
            mean_norm_m, vol_norm_m, MU_num, MU_denom, lmu_top, lmu_bottom, defInd, tom_qd, marg_return
        double precision :: new_q1

        integer :: ix, im, iy ! Loop counters

        double precision, dimension(quad_size) :: z_gc_points, mvec_top, mvec_bottom, yvec_top, yvec_bottom, xvec_top, & 
            xvec_bottom
        double precision, dimension(nDims-1) :: curr_info_state, tom_borr_state, tom_info_state
        double precision, dimension(nDims) :: tom_q_state, tom_qd_state

        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)

        curr_info_state(1) = current_y 
        curr_info_state(2) = issued_b

        ! Check what the information policy is for tomorrow; update expectations accordingly
        chosen_rho = infoPol(curr_info_state)

        mean_norm_m = norm_mu_est + chosen_rho*(current_x - norm_mu_est) ! Mean and volatility of x and m are the same
        vol_norm_m = sqrt(1.0-chosen_rho**2.0)*norm_sigma_est

        do iy=1,quad_size 
        ! This expectation is written to minimize the number of equilibrium function calls, since they are GPR calls and 
        ! can be more expensive than storing scalars
            tom_borr_state(1) = exp(z_gc_points(iy))
            tom_borr_state(2) = issued_b 

            tom_info_state(1) = exp(z_gc_points(iy))
            tom_info_state(2) = APol(tom_borr_state)

            chosen_rho_tom = infoPol(tom_info_state)

            tom_q_state(1) = exp(z_gc_points(iy))
            tom_q_state(2) = tom_info_state(2)

            tom_qd_state(1) = exp(z_gc_points(iy))
            tom_qd_state(2) = issued_b ! Tomorrow, if the country defaults it is nominally on the level it issued today

            do im=1,quad_size 
                defInd = defFun(exp(z_gc_points(iy)),issued_b,m_GCnodes(im))

                tom_qd_state(3) = m_GCnodes(im)

                tom_qd = qdRN(tom_qd_state)

                if(chosen_rho_tom > 0.01) then 
                ! If this is not the case, we can save a dimension since we know it will be constant in x

                    do ix=1,quad_size ! Here, it will not be constant in x, so we take an expectation over x

                        tom_q_state(3) = x_GCnodes(ix)

                        marg_return = (1.0-defInd)*(llambda + coup + (1.0-llambda)*qRN(tom_q_state) ) + defInd*tom_qd 

                        lmu_bottom = 1.0

                        lmu_top = lmu_bottom*marg_return

                        xvec_top(ix) = GCweights(ix)*lmu_top*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                        xvec_bottom(ix) = GCweights(ix)*lmu_bottom*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                    end do 
                    mtempObj_top = (xH-xL)/2.0*sum(xvec_top)
                    mtempObj_bottom = (xH-xL)/2.0*sum(xvec_bottom)

                else ! Here, it is constant in x, so we evaluate it at roughly the midpoint of x;
                ! If quad_size is odd, it will be exactly the midpoint
                
                    tom_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)

                    marg_return = (1.0-defInd)*(llambda + coup + (1.0-llambda)*qRN(tom_q_state) ) + defInd*tom_qd 

                    mtempObj_bottom = 1.0
                    mtempObj_top = mtempObj_bottom*marg_return

                end if 

                mvec_top(im) = GCweights(im)*mtempObj_top*sqrt(1.0-GCnodes(im)**2.0)* & 
                    logit_pdf(m_GCnodes(im),vol_norm_m,mean_norm_m)
                mvec_bottom(im) = GCweights(im)*mtempObj_bottom*sqrt(1.0-GCnodes(im)**2.0)* & 
                    logit_pdf(m_GCnodes(im),vol_norm_m,mean_norm_m)

            end do 
            ytempObj_top = (mH-mL)/2.0*sum(mvec_top)
            ytempObj_bottom = (mH-mL)/2.0*sum(mvec_bottom)

            yvec_top(iy) = GCweights(iy)*ytempObj_top*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 
            yvec_bottom(iy) = GCweights(iy)*ytempObj_bottom*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

        end do 
        MU_num = (zH1-zL1)/2.0*sum(yvec_top)
        MU_denom = (zH1-zL1)/2.0*sum(yvec_bottom)

        new_q1 = MU_num/((1.0+r)*MU_denom)

    end function new_qRN

    ! This updates the pricing function in default for the risk-neutral counterfactual lender
    ! It is not used in the equilibrium but only to compute risk-premia after the fact
    function new_qdRN(current_y,current_b,current_m) result(new_qd1)
        implicit none 

        double precision, intent(in) :: current_y, current_b, current_m ! Current state 
        double precision :: zL1, zH1, ytempObj_top, ytempObj_bottom, chosen_rho_tom, MU_num, & 
            MU_denom, lmu_top, lmu_bottom, defInd, tom_q, tom_qdR, tom_qdD, marg_return_R, marg_return_D
        double precision :: new_qd1

        integer :: ix, iy ! Loop counters

        double precision, dimension(quad_size) :: z_gc_points, yvec_top, yvec_bottom, xvec_top, xvec_bottom
        double precision, dimension(nDims-1) :: future_info_state, future_borrow_state
        double precision, dimension(nDims) :: tom_qd_state, tom_q_state
        
        zL1 = rrho_z*log(current_y) - stds_out*ssigma_z
        zH1 = rrho_z*log(current_y) + stds_out*ssigma_z

        z_gc_points =  zL1 + ( 1.0 + GCnodes)/2 *(zH1-zL1)


        do iy=1,quad_size 

            future_borrow_state(1) = exp(z_gc_points(iy)) 
            future_borrow_state(2) = current_m*current_b

            future_info_state(1) = future_borrow_state(1)
            future_info_state(2) = APol(future_borrow_state)

            chosen_rho_tom = infoPol(future_info_state)

            tom_q_state(1) = future_borrow_state(1)
            tom_q_state(2) = future_info_state(2)

            defInd = defFun(exp(z_gc_points(iy)),current_b*current_m,current_m)

            tom_qd_state(1) = exp(z_gc_points(iy))
            tom_qd_state(2) = current_b*current_m
            tom_qd_state(3) = current_m

            tom_qdR = qdRN(tom_qd_state)

            tom_qd_state(2) = current_b

            tom_qdD = qdRN(tom_qd_state)


            if(chosen_rho_tom > 0.01) then 

                do ix=1,quad_size 

                    tom_q_state(3) = x_GCnodes(ix)

                    marg_return_R = current_m*( (1.0-defInd)*(llambda + coup + (1.0-llambda)*qRN(tom_q_state) ) + defInd*tom_qdR)
                    marg_return_D = tom_qdD

                    lmu_bottom = 1.0
                    
                    lmu_top = pphi*marg_return_R + (1.0-pphi)*marg_return_D

                    xvec_top(ix) = GCweights(ix)*lmu_top*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                    xvec_bottom(ix) = GCweights(ix)*lmu_bottom*sqrt(1.0-GCnodes(ix)**2.0)* &
                            normPDF( (x_GCnodes(ix) - norm_mu_est)/norm_sigma_est )/norm_sigma_est/pdf_adjust 

                end do 

                ytempObj_top = (xH-xL)/2.0*sum(xvec_top)
                ytempObj_bottom = (xH-xL)/2.0*sum(xvec_bottom)

            else 

                tom_q_state(3) = x_GCnodes(FLOOR(quad_size/2.0)+1)

                marg_return_R = current_m*((1.0-defInd)*(llambda + coup + (1.0-llambda)*qRN(tom_q_state) ) + defInd*tom_qdR )              
                marg_return_D = tom_qdD

                ytempObj_bottom = 1.0
                    
                ytempObj_top = pphi*marg_return_R + (1.0-pphi)*marg_return_D

            end if

            yvec_top(iy) = GCweights(iy)*ytempObj_top*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 
            yvec_bottom(iy) = GCweights(iy)*ytempObj_bottom*sqrt(1.0-GCnodes(iy)**2.0)* & 
                normPDF( (z_gc_points(iy) - rrho_z*log(current_y))/ssigma_z )/ssigma_z/pdf_adjust 

        end do 

        MU_num = (zH1-zL1)/2.0*sum(yvec_top)
        MU_denom = (zH1-zL1)/2.0*sum(yvec_bottom)

        new_qd1 = MU_num/((1.0+r)*MU_denom)

    end function new_qdRN


    ! This subroutine draws Gauss-Chebyshev nodes and weights for integration
    subroutine GaussChebyshevQuad(wweights, nnodes)
        implicit none 

        double precision, intent(inout), dimension(quad_size) :: wweights, nnodes

        integer :: i ! Counter for loop

        do i=1,quad_size
            nnodes(i) = -cos( (2*i-1)/(2*DBLE(quad_size))*ppi)
            wweights(i) = ppi/DBLE(quad_size)
        end do

    end subroutine GaussChebyshevQuad

    ! Generate a randomly permuted list from 1 to number_of_values
    ! Used to generate training inputs drawn uniformly when multiple 
    ! uniform dimensions are required
    subroutine scramble( number_of_values, myArray )
        implicit none 

        integer, intent(in)  :: number_of_values
        integer, intent(inout), dimension(number_of_values)  :: myArray

        integer :: i, n, m, k, j, itemp
        double precision :: u 
        
        myArray=[(i,i=1,number_of_values)]
        n=1; m=number_of_values
        do k=1,2
            do i=1,m
                call random_number(u)
                j = n + FLOOR((m+1-n)*u)
                itemp=myArray(j); myArray(j)=myArray(i); myArray(i)=itemp
            end do
        end do

    end subroutine scramble

    ! This function draws training inputs (grid points) using a uniform shuffle techqnique
    subroutine drawInputs(trainingInputsX,trainingInputsM)
        implicit none 

        double precision, dimension(nDims,nInputs), intent(inout) :: trainingInputsX, trainingInputsM
        integer, dimension(nInputs) :: b_inds

        integer :: i ! Counter for loop 

        logical :: not_typical ! Boolean to check whether the set is unconditionally typical
        double precision :: x_draw_entropy, z_draw_entropy, x_draw, z_draw

        integer, parameter :: num_draws = 1 ! How many draws we do at any one time
        double precision, dimension(num_draws) :: unif_shock, std_norm_shock ! Used to draw standard normals

        do i=1,bint_thresh
            trainingInputsX(2,i) = bL + (b_thresh-bL)/(DBLE(bint_thresh)-1.0)*DBLE(i-1)
            trainingInputsM(2,i) = bL + (b_thresh-bL)/(DBLE(bint_thresh)-1.0)*DBLE(i-1)
        end do 
        do i=bint_thresh+1,nInputs
            trainingInputsX(2,i) = b_thresh + (bH-b_thresh)/(DBLE(nInputs-bint_thresh)-1.0)*DBLE(i-bint_thresh-1)
            trainingInputsM(2,i) = b_thresh + (bH-b_thresh)/(DBLE(nInputs-bint_thresh)-1.0)*DBLE(i-bint_thresh-1)
        end do 

        ! x and y and drawn to be unconditionally typical
        not_typical = .TRUE.
        do while(not_typical)
            z_draw_entropy = 0.0
            do i=1,nInputs 

                call random_number(unif_shock(1))
                call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)

                z_draw = ssigma_uncond*std_norm_shock(1)
                trainingInputsX(1,i) = exp( z_draw )
                z_draw_entropy = z_draw_entropy - log(normPDF(z_draw/ssigma_uncond)/ssigma_uncond)/nInputs
            end do 

            if(  abs(z_draw_entropy - z_entropy) < typ_thresh ) then 

                not_typical = .FALSE. 

            end if 

        end do 

        not_typical = .TRUE.
        do while(not_typical)
            x_draw_entropy = 0.0
            do i=1,nInputs 

                call random_number(unif_shock(1))
                call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)

                x_draw = norm_mu_est + norm_sigma_est*std_norm_shock(1)
                trainingInputsX(3,i) = x_draw
                x_draw_entropy = x_draw_entropy - log(normPDF((x_draw-norm_mu_est)/norm_sigma_est)/norm_sigma_est)/nInputs

            end do 

            if( abs(x_draw_entropy - x_entropy) < typ_thresh ) then 

                not_typical = .FALSE. 

            end if 

        end do 

        do i=1,nInputs
            trainingInputsM(1,i) = trainingInputsX(1,i) 
            trainingInputsM(3,i) = exp(trainingInputsX(3,i))/(1.0 + exp(trainingInputsX(3,i))) 
        end do 

    end subroutine drawInputs 

    ! This function draws gridpoints to evaluate the value function for convergence
    ! It draws uniformly from the whole state space
    subroutine drawConvergenceInputs(convergenceInputs)
        implicit none 

        double precision, dimension(nDims-1,N_validate), intent(inout) :: convergenceInputs

        integer :: i ! Counter for loop 
        double precision :: u 

        do i=1,N_validate

            call random_number(u)
            convergenceInputs(1,i) = yL + (yH-yL)*u

            call random_number(u)
            convergenceInputs(2,i) = bL + (bH-bL)*u

        end do 

    end subroutine drawConvergenceInputs 

    ! Define some logit-based functions to be used in the solution
    function logit1(x) result(logit1_1)
        implicit none 

        double precision, intent(in) :: x
        double precision :: logit1_1

        logit1_1 = log(x/(1.0-x))

    end function logit1 

    function logit_pdf(x,ssigma,mmu) result(logit_pdf1)
        implicit none 

        double precision, intent(in) :: x, ssigma, mmu
        double precision :: logit_pdf1

        logit_pdf1 = 1.0/(ssigma*sqrt(2.0*ppi))*1.0/(x*(1.0-x))*exp(-(logit1(x)-mmu)**2.0/(2.0*ssigma**2.0) )

    end function logit_pdf 

    function logit_pdf_est(x) result(logit_pdf_est1)
        implicit none 

        double precision, intent(in) :: x
        double precision :: logit_pdf_est1

        logit_pdf_est1 = logit_pdf(x,norm_sigma_est,norm_mu_est)

    end function logit_pdf_est 


    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    ! This is used in the construction of the log-likelihood function 
    ! in the GPR
    function inv(A) result(Ainv)

        implicit none

        double precision, dimension(:,:), intent(in) :: A
        double precision, dimension(size(A,1),size(A,2)) :: Ainv

        double precision, dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
            stop 'Matrix inversion failed!'
        end if

    end function inv

    ! Give the standard normal, i.e., mean=0, stdev = 1
    function normPDF(x) result(normPDF1)
        implicit none 

        double precision :: x
        double precision :: normPDF1

        normPDF1 = 1.0/sqrt(2.0*ppi)*exp(-0.5*x**2.0)

    end function normPDF

    ! This subroutine simulates the model at the current policy/price functions
    ! and fills the "key_moments vector with relevant model moments
    subroutine simulated_moments(key_moments) 

        implicit none 

        include "mkl_vml.f90"

        double precision, intent(inout), dimension(nMoments) :: key_moments
        double precision, dimension(Nsimul+1) :: sim_y, sim_b, sim_m, sim_x, sim_rho, sim_q, & 
            sim_qRN, sim_access, sim_issuing, sim_tb, sim_ann_sprd, sim_by_issuing, sim_tb_share_issuing, sim_ann_sprd_issuing, &
            sim_y_issuing, sim_info_issuing, sim_ann_sprdRN, sim_ann_sprd_issuingRN, sim_DP
        double precision ::  def_freq_quarterly, def_freq_annual, avg_by_repay, avg_tb_share_repay, avg_ann_sprd_repay, & 
            avg_ann_sprd_repayRN, tot_squares_y, tot_squares_tb_share, tot_squares_ann_sprd, tot_y_tb_share, tot_y_ann_sprd, & 
            var_y, var_tb_share, var_ann_sprd, cov_y_tb_share, cov_y_ann_sprd, corr_y_tb_share, corr_y_ann_sprd, &
            info_threshold, attn_share
        integer :: t, num_periods_access, num_periods_issuing, def_counter, high_attn_period
        double precision :: re_entry_shock, drawn_m_logit, mean_norm_x, vol_norm_x
        double precision, dimension(nDims) :: curr_def_state, curr_price_state
        double precision, dimension(nDims-1) :: curr_repay_state, curr_info_state

        integer, parameter :: num_draws = 1 ! How many draws we do at any one time

        double precision, dimension(num_draws) :: unif_shock, std_norm_shock ! To play nicely with MKL, these need to 
        ! be single-valued arrays, technically

        sim_issuing = 0.0 ! Set issuance to zero and only change it when we see issuance/buyback

        sim_y(1) =  1.0 ! start at steady state
        sim_b(1) = 0.0 ! start with no debt
        sim_m(1) = 0.651 ! Start at mean
        sim_access(1) = 1.0 ! start with access to credit markets
        sim_issuing(1) = 1.0 ! Start in a period with issuance taking place
        sim_DP(1) = 0.0 ! No default probability next period initially 

        num_periods_access = 0
        num_periods_issuing = 0
        def_counter = 0

        do t=1,Nsimul 

            call random_number(re_entry_shock)

            if ( (sim_access(t) > 0.5) .OR. ( (sim_access(t) < 0.5) .AND. (re_entry_shock < pphi) ) ) then 
                ! This is the block for when the country is in credit markets or when it has access to a deal
                num_periods_access = num_periods_access + 1

                if(sim_access(t) < 0.5) then 
                    ! Here, it was in default, but got a shot at a deal with lenders with recovery m_t;
                    ! This applies whether or not the terms of the deal are accepted
                    sim_b(t) = sim_m(t)*sim_b(t)
                    sim_access(t) = 1.0 ! Update access indicator since it regained access
                end if 
                
                if ( defFun(sim_y(t), sim_b(t), sim_m(t) ) > 0.5) then
                    ! Here, it defaulted even though it had access at the beginning of the period
                    def_counter = def_counter + 1

                    curr_def_state(1) = sim_y(t)
                    curr_def_state(2) = sim_b(t)
                    curr_def_state(3) = sim_m(t)

                    sim_b(t+1) = sim_b(t)
                    sim_access(t+1) = 0.0 
                    sim_DP(t) = 1.0-pphi ! Chance of remaining in default tomorrow
                    sim_q(t) = qd(curr_def_state)
                    sim_qRN(t) = qdRN(curr_def_state)
                    sim_rho(t) = 0.0 
                    sim_m(t+1) = sim_m(t)

                    ! Draw a random signal using noise
                    call random_number(unif_shock(1))
                    call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)
                    sim_x(t) = norm_mu_est + norm_sigma_est*std_norm_shock(1)

                    sim_tb(t) = 0.0 ! In default, no issuing of new debt but no servicing of old: Y = C

                else 
                    ! Here, it repaid
                    num_periods_issuing = num_periods_issuing + 1
                    sim_issuing(t) = 1.0

                    curr_repay_state(1) = sim_y(t)
                    curr_repay_state(2) = sim_b(t)

                    sim_b(t+1) = Apol(curr_repay_state)
                    sim_access(t+1) = 1.0 
                    
                    curr_info_state(1) = sim_y(t)
                    curr_info_state(2) = sim_b(t+1)

                    sim_rho(t) = infoPol(curr_info_state)


                    ! Draw a random recovery rate using noise (r)
                    call random_number(unif_shock(1)) 
                    call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)

                    drawn_m_logit = norm_mu_est + norm_sigma_est*std_norm_shock(1)

                    sim_m(t+1) = exp(drawn_m_logit)/(1.0+exp(drawn_m_logit))
                    
                    mean_norm_x = norm_mu_est + sim_rho(t)*(drawn_m_logit - norm_mu_est)
                    vol_norm_x = sqrt(1.0 - sim_rho(t)**2.0)*norm_sigma_est


                    ! Now draw a random signal of that rate using implied correlation
                    call random_number(unif_shock(1))
                    call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)
                    
                    sim_x(t) = mean_norm_x + vol_norm_x*std_norm_shock(1)

                    curr_price_state(1) = sim_y(t)
                    curr_price_state(2) = sim_b(t+1)
                    curr_price_state(3) = sim_x(t)

                    sim_q(t) = q(curr_price_state)
                    sim_qRN(t) = qRN(curr_price_state)
                    sim_DP(t) = defProb(curr_price_state)

                    sim_tb(t) = (llambda + coup)*sim_b(t) - sim_q(t)*(sim_b(t+1) - (1.0-llambda)*sim_b(t))

                    ! Set aside objects to compute repayment-only moments
                    sim_ann_sprd_issuing(num_periods_issuing) = ( (llambda + coup + (1.0-llambda)*sim_q(t))/sim_q(t) )**4.0 - (1.0 + r)**4.0
                    sim_ann_sprd_issuingRN(num_periods_issuing) = ( (llambda + coup + (1.0-llambda)*sim_qRN(t))/sim_qRN(t) )**4.0 - (1.0 + r)**4.0
                    sim_by_issuing(num_periods_issuing) = sim_b(t)/sim_y(t)
                    sim_tb_share_issuing(num_periods_issuing) = sim_tb(t)/sim_y(t)
                    sim_y_issuing(num_periods_issuing) = sim_y(t) 
                    sim_info_issuing(num_periods_issuing) = log(1.0/sqrt(1.0 - sim_rho(t)**2.0))

                end if 

            else 
                ! Here, it entered in a state of default and remained there

                curr_def_state(1) = sim_y(t)
                curr_def_state(2) = sim_b(t)
                curr_def_state(3) = sim_m(t)

                sim_b(t+1) = sim_b(t)
                sim_access(t+1) = 0.0 
                sim_DP(t) = 1.0-pphi ! Chance of remaining in default tomorrow
                sim_q(t) = qd(curr_def_state)
                sim_qRN(t) = qdRN(curr_def_state)
                sim_rho(t) = 0.0 
                sim_m(t+1) = sim_m(t) 

                ! Draw a random signal using noise 
                call random_number(unif_shock(1))
                call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)
                sim_x(t) = norm_mu_est + norm_sigma_est*std_norm_shock(1)

                sim_tb(t) = 0.0

            end if

            ! Convert prices to spreads
            sim_ann_sprd(t) = ( (llambda + coup + (1.0-llambda)*sim_q(t))/sim_q(t) )**4.0 - (1.0 + r)**4.0
            sim_ann_sprdRN(t) = ( (llambda + coup + (1.0-llambda)*sim_qRN(t))/sim_qRN(t) )**4.0 - (1.0 + r)**4.0

            ! Now update the output process
            call random_number(unif_shock(1))
            call vdcdfnorminv(num_draws,unif_shock,std_norm_shock)
            sim_y(t+1) = min( max( exp( rrho_z*log(sim_y(t)) + ssigma_z*std_norm_shock(1)), yL), yH)

        end do 

        ! With simulations in hand, we begin to compute key moments
        avg_by_repay = sum(sim_by_issuing(1:num_periods_issuing))/DBLE(num_periods_issuing)
        avg_ann_sprd_repay = sum(sim_ann_sprd_issuing(1:num_periods_issuing))/DBLE(num_periods_issuing)
        avg_ann_sprd_repayRN = sum(sim_ann_sprd_issuingRN(1:num_periods_issuing))/DBLE(num_periods_issuing)
        avg_tb_share_repay = sum(sim_tb_share_issuing(1:num_periods_issuing))/DBLE(num_periods_issuing)

        ! First the default frequency
        def_freq_quarterly = DBLE(def_counter)/DBLE(num_periods_access)
        def_freq_annual = 1.0 - (1.0-def_freq_quarterly)**4.0

        ! Next those moments related only to repayment periods
        info_threshold =  min(info_cap,maxval(sim_info_issuing(1:num_periods_issuing)))*0.5

        tot_squares_y = 0.0 
        tot_squares_tb_share = 0.0 
        tot_squares_ann_sprd = 0.0
        tot_y_tb_share = 0.0
        tot_y_ann_sprd = 0.0
        high_attn_period = 0
        do t=1,num_periods_issuing

            tot_squares_y = tot_squares_y + (sim_y_issuing(t) - 1.0)**2.0
            tot_squares_tb_share = tot_squares_tb_share + (sim_tb_share_issuing(t) - avg_tb_share_repay)**2.0
            tot_squares_ann_sprd = tot_squares_ann_sprd + (sim_ann_sprd_issuing(t) - avg_ann_sprd_repay)**2.0
            tot_y_tb_share = tot_y_tb_share + (sim_tb_share_issuing(t) - avg_tb_share_repay)*(sim_y_issuing(t) - 1.0)
            tot_y_ann_sprd = tot_y_ann_sprd + (sim_ann_sprd_issuing(t) - avg_ann_sprd_repay)*(sim_y_issuing(t) - 1.0)

            if(sim_info_issuing(t) > info_threshold) high_attn_period = high_attn_period + 1

        end do 

        var_y = tot_squares_y/(DBLE(num_periods_issuing) - 1.0)
        var_tb_share = tot_squares_tb_share/(DBLE(num_periods_issuing) - 1.0)
        var_ann_sprd = tot_squares_ann_sprd/(DBLE(num_periods_issuing) - 1.0)
        cov_y_tb_share = tot_y_tb_share/(DBLE(num_periods_issuing) - 1.0)
        cov_y_ann_sprd = tot_y_ann_sprd/(DBLE(num_periods_issuing) - 1.0)

        if(sqrt(var_y*var_tb_share) > .00001) then 
            corr_y_tb_share = cov_y_tb_share/sqrt(var_y*var_tb_share) 
        else 
            corr_y_tb_share = 0.0
        end if 

        if(sqrt(var_y*var_ann_sprd) > .00001) then 
            corr_y_ann_sprd = cov_y_ann_sprd/sqrt(var_y*var_ann_sprd)
        else 
            corr_y_ann_sprd = 0.0
        end if 

        attn_share = DBLE(high_attn_period)/DBLE(num_periods_issuing)

        ! Stack the key moments into the relevant vector and pass it back
        key_moments(1) = avg_by_repay ! Average debt-to-GDP in repayment
        key_moments(2) = def_freq_annual ! Annualized default frequency
        key_moments(3) = avg_ann_sprd_repay ! Average annualized spreads
        key_moments(4) = sqrt(var_ann_sprd) ! Annualized spread volatility
        key_moments(5) = attn_share ! Attention fraction 
        key_moments(6) = corr_y_ann_sprd ! cyclicality of spreads
        key_moments(7) = corr_y_tb_share ! cyclicality of the trade share 
        key_moments(8) = avg_ann_sprd_repayRN/avg_ann_sprd_repay ! risk-neutral spread share

        open(3,file="simulations.txt")

            do t=1,Nsimul 

                write(3,"(11F9.6)") sim_y(t), sim_b(t), sim_m(t), sim_x(t), sim_access(t), sim_q(t), sim_ann_sprd(t), & 
                    sim_rho(t), sim_ann_sprdRN(t), sim_issuing(t), sim_DP(t)

            end do 

        close(3)

        open(3,file="simulations_info.txt")

            write(3,*) Nsimul 
            write(3,*) num_periods_access 
            write(3,*) num_periods_issuing

        close(3)

    end subroutine simulated_moments


    ! Compute the determinant of a matrix. Used in the construction of the 
    ! log-likelihood in the GPR
    function DETERMINANT(aa) result(det)
        implicit none 

        double precision :: aa(:,:), det
        double precision :: tmp,c(size(aa,dim=1),size(aa,dim=2))
        double precision :: max
	    integer :: i,j,k,l,m,num(size(aa,dim=1)),n

	    n=size(aa,dim=1)
	    det=1.	
	    do k=1,n
		    max=aa(k,k);num(k)=k;
		    do i=k+1,n 
			    if(abs(max)<abs(aa(i,k))) then
				    max=aa(i,k)
				    num(k)=i
			    endif
		    enddo
		    if (num(k)/=k) then
			    do l=k,n 
				    tmp=aa(k,l)
			    	aa(k,l)=aa(num(k),l)
				    aa(num(k),l)=tmp
			    enddo
			    det=-1.*det
		    endif
		    do m=k+1,n
			    c(m,k)=aa(m,k)/aa(k,k)
			    do l=k,n 
				    aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
			    enddo
		    enddo !There we made matrix triangular!	
	    enddo

	    do i=1,n
	        det=det*aa(i,i)
	    enddo
	
    end function DETERMINANT


    ! These are the likelihood functions for the different 
    ! GPs used in the solution. There is one for each 
    ! equilibrium object
    function negLL_info(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_I,1)

    end function negLL_info

    function negLL_def(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_def,2)

    end function negLL_def

    function negLL_q(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_q,2)

    end function negLL_q

    function negLL_qd(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_qd,3)

    end function negLL_qd

    function negLL_qRN(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_qRN,2)

    end function negLL_qRN

    function negLL_qdRN(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_qdRN,3)

    end function negLL_qdRN

    function negLL_VR(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_VR,1)

    end function negLL_VR

    function negLL_A(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_A,1)

    end function negLL_A

    function negLL_VD(log_hyper_params) result(negLL1)
        implicit none

        double precision :: negLL1
        double precision, dimension(2) :: log_hyper_params

        negLL1 = negLL(log_hyper_params,tset_VD,3)

    end function negLL_VD

    ! The log-likelihood function for given observations t_set and a given set of 
    ! log-hyper-parameters.
    ! type_ind tells us whether to do 2D (1), 3D with x (2), or 3D with m (3)
    function negLL(log_hyper_params,t_set,type_ind) result(negLL1)
        implicit none 

        double precision :: negLL1, temp_num
        double precision, dimension(2), intent(in) :: log_hyper_params
        double precision, dimension(2) :: hyper_params
        double precision, dimension(nInputs), intent(in) :: t_set 
        double precision, dimension(nInputs,nInputs) :: currK, eyeN, Ssigma, log_Ssigma

        integer :: ik1, i1, type_ind, inf_flag

        if( (log_hyper_params(1) > kernSSLogLB) .AND. (log_hyper_params(1) < kernSSLogUB) & 
            .AND. (log_hyper_params(2) > kernLogLB) .AND. (log_hyper_params(2) < kernLogUB) ) then

            hyper_params = 10.0**log_hyper_params

            call createK(currK,hyper_params,type_ind)

            ! Make a quick identity matrix
            eyeN = 0.0 
            do ik1=1,nInputs 
                eyeN(ik1,ik1) = 1.0 
            end do 

            Ssigma = currK + sn_star**2.0*eyeN

            temp_num = DETERMINANT(Ssigma)

            if(temp_num > 0.0) then 
                negLL1 = log(temp_num) + sum(matmul(t_set,inv(Ssigma))*t_set)
            else 
                negLL1 = 1.0e10
            end if 

        else ! If it violates the bounds, return a big number
            negLL1 = 1.0e10 

        end if 

    end function negLL

    ! This creates the covariance matrix across the training inputs for a given 
    ! set of hyperparameters.
    ! type_ind tells us whether to do 2D (1), 3D with x (2), or 3D with m (3)
    subroutine createK(currK,hyper_params,type_ind)
        implicit none

        double precision, dimension(nInputs,nInputs), intent(inout) :: currK
        double precision, dimension(2), intent(in) :: hyper_params
        integer, intent(in) :: type_ind 
        integer :: ik1, ik2

        do ik1=1,nInputs 
            do ik2=1,ik1 

                if(type_ind .EQ. 1) then 
                    currK(ik1,ik2) = kernelFunc(scaledGridPointsX(1:2,ik1),scaledGridPointsX(1:2,ik2),hyper_params,2)
                elseif(type_ind .EQ. 2) then 
                    currK(ik1,ik2) = kernelFunc(scaledGridPointsX(:,ik1),scaledGridPointsX(:,ik2),hyper_params,3)
                else 
                    currK(ik1,ik2) = kernelFunc(scaledGridPointsM(:,ik1),scaledGridPointsM(:,ik2),hyper_params,3)
                end if 

                if(ik1 .NE. ik2) then 
                    currK(ik2,ik1) = currK(ik1,ik2)
                end if 

            end do 
        end do 

    end subroutine createK


    ! The kernel function
    ! x and xprime are the inputs of interest and ttheta_vec is the coefficient matrix
    ! The parameter vector here is comprised of two elements: signal strength and lengthscale,
    ! the latter of which applies to each dimension equally
    function kernelFunc(x,xprime,ttheta_vec,dim_size) result(kernelFunc_1)
        implicit none

        integer, intent(in) :: dim_size
        double precision, intent(in), dimension(2) :: ttheta_vec 
        double precision, intent(in), dimension(dim_size) :: x, xprime
        double precision :: kernelFunc_1
    
        double precision :: s, scaled_d ! s = signal_size
        double precision :: h ! for bringing scale into exponent
        double precision, dimension(nDims) :: l


        s = ttheta_vec(1) ! Signal strength parameter
        l = ttheta_vec(2) ! Scale length parameter

        h = log(s) ! h = log(s) is s=sigma^2; h = 2*log(s) is s=sigma

        ! Scale each dimension using the scale-length l to compute distance between x and xprime
        scaled_d = sqrt( sum( ((x - xprime)/l)**2.0) ) 

        ! Now apply the SE formula
        kernelFunc_1 = exp(-.5*scaled_d + h)

    end function kernelFunc

    ! This function delivers the GPR posterior mode, i.e., functional approximation
    ! at a point x using optimal interpolation from other known points xx and coefficients
    ! AA and kernel parameters kern_params
    function GPR_approx(x,AA,kern_params,xx,dim_size) result(GPR_approx1)
        implicit none 

        integer, intent(in) :: dim_size 
        double precision, intent(in), dimension(dim_size) :: x ! current point
        double precision, intent(in), dimension(dim_size,nInputs) :: xx ! grid
        double precision, intent(in), dimension(2) :: kern_params ! hyperparameters 
        double precision, intent(in), dimension(nInputs) :: AA ! GPR Coefficients
        double precision :: GPR_approx1

        integer :: i ! Counter for loop

        GPR_approx1 = 0.0

        do i=1,nInputs
            GPR_approx1 = GPR_approx1 + AA(i)*kernelFunc(x,xx(:,i),kern_params,dim_size)
        end do

    end function GPR_approx    

end module aux_functions

