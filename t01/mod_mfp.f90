module mod_mfp
    !
    ! This module handles the calculation of the distance to the next interaction.
    !
    use mod_rng

implicit none

contains

    real function pl(sigma)
        !
        ! This function samples the free flight, the path length of the particle 
        ! before suffering an interaction. The analytical inversion method is implemented.
        !
        real, intent(in) :: sigma ! total interaction cross-section (cm-1)
        real :: rnno

        rnno = rng_set()
        pl = -log(1.0-rnno)/sigma
        
    end function pl

    real function mfp(sigma,bj,bm)
        !
        ! This function samples the number of the mean-free-paths (mfp), the path length of the particle 
        ! before suffering an interaction. The analytical inversion method is implemented.
        !
        real, intent(in) :: sigma ! total interaction cross-section (cm-1)
        real, intent(in) :: bj
        real, intent(in) :: bm
        ! real :: rnno

        ! rnno = rng_set()
        ! mfp = -(log(1.0-rnno)+bm)/sigma
        mfp = (bj-bm)/sigma
        ! write(*,'(A,F15.5)') 'mfp : ', mfp
        
    end function mfp

    real function bk()
        ! hjj
        real :: rnno

        rnno = rng_set()
        bk = -log(1.0-rnno)

    end function bk

    function bm_calc(u, nreg, ir, xbound, x, sigma_t)

        integer, intent(in) :: nreg, ir
        real, dimension(nreg), intent(in) :: xbound, sigma_t
        real, intent(in) :: u, x

        integer i
        real, dimension(nreg) :: r
        real, dimension(0:nreg+1) :: bm
        real, dimension(0:nreg+1) :: bm_calc

        if(u .gt. 0.0) then
            bm = 0.0

            if(ir == nreg) then
                r(ir) = (xbound(ir+1)-x)/u
            else
                r(ir) = (xbound(ir+1)-x)/u
                r_recalc_xpos: do i = (ir+1),nreg,1
                    r(i) = (xbound(i+1)-x)/u
                enddo r_recalc_xpos
            endif
            if(ir == nreg) then
                bm(ir) = sigma_t(ir)*r(ir)
            else
                bm(ir) = sigma_t(ir)*r(ir)
                bm_recalc_xpos: do i = (ir+1),nreg,1
                    bm(i) = sigma_t(i)*(r(i)-r(i-1)) + bm(i-1)
                enddo bm_recalc_xpos
            endif

        else if(u .lt. 0.0) then
            bm = 0.0

            if(ir == 1) then
                r(ir) = (xbound(ir)-(x))/u
                ! r(ir) = (x)/u
            else
                r(ir) = (xbound(ir)-(x))/u
                ! r(ir) = (x)/u
                ! write(*,'(A, F15.5)') 'r ir:', r(ir)
                r_recalc_xneg: do i = (ir-1),1,-1
                    r(i) = (xbound(i)-(x))/u
                    ! r(i) = (x)/u
                    ! write(*,'(A, F15.5)') 'r :', r(i)
                enddo r_recalc_xneg
            endif
            ! r = -r
            if(ir == 1) then
                bm(ir) = sigma_t(ir)*r(ir)
            else
                bm(ir) = sigma_t(ir)*r(ir)
                bm_recalc_xneg: do i = (ir-1),1,-1
                    bm(i) = sigma_t(i)*(r(i)-r(i+1)) + bm(i+1)
                enddo bm_recalc_xneg
            endif
            ! do i = 1,ir
                ! write(*,'(A, F15.5)') 'bm_pos (): ', bm_pos(i)
                ! write(*,'(A, F15.5)') 'bm (): ', bm(i)
            ! enddo
        endif

        ! if(u .gt. 0.0) then
        !     bm = 0.0
        !     ! bm(1) = sigma_t(1)*xthick(1)                
        !     ! bm_recalc_xpos: do i = 2,nreg
        !     !     bm(i) = sigma_t(i)*xthick(i) + bm(i-1)
        !     ! enddo bm_recalc_xpos
        !     if(ir == nreg) then
        !         bm(ir) = sigma_t(ir)*((sum(xthick(1:ir)))-x)
        !     else
        !         bm(ir) = sigma_t(ir)*((sum(xthick(1:ir)))-x)

        !         bm_recalc_xpos: do i = (ir+1),nreg
        !             bm(i) = sigma_t(i)*xthick(i) + bm(i-1)
        !         enddo bm_recalc_xpos
        !     endif
        ! else if(u .lt. 0.0) then
        !     bm = 0.0
        !     ! bm(nreg) = sigma_t(nreg)*xthick(nreg)                
        !     ! bm_recalc_xneg: do i = (nreg-1),1
        !     !     bm(i) = sigma_t(i)*xthick(i) + bm(i+1)
        !     ! enddo bm_recalc_xneg
        !     if(ir == 1) then
        !         bm(ir) = sigma_t(ir)*x
        !     else
        !         bm(ir) = sigma_t(ir)*(x-(sum(xthick(1:(ir-1)))))
                
        !         bm_recalc_xneg: do i = (ir-1),1
        !             bm(i) = sigma_t(i)*xthick(i) + bm(i+1)
        !         enddo bm_recalc_xneg
        !     endif
        ! endif
        ! do i = 1,nreg
        !     write(*,'(F15.5)') bm_calc
        ! enddo
        bm_calc = bm

    end function bm_calc

end module mod_mfp