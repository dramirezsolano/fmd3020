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
        mfp = -(bj+bm)/sigma
        
    end function mfp

    real function bk()
        ! hjj
        real :: rnno

        rnno = rng_set()
        bk = -log(1.0-rnno)

    end function bk

    function bm_calc(u, nreg, ir, xthick, x, sigma_t)

        integer, intent(in) :: nreg, ir
        real, dimension(nreg), intent(in) :: xthick, sigma_t
        real, intent(in) :: u, x

        integer i
        real, dimension(nreg) :: bm
        real, dimension(nreg) :: bm_calc

        if(u .gt. 0.0) then
            bm = 0.0

            if(ir == nreg) then
                bm(ir) = sigma_t(ir)*((sum(xthick(1:ir)))-x)
            else
                bm(ir) = sigma_t(ir)*((sum(xthick(1:ir)))-x)

                bm_recalc_xpos: do i = (ir+1),nreg
                    bm(i) = sigma_t(i)*xthick(i) + bm(i-1)
                enddo bm_recalc_xpos
            endif
        else
            bm = 0.0                

            if(ir == 1) then
                bm(ir) = sigma_t(ir)*x
            else
                bm(ir) = sigma_t(ir)*(x-(sum(xthick(1:(ir-1)))))
                
                bm_recalc_xneg: do i = (ir-1),1
                    bm(i) = sigma_t(i)*xthick(i) + bm(i-1)
                enddo bm_recalc_xneg
            endif
        endif

        bm_calc = bm

    end function bm_calc

end module mod_mfp