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

    real function mfp(sigma)
        !
        ! This function samples the number of the mean-free-paths (mfp), the path length of the particle 
        ! before suffering an interaction. The analytical inversion method is implemented.
        !
        real, intent(in) :: sigma ! total interaction cross-section (cm-1)
        real :: rnno

        rnno = rng_set()
        mfp = -log(1.0-rnno)/sigma
        
    end function mfp

end module mod_mfp