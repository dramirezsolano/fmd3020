
module mod_mfp
    !
    ! This module handles the calculation of the distance to the next interaction.
    !
    use iso_fortran_env, only: int32, real64
    
    use mod_rng

implicit none

contains

    real(kind=real64) function mfp(sigma)
        !
        ! This function samples the free flight, the path length of the particle 
        ! before suffering an interaction. The analytical inversion method is implemented.
        !
        real(kind=real64), intent(in) :: sigma ! total interaction cross-section (cm-1)
        real(kind=real64) :: rnno

        rnno = rng_set()
        mfp = -log(1.0_real64-rnno)/sigma
        
    end function mfp

    real(kind=real64) function et_mfp(sigma,c,u)
        !
        ! This function samples the free flight, the path length of the particle 
        ! before suffering an interaction (biasing the pdf through Exponencial Transform
        ! as a Variance Reduction Technique). The analytical inversion method is implemented.
        !
        real(kind=real64), intent(in) :: sigma ! total interaction cross-section (cm-1)
        real(kind=real64), intent(in) :: c
        real(kind=real64), intent(in) :: u
        real(kind=real64) :: rnno

        rnno = rng_set()

        ! mfp due to the modified pdf
        et_mfp = -log(1.0_real64-rnno)/(sigma-(c*u))
        
    end function et_mfp

end module mod_mfp