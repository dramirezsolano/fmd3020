module mod_mfp
    !
    ! El modulo contiene dos metodos de transporte: vuelo libre (pl) y numero de 
    ! caminos libres medio (mfp). Ademas, contiene el calculo del parametro b y bm,
    ! requeridos para el metodo mfp.
    !
    ! Modulo basado en (Doerner, 2019)
    !
    use mod_rng

implicit none

contains

    real function pl(sigma) ! (Doerner, 2019)
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
        ! Funcion que muestrea el metodo mfp. (Haghighat, 2014)
        !
        real, intent(in) :: sigma  ! total interaction cross-section (cm-1)
        real, intent(in) :: bj     ! Parametro b = -log(1-rnno)
        real, intent(in) :: bm     ! Parametros bm para cada region. bm = sum(sigma(i)*r(i))

        mfp = (bj-bm)/sigma
        
    end function mfp

    real function bk()
        ! 
        ! Funcion que calcula el parametro b para cada transporte
        !
        real :: rnno ! numero aleatorio

        rnno = rng_set()
        bk = -log(1.0-rnno)

    end function bk

    function bm_calc(u, nreg, ir, xbound, x, sigma_t)
        !
        ! Funcion que calcula el parametro bm para cada region
        !

        ! Parametros de la geometria
        integer, intent(in) :: nreg, ir
        real, dimension(nreg), intent(in) :: xbound, sigma_t
        real, intent(in) :: u, x

        integer i
        real, dimension(nreg) :: r           ! Longitud de trayectoria en cada region en la direccion de movimiento
        real, dimension(0:nreg+1) :: bm      ! vector que almacena los bm calculados
        real, dimension(0:nreg+1) :: bm_calc ! salida de la funcion como vector

        ! Calculo en direccion positiva
        if(u .gt. 0.0) then
            bm = 0.0

            ! Calculo de longitud de trayectoria
            if(ir == nreg) then
                r(ir) = (xbound(ir+1)-x)/u
            else
                r(ir) = (xbound(ir+1)-x)/u
                r_recalc_xpos: do i = (ir+1),nreg,1
                    r(i) = (xbound(i+1)-x)/u
                enddo r_recalc_xpos
            endif
            
            ! Calculo de los bm
            if(ir == nreg) then
                bm(ir) = sigma_t(ir)*r(ir)
            else
                bm(ir) = sigma_t(ir)*r(ir)
                bm_recalc_xpos: do i = (ir+1),nreg,1
                    bm(i) = sigma_t(i)*(r(i)-r(i-1)) + bm(i-1)
                enddo bm_recalc_xpos
            endif

        ! Calculo en direccion negativa
        else if(u .lt. 0.0) then
            bm = 0.0

            ! Calculo de longitud de trayectoria
            if(ir == 1) then
                r(ir) = (xbound(ir)-(x))/u
            else
                r(ir) = (xbound(ir)-(x))/u
                r_recalc_xneg: do i = (ir-1),1,-1
                    r(i) = (xbound(i)-(x))/u
                enddo r_recalc_xneg
            endif

            ! Calculo de los bm
            if(ir == 1) then
                bm(ir) = sigma_t(ir)*r(ir)
            else
                bm(ir) = sigma_t(ir)*r(ir)
                bm_recalc_xneg: do i = (ir-1),1,-1
                    bm(i) = sigma_t(i)*(r(i)-r(i+1)) + bm(i+1)
                enddo bm_recalc_xneg
            endif
        endif

        ! Retorno de los valores de bm
        bm_calc = bm

    end function bm_calc

end module mod_mfp