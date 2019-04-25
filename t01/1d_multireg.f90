program shield_1d
    !
    ! Programa que simula el transporte de partículas a través de un material con diferentes
    ! secciones eficaces de interacción, mediante dos métodos de transporte: vuelo libre (pl)
    ! y numero de caminos libres medio (mfp), basado en el código 1d_shield.f90 version multi
    ! region (Doerner, 2019).
    !
    ! Prof.: Edgardo Doerner, edoerner@fis.puc.cl
    ! Autores: Jessica Hernandez
    !          Daniel Ramirez
    !
    use mod_rng ! (Doerner, 2019)
    use mod_mfp ! (Doerner, 2019), modificaciones agregadas
    use mod_scatt ! (Doerner, 2019)

    implicit none

    ! Parametros geometricos
    integer, parameter :: nreg = 6           ! numero de regiones
    real, dimension(nreg) :: xthick = 1.0    ! espesor de cada region (cm)
    real, dimension(nreg+1) :: xbounds = 0.0      ! bordes de cada region (cm)

    ! Parametros de transporte
    ! real, parameter :: c = 0.5 ! Caso 1
    ! real, parameter :: c = 0.7 ! Caso 2
    real, parameter :: c = 0.9 ! Caso 3
    real, parameter :: sigma_tin = 1.0
    real, dimension(0:nreg+1) :: bm_pos = 0.0  ! almacena los bm para u > 0.0
    real, dimension(0:nreg+1) :: bm_neg = 0.0   ! almacena los bm para u < 0.0
    real, dimension(nreg) :: sigma_t      ! seccion eficaz total de interaccion (cm-1)
    real, dimension(nreg) :: sigma_s      ! seccion eficaz de dispersion (cm-1)
    real, dimension(nreg) :: sigma_a      ! seccion eficaz de absorcion (cm-1)

    ! Parametros de la particula
    real, parameter :: xin = 0.0        ! posicion inicial (cm)
    real, parameter :: uin = 1.0        ! direccion inicial
    real, parameter :: wtin = 1.0       ! peso estadistico
    integer, parameter :: irin = 1      ! region inicial

    integer :: ir
    real :: x, u, wt

    ! Parametros de la simulacion
    integer, parameter :: nhist = 100000    ! numero de historias
    character (len = 10) :: mode = 'pl'              ! metodo de transporte
    ! character (len = 10) :: mode = 'mfp'              ! metodo de transporte

    ! Scoring variables
    real, dimension(0:nreg+1) :: score = 0.0    ! score(0) : reflexion
                                                ! score(1:nreg) : absorcion
                                                ! score(nreg+1) : transmision

    integer :: i, ihist   
    logical :: pdisc          ! marcador para descartar la particula   
    integer :: irnew          ! indice de la nueva region
    real :: pstep, dist, b    ! distancia a la siguiente interaccion/ distancia al borde proximo/ parametro b
    real :: rnno 
    real :: start_time, end_time

    call cpu_time(start_time)

    ! Inicializacion RNG
    call rng_init(20180815)

    ! Calculo de las secciones eficaces de cada region
    sigma_t(1) = sigma_tin
    sigma_s(1) = c*sigma_t(1)
    sigma_a(1) = sigma_t(1) - sigma_s(1)
    do i = 2,nreg
        sigma_t(i) = sigma_t(i-1) - 0.1
        sigma_s(i) = c*sigma_t(i)
        sigma_a(i) = sigma_t(i) - sigma_s(i)
    enddo

    ! Impresion en pantalla de parametros de la simulacion
    write(*,'(A, I15)') 'Number of histories ', nhist
    do i = 1,nreg
        write(*,'(A, I5, A)') 'For region ', i, ':'
        write(*,'(A, F15.5)') 'Total interaction cross-section (cm-1): ', sigma_t(i)
        write(*,'(A, F15.5)') 'Absoption interaction cross-section (cm-1): ', sigma_a(i)
        write(*,'(A, F15.5)') 'Scattering interaction cross-section (cm-1): ', sigma_s(i)
        write(*,'(A, F15.5)') '1D shield thickness (cm): ', xthick(i)
    enddo
   
    ! Inicializacion de la geometria
    write(*,'(A)') 'Region boundaries (cm):'
    do i = 1,nreg
        xbounds(i+1) = xbounds(i) + xthick(i)
        write(*,'(F15.5)') xbounds(i+1)
    enddo
    
    ! Iniciamos el metodo de transporte elegido
    ! El primero corresponde al enfoque vuelo libre (pl)    
    if(mode == 'pl') then
        write(*,'(A)') 'MODO: path-length'

        ihist_loop: do ihist = 1,nhist
            
            ! Inicializacion de la historia
            x = xin
            u = uin
            wt = wtin
            ir = irin
            
            ! Marcador para no descartar la particula
            pdisc = .false.

            ! Inicia el proceso de transporte
            particle_loop: do

                ptrans_loop: do

                    ! Distancia a la siguiente interaccion
                    if(ir == 0 .or. ir == nreg+1) then
                        ! pstep fuera de la geometría
                        pstep = 1.0E8
                    else
                        pstep = pl(sigma_t(ir))
                    endif

                    ! Actualizamos la region
                    irnew = ir

                    ! Verificacion de la continuidad de la particula en la geometria
                    if(u .lt. 0.0) then
                        if(irnew == 0) then
                            ! Abandono la geometria, particula descartada
                            pdisc = .true.
                        else
                            ! Traslado de la particula a la siguiente frontera
                            dist = (xbounds(ir) - x)/u

                            ! Verificacion de si salio de la region
                            if(dist < pstep) then
                                pstep = dist                      
                                irnew = irnew - 1
                            endif
                        endif
                        
                    else if(u .gt. 0.0) then
                        if(irnew == nreg+1) then
                            ! Abandono la geometria, particula descartada
                            pdisc = .true.
                        else
                            ! Traslado de la particula a la frontera anterior
                            dist = (xbounds(ir+1) - x)/u

                            ! Verificacion de si salio de la region
                            if(dist < pstep) then
                                pstep = dist
                                irnew = irnew + 1
                            endif
                        endif
                        
                    endif

                    ! Particula descartada, sale del transporte
                    if(pdisc .eqv. .true.) then
                        exit
                    endif

                    ! Transportamos a la nueva posicion
                    x = x + pstep*u

                    ! Si la particula permanece en la misma region, interactua
                    if(ir == irnew) then
                        exit
                    else
                        ! Actualizacion de la region actual
                        ir = irnew
                    endif

                enddo ptrans_loop
                
                if(pdisc .eqv. .true.) then
                    ! Particula descartada, se incrementa marcador. Se termina su historia
                    score(ir) = score(ir) + wt
                    exit
                endif

                ! Determinacion del tipo de interaccion
                rnno = rng_set()
                if(rnno .le. sigma_a(ir)/sigma_t(ir)) then
                    ! Marcador de absorcion de particulas
                    score(ir) = score(ir) + wt
                    exit
                else
                    ! Particula dispersada, continua su historia en la nueva direccion
                    u = scatt(u)
                endif
                
            enddo particle_loop

        enddo ihist_loop

        ! Impresion de resultados en pantalla
        write(*,'(A,F15.5)') 'Reflection : ', score(0)/real(nhist)
        write(*,'(A,F15.5)') 'Absorption : ', sum(score(1:nreg))/real(nhist)
        write(*,'(A,F15.5)') 'Transmission : ', score(nreg+1)/real(nhist)

        ! Tiempo de simulacion
        call cpu_time(end_time)
        write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time

! Enfoque de numero de caminos libres medio (mfp)
! (son extrapolables los comentarios de la seccion anterior)
    elseif (mode == 'mfp') then
        write(*,'(A)') 'MODO: mfp'

        ihist_loop_mfp: do ihist = 1,nhist

            x = xin
            u = uin
            wt = wtin
            ir = irin
            
            pdisc = .false.

            particle_loop_mfp: do
                
                pmfp_loop: do

                    if(ir == 0 .or. ir == nreg+1) then
                        pstep = 1.0E8
                    else
                        b = bk() ! calculo del parametro b = -log(1-rnno)
                    endif                
                    
                    irnew = ir

                    ! Particula va en direccion negativa
                    if (u .lt. 0.0) then

                        ! Calculo de los bm para las regiones siguientes a su direccion
                        ! de movimiento
                        bm_neg = bm_calc(u, nreg, ir, xbounds, x, sigma_t)

                        if(irnew == 0) then
                            pdisc = .true.
                            pstep = 1.0E8
                        else
                            ! Comprobacion de bm-1 > b >= bm para conocer la region
                            if(b > bm_neg(1)) then  ! particula sale de la geometria
                                irnew = 0
                                pdisc = .true.
                                pstep = 1.0E8
                            elseif(b <= bm_neg(irnew)) then ! particula en la misma region
                                pstep = mfp(sigma_t(irnew),b,bm_neg(irnew+1))
                            else
                                bm_loop_xneg: do i = irnew,2,-1 ! loop para conocer la region
                                    if(bm_neg(i-1) > b .and. b >= bm_neg(i)) then
                                        irnew = i-1
                                        pstep = mfp(sigma_t(irnew),b,bm_neg(irnew+1))
                                        exit
                                    endif
                                enddo bm_loop_xneg
                            endif                                        
                        endif

                    ! Particula va en direccion positiva
                    else if (u .gt. 0.0) then

                        ! Calculo de los bm para las regiones siguientes a su direccion
                        ! de movimiento
                        bm_pos = bm_calc(u, nreg, ir, xbounds, x, sigma_t)

                        if(irnew == nreg+1) then
                            pdisc = .true.
                            pstep = 1.0E8
                        else
                            ! Comprobacion de bm-1 <= b < bm para conocer la region
                            if(b > bm_pos(nreg)) then   ! particula sale de la geometria
                                irnew = nreg+1
                                pdisc = .true.
                                pstep = 1.0E8
                            elseif(b <= bm_pos(irnew)) then ! particula en la misma region
                                pstep = mfp(sigma_t(irnew),b,bm_pos(irnew-1))
                            else
                                bm_loop_xpos: do i = irnew,(nreg-1),1 ! loop para conocer la region                                                                                           
                                    if(bm_pos(i) <= b .and. b < bm_pos(i+1)) then
                                        irnew = i+1
                                        pstep = mfp(sigma_t(irnew),b,bm_pos(irnew-1))
                                        exit
                                    endif
                                enddo bm_loop_xpos
                            endif
                        endif
                    endif

                    ir = irnew                

                    ! Transporte de la particula segun su direccion de movimiento
                    if(u .gt. 0.0) then
                        x = xbounds(irnew) + pstep*u
                    else if (u .lt. 0.0) then 
                        x = xbounds(irnew+1) + pstep*u
                    endif
                    
                    ! Verificacion de transporte de la particula fuera de la geometria
                    if(x > xbounds(nreg+1)) then
                        pdisc = .true.
                        ir = nreg+1
                    elseif(x < 0.0) then
                        pdisc = .true.
                        ir = 0
                    endif

                    ! Verificar descarte de la particula
                    if(pdisc .eqv. .true.) then
                        exit
                    endif
                    
                    if(ir == irnew) then
                        exit
                    else
                        ir = irnew 
                        exit    ! particula forzada a interactuar
                    endif

                enddo pmfp_loop    
                
                if(pdisc .eqv. .true.) then
                    score(ir) = score(ir) + wt
                    exit
                endif

                ! Seleccion del tipo de interaccion
                rnno = rng_set()
                if(rnno .le. sigma_a(ir)/sigma_t(ir)) then
                    score(ir) = score(ir) + wt
                    exit
                else
                    ! Dispersion de la particula
                    u = scatt(u)
                endif
                
            enddo particle_loop_mfp

        enddo ihist_loop_mfp

        write(*,'(A,F15.5)') 'Reflection : ', score(0)/real(nhist)
        write(*,'(A,F15.5)') 'Absorption : ', sum(score(1:nreg))/real(nhist)
        write(*,'(A,F15.5)') 'Transmission : ', score(nreg+1)/real(nhist)

        call cpu_time(end_time)
        write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time
        
    endif

end program shield_1d