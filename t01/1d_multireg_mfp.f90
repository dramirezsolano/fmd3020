program shield_1d
    !
    ! Program that simulates the particle transport through a 1D shield using Monte Carlo methods.
    ! Multi-region version
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_mfp
    use mod_scatt

implicit none

    ! Geometry parameters
    integer, parameter :: nreg = 6           ! number of region in the 1D shield
    real, dimension(nreg) :: xthick = 1.0    ! thickness of the 1D shield (cm)
    real, dimension(nreg+1) :: xbounds_pos = 0.0      ! boundaries of the shield regions (cm)
    real, dimension(nreg+1) :: xbounds_neg = 0.0      ! boundaries of the shield regions (cm)

    ! Transport parameters
    ! real, parameter :: c = 0.5 ! Case 1
    ! real, parameter :: c = 0.7 ! Case 2
    real, parameter :: c = 0.9 ! Case 3
    real, parameter :: sigma_tin = 1.0
    real, dimension(0:nreg+1) :: bmin_pos, bm_pos       ! Allocate the number of the intial mean-free-paths
    real, dimension(0:nreg+1) :: bmin_neg, bm_neg         ! Allocate the number of the intial mean-free-paths
    real, dimension(nreg) :: sigma_t    ! Total interaction cross section (cm-1)
    real, dimension(nreg) :: sigma_s    ! Scattering cross section (cm-1)
    real, dimension(nreg) :: sigma_a    ! Absorption cross section (cm-1)

    ! Particle parameters
    real, parameter :: xin = 0.0        ! initial position (cm)
    real, parameter :: uin = 1.0        ! initial direction
    real, parameter :: wtin = 1.0       ! statistical weight
    integer, parameter :: irin = 1         ! initial region

    integer :: ir
    real :: x, u, wt

    ! Simulation parameters
    integer, parameter :: nhist = 500

    ! Scoring variables
    real, dimension(0:nreg+1) :: score = 0.0    ! score(0) : reflection
                                                ! score(1:nreg) : absorption
                                                ! score(nreg+1) : transmission

    integer :: i, ihist 
    logical :: pdisc    ! flag to discard a particle
    integer :: irnew    ! index of new region 
    real :: pstep, r, b    ! distance to next interaction 
    real :: rnno 
    real :: start_time, end_time

    call cpu_time(start_time)

    ! Initialize RNG
    call rng_init(20180815)

    ! Print some info to console
    write(*,'(A, I15)') 'Number of histories ', nhist

    bm_neg = 0.0
    bm_pos = 0.0

    sigma_t(1) = sigma_tin
    sigma_s(1) = c*sigma_t(1)
    sigma_a(1) = sigma_t(1) - sigma_s(1)
    do i = 2,nreg
        sigma_t(i) = sigma_t(i-1) - 0.1
        sigma_s(i) = c*sigma_t(i)
        sigma_a(i) = sigma_t(i) - sigma_s(i)
    enddo

    bmin_pos(1) = sigma_t(1)*xthick(1)
    bmin_neg(nreg) = sigma_t(nreg)*xthick(nreg)
    do i = 2,nreg
        bmin_pos(i) = bmin_pos(i-1) + sigma_t(i)*xthick(i)
        bmin_neg(nreg-i+1) = bmin_neg(nreg-i+2) + sigma_t(nreg-i+1)*xthick(nreg-i+1)
    enddo

    do i = 1,nreg
        write(*,'(A, I5, A)') 'For region ', i, ':'
        write(*,'(A, F15.5)') 'Total interaction cross-section (cm-1): ', sigma_t(i)
        write(*,'(A, F15.5)') 'Absoption interaction cross-section (cm-1): ', sigma_a(i)
        write(*,'(A, F15.5)') 'Scattering interaction cross-section (cm-1): ', sigma_s(i)
        write(*,'(A, F15.5)') '1D shield thickness (cm): ', xthick(i)
        write(*,'(A, F15.5)') 'bm_pos (): ', bmin_pos(i)
        write(*,'(A, F15.5)') 'bm_neg (): ', bmin_neg(i)
    enddo
   
    ! Initialize simulation geometry.
    write(*,'(A)') 'Region boundaries (cm):'
    ! write(*,'(F15.5,    F15.5)') xbounds_pos(1), xbounds_neg(1)
    do i = 1,nreg
        xbounds_pos(i+1) = xbounds_pos(i) + xthick(i)
        xbounds_neg(nreg-i+1) = xbounds_neg(nreg-i+2) + xthick(nreg-i+1)
        ! write(*,'(F15.5,    F15.5)') xbounds_pos(i+1), xbounds_neg(nreg-i+1)
    enddo
    write(*,'(F15.5,    F15.5)') xbounds_pos(1), xbounds_neg(1)
    do i = 1,nreg
        ! xbounds_pos(i+1) = xbounds_pos(i) + xthick(i)
        ! xbounds_neg(nreg-i+1) = xbounds_neg(nreg-i+2) + xthick(nreg-i+1)
        write(*,'(F15.5,    F15.5)') xbounds_pos(i+1), xbounds_neg(i+1)
    enddo
    

    ihist_loop: do ihist = 1,nhist
        write(*,'(A,    I5)') 'Historia', ihist
        ! Initialize history
        x = xin
        u = uin
        wt = wtin
        ir = irin
        
        ! Set flag used for particle discard.
        pdisc = .false.

        ! Enter transport process
        particle_loop: do
            
            pmfp_loop: do

                ! Distance to the next interaction.
                if(ir == 0 .or. ir == nreg+1) then
                    ! Vacuum step
                    pstep = 1.0E8
                else
                    b = bk()
                endif                
                
                ! b = bk()
                write(*,'(F15.5,    I5, F15.5, F15.5)') x, ir, u, b
                
                ! Save particle current region.
                irnew = ir

                if (u .lt. 0.0) then
                    bm_neg = bm_calc(u, nreg, ir, xbounds_pos, x, sigma_t)
                    do i = 1,nreg
                        ! write(*,'(A, F15.5)') 'bm_pos (): ', bm_pos(i)
                        write(*,'(A, F15.5)') 'bm_neg (): ', bm_neg(i)
                    enddo

                    if(irnew == 0) then
                        ! The particle is leaving the geometry, discard it.
                        pdisc = .true.
                        ! Vacuum step
                        pstep = 1.0E8
                    else
                        if(b > bm_neg(1)) then 
                            irnew = 0
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                            ! Vacuum step
                            pstep = 1.0E8
                        elseif(b <= bm_neg(irnew)) then
                            pstep = mfp(sigma_t(irnew),b,bm_neg(irnew+1))
                            ! r = mfp(sigma_t(irnew),b,bm_neg(irnew))
                        else
                            bm_loop_xneg: do i = irnew,2,-1
                                if(bm_neg(i-1) > b .and. b >= bm_neg(i)) then
                                    irnew = i-1
                                    pstep = mfp(sigma_t(irnew),b,bm_neg(irnew+1))
                                    ! write(*,'(A)') 'ENTROOOOO LT'
                                ! if(bm_neg(i-1) >= b .and. b > bm_neg(i)) then
                                !     irnew = i-1
                                !     r = mfp(sigma_t(irnew),b,bm_neg(irnew))
                                    exit
                                endif
                            enddo bm_loop_xneg
                        endif                                        
                    endif

                else if (u .gt. 0.0) then
                    bm_pos = bm_calc(u, nreg, ir, xbounds_pos, x, sigma_t)
                    do i = 1,nreg
                        write(*,'(A, F15.5)') 'bm_pos (): ', bm_pos(i)
                        ! write(*,'(A, F15.5)') 'bm_neg (): ', bm_neg(i)
                    enddo

                    if(irnew == nreg+1) then
                        ! The particle is leaving the geometry, discard it.
                        pdisc = .true.
                        ! Vacuum step
                        pstep = 1.0E8
                    else
                        if(b > bm_pos(nreg)) then 
                            irnew = nreg+1
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                            ! Vacuum step
                            pstep = 1.0E8
                        elseif(b <= bm_pos(irnew)) then
                            pstep = mfp(sigma_t(irnew),b,bm_pos(irnew-1))
                            ! r = mfp(sigma_t(irnew),b,bm_pos(irnew))
                        else
                            bm_loop_xpos: do i = irnew,(nreg-1),1                                                                                            
                                if(bm_pos(i) <= b .and. b < bm_pos(i+1)) then
                                    irnew = i+1
                                    pstep = mfp(sigma_t(irnew),b,bm_pos(irnew-1))
                                    ! write(*,'(A)') 'ENTROOOOO GT'
                                ! if(bm_pos(i) < b .and. b <= bm_pos(i+1)) then
                                !     irnew = i+1
                                !     r = mfp(sigma_t(irnew),b,bm_pos(irnew))
                                    exit
                                endif
                            enddo bm_loop_xpos
                        endif
                    endif
                endif

                ! Save particle current region.
                ir = irnew                
                write(*,'(A,F10.5)') 'pstep :', pstep

                ! Update position of particle.
                if(u .gt. 0.0) then
                    x = xbounds_pos(irnew) + pstep*u
                !     write(*,'(A, F15.5)') 'x_pos (): ', x
                else if (u .lt. 0.0) then 
                    x = xbounds_pos(irnew+1) + pstep*u
                !     write(*,'(A, F15.5)') 'x_neg (): ', x
                endif
                ! x = xbounds_pos(irnew) + pstep*u
                write(*,'(A, F15.5)') 'x (): ', x
                ! if(u .gt. 0.0) then
                !     pstep = xbounds_pos(irnew+1) - r*u
                !     x = x + pstep
                ! else if (u .lt. 0.0) then 
                !     pstep = xbounds_pos(irnew) + r*u
                !     x = x - pstep
                ! endif
                ! write(*,'(F15.5,    I5, F15.5)') x, ir, u

                if(x > xbounds_pos(nreg+1)) then
                    ! The particle is leaving the geometry, discard it.
                    pdisc = .true.
                    ir = nreg+1
                    write(*,'(A, F15.5)') 'TRANSMITIDO x :', x
                elseif(x < 0.0) then
                    ! The particle is leaving the geometry, discard it.
                    pdisc = .true.
                    ir = 0
                    write(*,'(A, F15.5)') 'REFLEJADO x :', x
                endif

                ! Check if particle has been discarded.
                if(pdisc .eqv. .true.) then
                    exit
                endif
                
                ! Reached this point, if the particle has not changed region, it must interact.
                if(ir == irnew) then
                    exit
                else
                    ! Update particle region index and continue transport process.
                    ir = irnew
                    exit
                endif

            enddo pmfp_loop    
            
            if(pdisc .eqv. .true.) then
                ! The particle has been discarded. Score it and end tracking.
                score(ir) = score(ir) + wt
                exit
            endif

            ! Determine interaction type.
            rnno = rng_set()
            if(rnno .le. sigma_a(ir)/sigma_t(ir)) then
                ! The particle was absorbed.
                score(ir) = score(ir) + wt
                write(*,'(A)') 'ABSORBIDO '
                exit
            else
                ! The particle was scattered, get new direction to re-enter transport loop.
                u = scatt(u)
                ! bm_pos = bm_calc(u, nreg, ir, xbounds_pos, x, sigma_t)
                ! bm_neg = bm_calc(u, nreg, ir, xbounds_neg, x, sigma_t)
            endif
            
        enddo particle_loop

    enddo ihist_loop

    ! Print results to console
    write(*,'(A,F15.5)') 'Reflection : ', score(0)/real(nhist)
    write(*,'(A,F15.5)') 'Absorption : ', sum(score(1:nreg))/real(nhist)
    write(*,'(A,F15.5)') 'Transmission : ', score(nreg+1)/real(nhist)

    ! Get elapsed time
    call cpu_time(end_time)
    write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time

end program shield_1d