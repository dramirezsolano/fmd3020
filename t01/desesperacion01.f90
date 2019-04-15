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
    integer, parameter :: nreg = 6                ! number of region in the 1D shield
    real, dimension(nreg) :: xthick = 1.0         ! thickness of the 1D shield (cm)
    real, dimension(nreg+1) :: xbounds = 0.0      ! boundaries of the shield regions (cm)

    ! Transport parameters
    real, parameter :: c = 0.5 ! Case 1
    ! real, parameter :: c = 0.7 ! Case 2
    ! real, parameter :: c = 0.9 ! Case 3
    real, parameter :: sigma_tin = 1.0
    real, dimension(nreg) :: bmin       ! Allocate the number of the intial mean-free-paths
    real, dimension(nreg) :: bm         ! Allocate the number of the intial mean-free-paths
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
    integer, parameter :: nhist = 100000

    ! Scoring variables
    real, dimension(0:nreg+1) :: score = 0.0    ! score(0) : reflection
                                                ! score(1:nreg) : absorption
                                                ! score(nreg+1) : transmission

    integer :: i, ihist 
    logical :: pdisc    ! flag to discard a particle
    integer :: irnew    ! index of new region 
    real :: pstep, b    ! distance to next interaction 
    real :: rnno 
    real :: start_time, end_time

    call cpu_time(start_time)

    ! Initialize RNG
    call rng_init(20180815)

    ! Print some info to console
    write(*,'(A, I15)') 'Number of histories COÑÑOOOOOOOOOOOO: ', nhist

    do i = 1,nreg
        if(i == 1) then
            sigma_t(i) = sigma_tin
            sigma_s(i) = c*sigma_t(i)
            sigma_a(i) = sigma_t(i) - sigma_s(i)
            bmin(i) = sigma_t(i)*xthick(i)
        else
            sigma_t(i) = sigma_t(i-1) - 0.1
            sigma_s(i) = c*sigma_t(i)
            sigma_a(i) = sigma_t(i) - sigma_s(i)
            bmin(i) = bmin(i-1) + sigma_t(i)*xthick(i)
        endif

        write(*,'(A, I5, A)') 'For region ', i, ':'
        write(*,'(A, F15.5)') 'Total interaction cross-section (cm-1): ', sigma_t(i)
        write(*,'(A, F15.5)') 'Absoption interaction cross-section (cm-1): ', sigma_a(i)
        write(*,'(A, F15.5)') 'Scattering interaction cross-section (cm-1): ', sigma_s(i)
        write(*,'(A, F15.5)') '1D shield thickness (cm): ', xthick(i)
        write(*,'(A, F15.5)') 'b (): ', bmin(i)

    enddo
   
    ! Initialize simulation geometry.
    write(*,'(A)') 'Region boundaries (cm):'
    write(*,'(F15.5)') xbounds(1)
    do i = 1,nreg
        xbounds(i+1) = xbounds(i) + xthick(i)
        write(*,'(F15.5)') xbounds(i+1)
    enddo
    

    ihist_loop: do ihist = 1,nhist
        
        ! Initialize history
        x = xin
        u = uin
        wt = wtin
        ir = irin
        bm = bmin
        
        ! Set flag used for particle discard.
        pdisc = .false.

        ! Enter transport process
        particle_loop: do
            
            pmfp_loop: do

                ! Distance to the next interaction.
                b = bk()

                if(u .gt. 0.0) then
                    if(ir < nreg) then
                        if(b <= bm(ir)) then 
                            ! irnew = ir
                            pstep = mfp(sigma_t(ir),b,bm(ir))
                        else
                            bm_loop_xpos: do i = ir,(nreg-1)                                                                                            
                                if(bm(i) < b .and. b <= bm(i+1)) then
                                    irnew = i+1
                                    pstep = mfp(sigma_t(irnew),b,bm(irnew))
                                    exit
                                else
                                    irnew = nreg+1
                                    ! Vacuum step
                                    pstep = 1.0E8
                                    ! The particle is leaving the geometry, discard it.
                                    pdisc = .true.
                                    exit
                                endif
                            enddo bm_loop_xpos
                        endif
                    else if(ir == nreg) then
                        if(b <= bm(nreg)) then
                            irnew = nreg  
                            pstep = mfp(sigma_t(irnew),b,bm(irnew))
                        else
                            ir = nreg+1
                            ! Vacuum step
                            pstep = 1.0E8
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                            exit
                        endif
                    endif
                ! else if (u .lt. 0.0) then
                !     if(ir > 1) then
                !         if(b <= bm(ir)) then 
                !             ! irnew = ir
                !             pstep = mfp(sigma_t(ir),b,bm(ir))
                !         else
                !             bm_loop_xneg: do i = ir,2                                                                                            
                !                 if(bm(i-1) < b .and. b <= bm(i)) then
                !                     irnew = i
                !                     pstep = mfp(sigma_t(irnew),b,bm(irnew))
                !                     exit
                !                 else
                !                     irnew = 0
                !                     ! Vacuum step
                !                     pstep = 1.0E8
                !                     ! The particle is leaving the geometry, discard it.
                !                     pdisc = .true.
                !                     exit
                !                 endif
                !             enddo bm_loop_xneg
                !         endif
                !     else if(ir == 1) then
                !         if(b <= bm(1)) then
                !             irnew = 1 
                !             pstep = mfp(sigma_t(irnew),b,bm(irnew))
                !         else
                !             ir = 0
                !             ! Vacuum step
                !             pstep = 1.0E8
                !             ! The particle is leaving the geometry, discard it.
                !             pdisc = .true.
                !             exit
                !         endif
                !     endif
                endif

                ! Save particle current region.
                ir = irnew

                ! Check if particle has been discarded.
                if(pdisc .eqv. .true.) then
                    exit
                endif

                ! Update position of particle.
                x = x + pstep*u
                
                ! Reached this point, if the particle has not changed region, it must interact.
                if(ir == irnew) then
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
                exit
            else
                ! The particle was scattered, get new direction to re-enter transport loop.
                u = scatt(u)
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