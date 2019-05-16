program mfp_sampling
    !
    ! Program that samples the free flight, i.e. distance to the next interaction. and estimates its 
    ! mean value and variance for different number of histories.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_mfp

    implicit none

    ! Transport parameters
    integer, parameter :: ncase = 9
    
    ! Caso a)
    integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]
    integer, dimension(ncase) :: nbatch = 1000/nperbatch   ! number of statistical batches
    
    ! Caso b)
    ! integer, parameter, dimension(ncase) :: nperbatch = 10
    ! integer, parameter, dimension(ncase) :: nbatch = [10, 20, 30, 40, 50, 100, 200, 400, 500]   ! number of statistical batches
        
    ! Caso c)
    ! integer, parameter, dimension(ncase) :: nbatch = 10   ! number of statistical batches
    ! integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]
    
    ! Caso d)
    ! integer, parameter, dimension(ncase) :: nbatch = 20   ! number of statistical batches
    ! integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]

    ! Geometry parameters
    real :: sigma = 2.0   ! total interaction cross section (cm-1)

    integer :: icase, ibatch, isample
    real :: path, score, score2, desv           ! auxiliary variables used for scoring
    real, allocatable :: step(:)    ! distance to next interaction (cm)
    real, allocatable :: step2(:)   ! distance to next interaction (squared)
    real, allocatable :: mean(:)    ! distance to next interaction (cm)
    real, allocatable :: var(:)   ! distance to next interaction (squared)

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Allocate the arrays that will hold the samples.
    allocate(step(size(nperbatch)))
    allocate(step2(size(nperbatch)))
    ! Allocate the arrays that will hold the mean and var.
    allocate(mean(size(nperbatch)))
    allocate(var(size(nperbatch)))

    ! do icase = 1,ncase
    !     write(unit=6, fmt='(I5)') nbatch(icase)       
    ! enddo

    step = 0.0
    step2 = 0.0
    mean = 0.0
    var = 0.0

    ! Start the sampling process. We calculate the mean and variance of the sample for 
    ! each nperbatch value given by the user.
    ncase_loop: do icase = 1,ncase
            
        ! Proceed with the sampling process. Start the batch loop
        batch_loop: do ibatch = 1,nbatch(icase)   

            ! Initialize scoring variable
            score = 0.0
            score2 = 0.0
            
            sampling_loop: do isample = 1,nperbatch(icase)
                ! Sampling process. Accumulate the sample and its squared value.
                path = mfp(sigma)
                score = score + path
                score2 = score2 + path**2
                
            enddo sampling_loop

            ! Accumulate results and proceed to next batch.
            step(icase) = step(icase) + score
            step2(icase) = step2(icase) + score2
            
        enddo batch_loop

        ! Statistical analysis.
        mean(icase) = step(icase)/(nbatch(icase)*nperbatch(icase))
        var(icase) = (step2(icase) - (nbatch(icase)*nperbatch(icase))*mean(icase)**2)/(((nbatch(icase)*nperbatch(icase))-1))

        write(unit=6, fmt='(I2, I5, I5, F10.5, F10.5)') icase, nperbatch(icase), nbatch(icase), mean(icase), var(icase)

    enddo ncase_loop

    ! Save results to file
    open(unit=1, file='mfp_sampling_case_a.txt')
    do icase = 1,ncase
        desv = sqrt(var(icase))
        write(unit=1, fmt='(I5, I5, F10.5, F10.5, F10.5)') nperbatch(icase), nbatch(icase), mean(icase), var(icase), desv        
    enddo
    close(1)

end program mfp_sampling