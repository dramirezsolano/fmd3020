program mfp_sampling
    !
    ! Programa que muestrea la longitud de camino que recorre un partícula en un medio
    ! empleando el método de vuelo libre. Además estima el valor medio (promedio) y la varianza 
    ! para diferentes combinaciones entre número de lotes estadísticos e historias por lote.
    !
    ! mod_rng.f90 y mod_mfp.f90 no sufrieron modificaciones. Ambos son de (Doerner, 2019)
    !
    ! Para compilar ejecute: gfortran mfp.f90 mod_rng.f90 mod_mfp.f90 -o mfp.exe
    ! Para correr ejecute: ./mfp.exe
    !
    ! Autores: Jessica Hernández
    !          Daniel Ramirez
    !
    ! Basado en (Doerner, 2019)
    !
    use mod_rng
    use mod_mfp

    implicit none

    ! Número de casos
    integer, parameter :: ncase = 9

    ! Diferentes combinaciones entre número de historias por lote, y cantidad de lotes estadísticos.
    ! Descomente la combinación de interés y comente la que no desea ejecutar.
    
    ! Combinación a)
    integer, dimension(ncase) :: nhist = 1000
    integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]
    integer, dimension(ncase) :: nbatch = 1000/nperbatch   ! number of statistical batches
    
    ! Combinación b)
    ! integer, parameter, dimension(ncase) :: nperbatch = 10    ! Número de historias por lote
    ! integer, parameter, dimension(ncase) :: nbatch = [10, 20, 30, 40, 50, 100, 200, 400, 500]   ! Número de lotes estadísticos
    ! integer, dimension(ncase) :: nhist = nbatch*nperbatch     ! Número de historias total
        
    ! Combinación c)
    ! integer, parameter, dimension(ncase) :: nbatch = 10 
    ! integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]
    ! integer, dimension(ncase) :: nhist = nbatch*nperbatch
    
    ! Combinación d)
    ! integer, parameter, dimension(ncase) :: nbatch = 20
    ! integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]
    ! integer, dimension(ncase) :: nhist = nbatch*nperbatch

    ! Parámetros geométricos
    real :: sigma = 2.0   ! Sección eficaz total de interacción (cm-1)

    integer :: icase, ibatch, isample
    real :: path, score, score2, desv, r, fom     ! Variables auxiliares
    real, allocatable :: step(:)    ! Distancia a la siguiente interacción (cm)
    real, allocatable :: step2(:)   ! Distancia a la siguiente interaccion (al cuadrado)
    real, allocatable :: mean(:)    ! Arreglo para almacenar el promedio
    real, allocatable :: var(:)     ! Arreglo para almacenar la varianza
    real, dimension(ncase) :: start_time
    real, dimension(ncase) :: end_time

    ! Inicializa PRNG
    call rng_init(20180815)

    ! Almacena las distancias, promedios y varianzas.
    allocate(step(size(nperbatch)))
    allocate(step2(size(nperbatch)))
    allocate(mean(size(nperbatch)))
    allocate(var(size(nperbatch)))

    step = 0.0
    step2 = 0.0
    mean = 0.0
    var = 0.0

    ! Inicia el proceso de muestreo. Se calcula el promedio y varianza para caso 
    ncase_loop: do icase = 1,ncase

        ! Tiempo de inicio para cada caso
        call cpu_time(start_time(icase))
            
        ! Inicia proceso de muestreo para cada lote estadístico
        batch_loop: do ibatch = 1,nbatch(icase)   

            ! Inicializa variables de conteo
            score = 0.0
            score2 = 0.0
            
            sampling_loop: do isample = 1,nperbatch(icase)
                ! Proceso de muestreo. Se almacena el valor y su cuadrado.
                path = mfp(sigma)
                score = score + path
                score2 = score2 + path**2
                
            enddo sampling_loop

            ! Valores acumulados por lote
            step(icase) = step(icase) + score
            step2(icase) = step2(icase) + score2
            
        enddo batch_loop

        ! Análisis estadístico
        mean(icase) = step(icase)/(nbatch(icase)*nperbatch(icase))
        var(icase) = (step2(icase) - ((nbatch(icase)*nperbatch(icase))*mean(icase)**2))/((nbatch(icase)*nperbatch(icase))-1)

        ! Variables a mostrar en pantalla
        write(unit=6, fmt='(I2, I5, I5, F10.5, F13.8)') icase, nperbatch(icase), nbatch(icase), mean(icase), var(icase)

        call cpu_time(end_time(icase))

    enddo ncase_loop

    ! Salvando los resultados en un archivo. Cambie el nombre según corresponda
    open(unit=1, file='mfp_sampling_case.txt')
    write(unit=1, fmt=*) 'nhist ', 'nperb ', 'batch ', 'mean ', '   var verd ', '  var prom', '    DS', '        R', &
            '         FOM', '        time'
    do icase = 1,ncase
        ! Cálculo de la desviación estándar, varianza relativa y FOM
        desv = sqrt(var(icase)/nperbatch(icase))
        r = sqrt(var(icase))/(mean(icase)*sqrt(nperbatch(icase)+0.0))
        fom = 1/((r**2)*(end_time(icase)-start_time(icase)))

        write(unit=1, fmt='(I5, I5, I5, F10.5, F10.5, F10.5, F10.5, F10.5, F13.2, F10.5)') &
             nhist(icase), nperbatch(icase), nbatch(icase), mean(icase), var(icase), var(icase)/nperbatch(icase), &
             desv, r, fom, (end_time(icase)-start_time(icase))         
    enddo

    close(1)

end program mfp_sampling