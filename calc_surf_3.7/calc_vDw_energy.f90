MODULE calc_vDw

    USE kinds

    IMPLICIT NONE

    CONTAINS

    !! This subroutine calculates the van der Waals energy between 
    !  liquid particles and solid particles using the formula of
    !  particle <-> infinetly thick solid substrate
    !! If used for non-periodic boundary then probably wrong
    !! This is the multi files version
    SUBROUTINE calc_vDw_liquid_solid(filename,input,Nliquid,output)
        ! filename for liquid
        CHARACTER(LEN=*), INTENT(IN) :: filename
        ! input: simulation boundaries (6),
        !        upper liquid limit, lower liquid limit,
        !        wall position
        REAL(REAL64), DIMENSION(9), INTENT(IN) :: input
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        REAL(REAL64), INTENT(OUT) :: output
        REAL(REAL64) :: tempx, tempy, tempz, trash1, trash2
        REAL(REAL64) :: zlen, posz, dist
        INTEGER :: i

        output = 0.0d0

        ! PRINT*, "OPEN FILE"
        OPEN(1, file=filename)

        ! header has 9 lines
        DO i=1,9
            READ(1,*)
        END DO

        ! Rescale positions
        ! xlen = input(2) - input(1)
        ! ylen = input(4) - input(3)
        zlen = input(6) - input(5)

        ! read data
        DO i=1,Nliquid
            READ(1,*) trash1, trash2, tempx, tempy, tempz
            posz = tempz * zlen + input(5)
            ! if the position of liquid particle is in the 
            ! consideration zone
            IF ((posz.GT.input(7)).AND.(posz.LT.input(8))) THEN
                dist = posz - input(9)
                output = output - 1/dist/dist/dist
            END IF
        END DO

        CLOSE(1)


    END SUBROUTINE

    !! Subroutine to read in bunch of times, and output energy
    !! Together with calc_vDw_liquid_solid
    SUBROUTINE calc_vDw_liquid_solid_mpi(path,input,T,Nliquid,time,output)
        !!!! path is not full filename, but filenames without timestep
        CHARACTER(LEN=*), INTENT(IN) :: path
        ! input: simulation boundaries (6),
        !        upper liquid limit, lower liquid limit,
        !        wall position
        REAL(REAL64), DIMENSION(9), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        INTEGER, INTENT(IN), DIMENSION(:) :: time
        REAL(REAL64), DIMENSION(T), INTENT(OUT) :: output
        INTEGER :: i
        CHARACTER(LEN=10) :: cTemp
        REAL(REAL64) :: energy
        CHARACTER(LEN=300) :: filename

        output = 0
        DO i=1,t
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            energy = 0
            CALL calc_vDw_liquid_solid(TRIM(ADJUSTL(filename)),input,Nliquid,energy)
            output(i) = energy
        END DO  

    END SUBROUTINE

    ! !! This subroutine calculates 
    ! SUBROUTINE calc_vDw_liquid_solid2()

    ! END SUBROUTINE





END MODULE