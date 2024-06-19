MODULE my_calc_velo

    USE kinds

    IMPLICIT NONE

    CONTAINS

    ! subroutine to calculate average velocity in x and y direction
    ! for zcolu layers in z direction
    SUBROUTINE calc_velo_xy(path,input,zcolu,n,T,time,output)
        CHARACTER(LEN=*), INTENT(IN) :: path
        !!!! input: simulation box boundaries (6), interested area boundaries (6)
        REAL(REAL64), DIMENSION(12), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: zcolu  ! number of layers in z direction
        INTEGER, INTENT(IN) :: n, T ! number of particles, number of timesteps
        INTEGER, DIMENSION(:), INTENT(IN) :: time
        !!!! output: [timestep, zcolu, (vx,vy)]
        REAL(REAL64), DIMENSION(T,zcolu,2), INTENT(OUT) :: output
        CHARACTER(LEN=10) :: cTemp
        INTEGER :: i, j
        CHARACTER(LEN=200) :: filename
        REAL(REAL64) :: trash1, trash2, trash3, trash4, tempz, tempvx, tempvy, tempvz
        REAL(REAL64) :: xlen, ylen, zlen
        INTEGER, DIMENSION(zcolu) :: count_n ! count number of particles in a layer
        REAL(REAL64) :: dz
        INTEGER :: iz

        xlen = input(2) - input(1)
        ylen = input(4) - input(3)
        zlen = input(6) - input(5)

        dz = (input(12)-input(11))/zcolu

        output = 0.0d0

        DO i=1,T
            count_n = 0
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            OPEN(1, file=filename)
            DO j=1,9
                READ(1,*)
            END DO
            DO j=1,n
                READ(1,*) trash1, trash2, trash3, trash4, tempz, tempvx, tempvy, tempvz
                tempz = tempz * zlen + input(5)
                iz = INT(CEILING((tempz-input(11))/dz))
                IF (iz.GT.zcolu) CYCLE
                output(i,iz,1) = output(i,iz,1) + tempvx
                output(i,iz,2) = output(i,iz,2) + tempvy
                count_n(iz) = count_n(iz) + 1
            END DO
            CLOSE(1)
            DO j=1,zcolu
                IF (count_n(j).EQ.0) THEN
                    output(i,j,:) = 0
                ELSE
                    output(i,j,:) = output(i,j,:)/count_n(j)
                END IF
            END DO
        END DO
        

    END SUBROUTINE

END MODULE
