MODULE calc_rup_solid2d

    USE kinds

    IMPLICIT NONE

    CONTAINS
    !! ** Subroutine to calculate number density of liquid particles for solid 
    !! wall particles on a solid plate, for quasi-2D MD simulations.
    !! ** Use a wedge shape.
    !! ** Solid particles are uniformlly distributed on a plane, perpendicular
    !! to z-direction
    !! ** Since we have a quasi-2d problem, it is always periodic in y-direction.
    !! ** Since we are considering solid plate, it is always fix in z-direction.
    !! ** The boundary conditions in x-direction are given by the xlo and xhi of 
    !! solid plate.
    !! ** hwedge is the height of the wedge, hplate is the height of the plate,
    !! hcut is between hwedge and hplate.
    SUBROUTINE calc_dens(filename,input,Nliquid,Nsolid,output)
        ! filename for liquid
        CHARACTER(LEN=*), INTENT(IN) :: filename
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the wedge, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        hcut
        REAL(REAL64), DIMENSION(13), INTENT(IN) :: input
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolid
        ! output, {{Nsolid},{x,dens}}
        REAL(REAL64), DIMENSION(Nsolid), INTENT(OUT) :: output
        REAL(REAL64) :: slope, hwedge, hplate, xloplate,xhiplate
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posx, posy, posz
        ! position of solid particles, x-direction
        ! REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posSolid
        REAL(REAL64) :: xlen, ylen, zlen
        REAL(REAL64) :: tempx, tempy, tempz, trash1, trash2
        REAL(REAL64) :: tempxhi,tempxlo
        REAL(REAL64) :: volume, hcut
        INTEGER :: i, mlo, mhi, m
        ! spacing between solid particles
        REAL(REAL64) :: dx
        ! list to count number of liquid particles
        INTEGER, DIMENSION(:), ALLOCATABLE :: list

        slope = input(7)
        hwedge = input(8)
        hplate = input(9)
        xloplate = input(10)
        xhiplate = input(11)
        hcut = input(13)

        ! PRINT*, "OPEN FILE"
        OPEN(1, file=filename)

        ALLOCATE(posx(Nliquid))
        ALLOCATE(posy(Nliquid))
        ALLOCATE(posz(Nliquid))

        ! header has 9 lines
        DO i=1,9
            READ(1,*)
        END DO

        ! Rescale positions
        xlen = input(2) - input(1)
        ylen = input(4) - input(3)
        zlen = input(6) - input(5)

        ! read data
        DO i=1,Nliquid
            READ(1,*) trash1, trash2, tempx, tempy, tempz
            posx(i) = tempx * xlen + input(1)
            posy(i) = tempy * ylen + input(3)
            posz(i) = tempz * zlen + input(5)
        END DO

        CLOSE(1)

        output = 0.0d0

        ! now generate position of solid particles
        ! ALLOCATE(posSolid(Nsolid))
        dx = (xhiplate-xloplate)/REAL(Nsolid)
        ! DO i=1,Nsolid
        !     posSolid(i) = xloplate + dx/2.0d0 + dx*REAL(i-1)
        ! END DO

        ! allocate list
        ALLOCATE(list(Nsolid))
        list = 0
        output = 0.0d0
        ! For each liquid particle, allocated it to solid particles
        ! according to the wedge shape and calculate dens
        ! periodic
        IF (input(12).EQ.1) THEN
            DO i=1,Nliquid
                IF (posz(i).GT.hwedge) THEN
                    CYCLE
                ELSE
                    tempxlo = posx(i) - (posz(i)-hplate)/slope
                    tempxhi = posx(i) + (posz(i)-hplate)/slope
                    mlo = CEILING((tempxlo-xloplate+dx/2)/dx)
                    mhi = FLOOR((tempxhi-xloplate+dx/2)/dx)
                    ! touch left boundary
                    IF ((mlo.LT.1).AND.(mhi.LT.Nsolid)) THEN
                        DO m = 1,mhi
                            list(m) = list(m) + 1
                        END DO
                        DO m = Nsolid-(ABS(mlo)),Nsolid
                            list(m) = list(m) + 1
                        END DO
                    ! touch right boundary
                    ELSE IF ((mlo.GT.0).AND.(mhi.GT.Nsolid)) THEN
                        DO m = mlo, Nsolid
                            list(m) = list(m) + 1
                        END DO
                        DO m = 1, (mhi-Nsolid)
                            list(m) = list(m) + 1
                        END DO
                    ! touch both boundary, WRONG
                    ELSE IF ((mlo.LT.1).AND.(mhi.GT.Nsolid)) THEN
                        PRINT*, "error, reduce hwedge or increase slope"
                        STOP "1"
                    ! sits in middle
                    ELSE
                        DO m = mlo, mhi
                            list(m) = list(m) + 1
                        END DO
                    END IF
                END IF
            END DO
            ! now fill output
            volume = ((hwedge-hplate)/slope+(hcut-hplate)/slope)*(hwedge-hcut)*&
            (input(4)-input(3))
            output(:) = list/volume
        ! fix
        ELSE
            DO i = 1, Nliquid
                IF (posz(i).GT.hwedge) THEN
                    CYCLE
                ELSE
                    tempxlo = posx(i) - (posz(i)-hplate)/slope
                    tempxhi = posx(i) + (posz(i)-hplate)/slope
                    mlo = CEILING((tempxlo-xloplate+dx/2)/dx)
                    mhi = FLOOR((tempxhi-xloplate+dx/2)/dx)
                    ! touch left boundary
                    IF ((mlo.LT.1).AND.(mhi.LT.Nsolid)) THEN
                        DO m = 1,mhi
                            list(m) = list(m) + 1
                        END DO
                    ! touch right boundary
                    ELSE IF ((mlo.GT.0).AND.(mhi.GT.Nsolid)) THEN
                        DO m = mlo, Nsolid
                            list(m) = list(m) + 1
                        END DO
                    ! touch both boundary, WRONG
                    ELSE IF ((mlo.LT.1).AND.(mhi.GT.Nsolid)) THEN
                        PRINT*, "error, reduce hwedge or increase slope"
                        STOP "1"
                    ! sits in middle
                    ELSE
                        DO m = mlo, mhi
                            list(m) = list(m) + 1
                        END DO
                    END IF
                END IF
            END DO
            ! now fill output
            !!!!!! not finished
        END IF
        
        DEALLOCATE(posx)
        DEALLOCATE(posy)
        DEALLOCATE(posz)
        DEALLOCATE(list)
        ! DEALLOCATE(posSolid)

    END SUBROUTINE

    ! ** Subroutine to read in time, and then output dens
    ! ** Suitiable for mpi
    SUBROUTINE calc_dens_mpi(path,input,T,Nliquid,Nsolid,time,output)
        !!!! path is not full filename, but filenames without timestep
        CHARACTER(LEN=*), INTENT(IN) :: path
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the wedge, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        hcut
        REAL(REAL64), DIMENSION(13), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
                ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolid
        INTEGER, INTENT(IN), DIMENSION(:) :: time
        REAL(REAL64), DIMENSION(T,Nsolid), INTENT(OUT) :: output
        INTEGER :: i
        CHARACTER(LEN=10) :: cTemp
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: dens
        CHARACTER(LEN=300) :: filename

        output = 0
        DO i=1,t
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            ALLOCATE(dens(Nsolid))
            dens = 0
            CALL calc_dens(TRIM(ADJUSTL(filename)),input,Nliquid,Nsolid,dens)
            output(i,:) = dens
            DEALLOCATE(dens)
        END DO
    END SUBROUTINE


    !! ** Subroutine to calculate number density of liquid particles for solid 
    !! wall particles on a solid plate, for quasi-2D MD simulations.
    !! ** Use a wedge shape.
    !! ** Solid particles are uniformlly distributed on a plane, perpendicular
    !! to z-direction
    !! ** Since we have a quasi-2d problem, it is always periodic in y-direction.
    !! ** Since we are considering solid plate, it is always fix in z-direction.
    !! ** The boundary conditions in x-direction are given by the xlo and xhi of 
    !! solid plate.
    !! ** hwedge is the height of the wedge, hplate is the height of the plate,
    !! hcut is between hwedge and hplate.
    !! ** This is for onefile
    SUBROUTINE calc_dens_onefile(posx,posz,input,Nliquid,Nsolid,output)
        ! filename for liquid
        REAL(REAL64), DIMENSION(Nliquid), INTENT(IN) :: posx, posz
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the wedge, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        hcut
        REAL(REAL64), DIMENSION(13), INTENT(IN) :: input
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolid
        ! output, {{Nsolid},{x,dens}}
        REAL(REAL64), DIMENSION(Nsolid), INTENT(OUT) :: output
        REAL(REAL64) :: slope, hwedge, hplate, xloplate,xhiplate
        ! position of solid particles, x-direction
        ! REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posSolid
        REAL(REAL64) :: tempxhi,tempxlo
        REAL(REAL64) :: volume, hcut
        INTEGER :: i, mlo, mhi, m
        ! spacing between solid particles
        REAL(REAL64) :: dx
        ! list to count number of liquid particles
        INTEGER, DIMENSION(:), ALLOCATABLE :: list

        slope = input(7)
        hwedge = input(8)
        hplate = input(9)
        xloplate = input(10)
        xhiplate = input(11)
        hcut = input(13)

        output = 0.0d0

        ! now generate position of solid particles
        ! ALLOCATE(posSolid(Nsolid))
        dx = (xhiplate-xloplate)/REAL(Nsolid)
        ! DO i=1,Nsolid
        !     posSolid(i) = xloplate + dx/2.0d0 + dx*REAL(i-1)
        ! END DO

        ! allocate list
        ALLOCATE(list(Nsolid))
        list = 0
        output = 0.0d0
        ! For each liquid particle, allocated it to solid particles
        ! according to the wedge shape and calculate dens
        ! periodic
        IF (input(12).EQ.1) THEN
            DO i=1,Nliquid
                IF (posz(i).GT.hwedge) THEN
                    CYCLE
                ELSE
                    tempxlo = posx(i) - (posz(i)-hplate)/slope
                    tempxhi = posx(i) + (posz(i)-hplate)/slope
                    mlo = CEILING((tempxlo-xloplate+dx/2)/dx)
                    mhi = FLOOR((tempxhi-xloplate+dx/2)/dx)
                    ! touch left boundary
                    IF ((mlo.LT.1).AND.(mhi.LT.Nsolid)) THEN
                        DO m = 1,mhi
                            list(m) = list(m) + 1
                        END DO
                        DO m = Nsolid-(ABS(mlo)),Nsolid
                            list(m) = list(m) + 1
                        END DO
                    ! touch right boundary
                    ELSE IF ((mlo.GT.0).AND.(mhi.GT.Nsolid)) THEN
                        DO m = mlo, Nsolid
                            list(m) = list(m) + 1
                        END DO
                        DO m = 1, (mhi-Nsolid)
                            list(m) = list(m) + 1
                        END DO
                    ! touch both boundary, WRONG
                    ELSE IF ((mlo.LT.1).AND.(mhi.GT.Nsolid)) THEN
                        PRINT*, "error, reduce hwedge or increase slope"
                        STOP "1"
                    ! sits in middle
                    ELSE
                        DO m = mlo, mhi
                            list(m) = list(m) + 1
                        END DO
                    END IF
                END IF
            END DO
            ! now fill output
            volume = ((hwedge-hplate)/slope+(hcut-hplate)/slope)*(hwedge-hcut)*&
            (input(4)-input(3))
            output(:) = list/volume
        ! fix
        ELSE
            DO i = 1, Nliquid
                IF (posz(i).GT.hwedge) THEN
                    CYCLE
                ELSE
                    tempxlo = posx(i) - (posz(i)-hplate)/slope
                    tempxhi = posx(i) + (posz(i)-hplate)/slope
                    mlo = CEILING((tempxlo-xloplate+dx/2)/dx)
                    mhi = FLOOR((tempxhi-xloplate+dx/2)/dx)
                    ! touch left boundary
                    IF ((mlo.LT.1).AND.(mhi.LT.Nsolid)) THEN
                        DO m = 1,mhi
                            list(m) = list(m) + 1
                        END DO
                    ! touch right boundary
                    ELSE IF ((mlo.GT.0).AND.(mhi.GT.Nsolid)) THEN
                        DO m = mlo, Nsolid
                            list(m) = list(m) + 1
                        END DO
                    ! touch both boundary, WRONG
                    ELSE IF ((mlo.LT.1).AND.(mhi.GT.Nsolid)) THEN
                        PRINT*, "error, reduce hwedge or increase slope"
                        STOP "1"
                    ! sits in middle
                    ELSE
                        DO m = mlo, mhi
                            list(m) = list(m) + 1
                        END DO
                    END IF
                END IF
            END DO
            ! now fill output
            !!!!!! not finished
        END IF
        
        DEALLOCATE(list)
        ! DEALLOCATE(posSolid)

    END SUBROUTINE

    ! ** Subroutine to read in time, and then output dens
    ! ** Suitiable for mpi
    SUBROUTINE calc_dens_mpi_onefile(filename,input,T,Nliquid,Nsolid,time,output)
        !!!! filename
        CHARACTER(LEN=*), INTENT(IN) :: filename
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the wedge, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        hcut
        REAL(REAL64), DIMENSION(13), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
                ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolid
        !!!! time: init_timestep -- initial timestep of the simulation
        !!!!       begin_timestep -- first timestep to process
        !!!!       end_timestep -- final timestep to process
        !!!!       inter_timestep -- interval between timesteps
        INTEGER, INTENT(IN), DIMENSION(4) :: time
        REAL(REAL64), DIMENSION(T,Nsolid), INTENT(OUT) :: output
        INTEGER :: i, j
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: dens
        REAL(REAL64) :: xlen, ylen, zlen
        INTEGER :: skip_lines
        REAL(REAL64) :: tempx, tempy, tempz, trash1, trash2
        REAL(REAL64), DIMENSION(Nliquid) :: posx, posz

        ! Rescale positions
        xlen = input(2) - input(1)
        ylen = input(4) - input(3)
        zlen = input(6) - input(5)

        OPEN(1, file=filename)

        !! skip lines
        skip_lines = (time(2)-time(1))/time(4)
        DO i=1,skip_lines
            DO j=1,(9+Nliquid)
                READ(1,*)
            END DO
        END DO

        output = 0
        DO i=1,T
            ! skip headers
            DO j=1,9
                READ(1,*)
            END DO
            !! allocate position vectors
            posx = 0
            posz = 0
            DO j=1,Nliquid
                READ(1,*) trash1, trash2, tempx, tempy, tempz
                posx(j) = tempx * xlen + input(1)
                posz(j) = tempz * zlen + input(5)
            END DO
            ALLOCATE(dens(Nsolid))
            dens = 0
            CALL calc_dens_onefile(posx,posz,input,Nliquid,Nsolid,dens)
            output(i,:) = dens
            DEALLOCATE(dens)
        END DO

        CLOSE(1)
    END SUBROUTINE

    


END MODULE