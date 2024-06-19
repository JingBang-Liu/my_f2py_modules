MODULE calc_rup_solid3d

    USE kinds

    IMPLICIT NONE

    CONTAINS
    !! ** Subroutine to calculate number density of liquid particles for solid 
    !! wall particles on a solid plate, for 3D MD simulations, rectangular domain.
    !! ** Use a up-side-down cone.
    !! ** Solid particles are uniformlly distributed on a plane, perpendicular
    !! to z-direction
    !! ** Since we have a 3d problem, always periodic in x and y-direction.
    !! ** Since we are considering solid plate, it is always fix in z-direction.
    !! ** The boundary conditions in x-direction are given by the xlo and xhi of 
    !! solid plate.
    !! ** The boundary conditions in y-direction are given by the ylo and yhi of 
    !! solid plate.
    !! ** hcone is the height of the cone, hplate is the height of the plate,
    !! hcut is between hcone and hplate.
    SUBROUTINE calc_dens3d(filename,input,Nliquid,Nsolidx,Nsolidy,output)
        ! filename for liquid
        CHARACTER(LEN=*), INTENT(IN) :: filename
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the wedge, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        ylo of solid plate, yhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        boundary condition in y-direction (0-f, 1-p),
        !        hcut
        REAL(REAL64), DIMENSION(16), INTENT(IN) :: input
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolidx, Nsolidy
        ! output, {{Nsolid},{x,dens}}
        REAL(REAL64), DIMENSION(Nsolidx,Nsolidy), INTENT(OUT) :: output
        REAL(REAL64) :: slope, hcone, hplate, xloplate,xhiplate, yloplate, yhiplate
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posx, posy, posz
        ! position of solid particles, x-direction
        ! REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posSolidx, posSolidy
        REAL(REAL64) :: xlen, ylen, zlen
        REAL(REAL64) :: tempx, tempy, tempz, trash1, trash2
        REAL(REAL64) :: tempxhi,tempxlo,tempyhi,tempylo
        REAL(REAL64) :: volume, hcut, tempr, tempdist
        INTEGER :: i, mlo, mhi, m, nlo, nhi, n
        ! spacing between solid particles
        REAL(REAL64) :: dx, dy
        ! list to count number of liquid particles
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: list

        slope = input(7)
        hcone = input(8)
        hplate = input(9)
        xloplate = input(10)
        xhiplate = input(11)
        yloplate = input(12)
        yhiplate = input(13)
        hcut = input(16)

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
        ! ALLOCATE(posSolidx(Nsolidx))
        dx = (xhiplate-xloplate)/REAL(Nsolidx)
        ! DO i=1,Nsolidx
        !     posSolidx(i) = xloplate + dx/2.0d0 + dx*REAL(i-1)
        ! END DO
        ! ALLOCATE(posSolidy(Nsolidy))
        dy = (yhiplate-yloplate)/REAL(Nsolidy)
        ! DO i=1,Nsolidy
        !     posSolidy(i) = yloplate + dy/2.0d0 + dy*REAL(i-1)
        ! END DO

        ! allocate list
        ALLOCATE(list(Nsolidx,Nsolidy))
        list = 0
        output = 0.0d0
        ! For each liquid particle, allocated it to solid particles
        ! according to the cone shape and calculate dens
        ! periodic
        IF ((input(14).EQ.1).AND.(input(15).EQ.1)) THEN
            DO i=1,Nliquid
                IF (posz(i).GT.hcone) THEN
                    CYCLE
                ELSE
                    ! first check which x lines are crossed
                    tempr = (posz(i)-hplate)/slope
                    tempxlo = posx(i) - tempr
                    tempxhi = posx(i) + tempr
                    mlo = CEILING((tempxlo-xloplate+dx/2)/dx)
                    mhi = FLOOR((tempxhi-xloplate+dx/2)/dx)
                    ! touch left x boundary
                    IF ((mlo.LT.1).AND.(mhi.LT.Nsolidx)) THEN
                        DO m = 1,mhi
                            tempdist = SQRT(tempr**2-(posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "2"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                        DO m = Nsolidx-(ABS(mlo)),Nsolidx
                            tempdist = SQRT(tempr**2-(input(2)+posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "3"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                    ! touch right x boundary
                    ELSE IF ((mlo.GT.0).AND.(mhi.GT.Nsolidx)) THEN
                        DO m = mlo, Nsolidx
                            tempdist = SQRT(tempr**2-(posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "4"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                        DO m = 1, (mhi-Nsolidy)
                            tempdist = SQRT(tempr**2-(input(2)-posx(i)+(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "5"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                    ! touch both x boundary, WRONG
                    ELSE IF ((mlo.LT.1).AND.(mhi.GT.Nsolidx)) THEN
                        PRINT*, "error, reduce hwedge or increase slope"
                        STOP "1"
                    ! sits in middle
                    ELSE
                        DO m = mlo, mhi
                            tempdist = SQRT(tempr**2-(posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "6"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                    END IF
                END IF
            END DO
            ! now fill output
            volume = 3.1415926d0/3.0d0*(hcone**3-hplate**3)/slope**2
            output(:,:) = list/volume
        ! fix
        ELSE
            !!!!!! not finished
        END IF
        
        DEALLOCATE(posx)
        DEALLOCATE(posy)
        DEALLOCATE(posz)
        DEALLOCATE(list)
        ! DEALLOCATE(posSolidx)
        ! DEALLOCATE(posSolidy)

    END SUBROUTINE

    ! ** Subroutine to read in time, and then output dens
    ! ** Suitiable for mpi
    SUBROUTINE calc_dens_mpi3d(path,input,T,Nliquid,Nsolidx,Nsolidy,time,output)
        !!!! path is not full filename, but filenames without timestep
        CHARACTER(LEN=*), INTENT(IN) :: path
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the cone, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        ylo of solid plate, yhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        boundary condition in y-direction (0-f, 1-p)
        !        hcut
        REAL(REAL64), DIMENSION(16), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolidx, Nsolidy
        INTEGER, INTENT(IN), DIMENSION(:) :: time
        REAL(REAL64), DIMENSION(T,Nsolidx,Nsolidy), INTENT(OUT) :: output
        INTEGER :: i
        CHARACTER(LEN=10) :: cTemp
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens
        CHARACTER(LEN=300) :: filename

        output = 0
        DO i=1,t
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            ALLOCATE(dens(Nsolidx,Nsolidy))
            dens = 0
            CALL calc_dens3d(TRIM(ADJUSTL(filename)),input,Nliquid,Nsolidx,Nsolidy,dens)
            output(i,:,:) = dens
            DEALLOCATE(dens)
        END DO
    END SUBROUTINE

    !! ** Subroutine to calculate number density of liquid particles for solid 
    !! wall particles on a solid plate, for 3D MD simulations, rectangular domain.
    !! ** Use a up-side-down cone.
    !! ** Solid particles are uniformlly distributed on a plane, perpendicular
    !! to z-direction
    !! ** Since we have a 3d problem, always periodic in x and y-direction.
    !! ** Since we are considering solid plate, it is always fix in z-direction.
    !! ** The boundary conditions in x-direction are given by the xlo and xhi of 
    !! solid plate.
    !! ** The boundary conditions in y-direction are given by the ylo and yhi of 
    !! solid plate.
    !! ** hcone is the height of the cone, hplate is the height of the plate,
    !! hcut is between hcone and hplate.
    !! ** This is the onefile version
    SUBROUTINE calc_dens3d_onefile(posx,posy,posz,input,Nliquid,Nsolidx,Nsolidy,output)
        REAL(REAL64), DIMENSION(Nliquid), INTENT(IN) :: posx, posy, posz
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the wedge, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        ylo of solid plate, yhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        boundary condition in y-direction (0-f, 1-p),
        !        hcut
        REAL(REAL64), DIMENSION(16), INTENT(IN) :: input
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolidx, Nsolidy
        ! output, {{Nsolid},{x,dens}}
        REAL(REAL64), DIMENSION(Nsolidx,Nsolidy), INTENT(OUT) :: output
        REAL(REAL64) :: slope, hcone, hplate, xloplate,xhiplate, yloplate, yhiplate
        ! position of solid particles, x-direction
        ! REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posSolidx, posSolidy
        REAL(REAL64) :: tempxhi,tempxlo,tempyhi,tempylo
        REAL(REAL64) :: volume, hcut, tempr, tempdist
        INTEGER :: i, mlo, mhi, m, nlo, nhi, n
        ! spacing between solid particles
        REAL(REAL64) :: dx, dy
        ! list to count number of liquid particles
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: list

        slope = input(7)
        hcone = input(8)
        hplate = input(9)
        xloplate = input(10)
        xhiplate = input(11)
        yloplate = input(12)
        yhiplate = input(13)
        hcut = input(16)


        output = 0.0d0

        ! now generate position of solid particles
        ! ALLOCATE(posSolidx(Nsolidx))
        dx = (xhiplate-xloplate)/REAL(Nsolidx)
        ! DO i=1,Nsolidx
        !     posSolidx(i) = xloplate + dx/2.0d0 + dx*REAL(i-1)
        ! END DO
        ! ALLOCATE(posSolidy(Nsolidy))
        dy = (yhiplate-yloplate)/REAL(Nsolidy)
        ! DO i=1,Nsolidy
        !     posSolidy(i) = yloplate + dy/2.0d0 + dy*REAL(i-1)
        ! END DO

        ! allocate list
        ALLOCATE(list(Nsolidx,Nsolidy))
        list = 0
        output = 0.0d0
        ! For each liquid particle, allocated it to solid particles
        ! according to the cone shape and calculate dens
        ! periodic
        IF ((input(14).EQ.1).AND.(input(15).EQ.1)) THEN
            DO i=1,Nliquid
                IF (posz(i).GT.hcone) THEN
                    CYCLE
                ELSE
                    ! first check which x lines are crossed
                    tempr = (posz(i)-hplate)/slope
                    tempxlo = posx(i) - tempr
                    tempxhi = posx(i) + tempr
                    mlo = CEILING((tempxlo-xloplate+dx/2)/dx)
                    mhi = FLOOR((tempxhi-xloplate+dx/2)/dx)
                    ! touch left x boundary
                    IF ((mlo.LT.1).AND.(mhi.LT.Nsolidx)) THEN
                        DO m = 1,mhi
                            tempdist = SQRT(tempr**2-(posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "2"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                        DO m = Nsolidx-(ABS(mlo)),Nsolidx
                            tempdist = SQRT(tempr**2-(input(2)+posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "3"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                    ! touch right x boundary
                    ELSE IF ((mlo.GT.0).AND.(mhi.GT.Nsolidx)) THEN
                        DO m = mlo, Nsolidx
                            tempdist = SQRT(tempr**2-(posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "4"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                        DO m = 1, (mhi-Nsolidy)
                            tempdist = SQRT(tempr**2-(input(2)-posx(i)+(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "5"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                    ! touch both x boundary, WRONG
                    ELSE IF ((mlo.LT.1).AND.(mhi.GT.Nsolidx)) THEN
                        PRINT*, "error, reduce hwedge or increase slope"
                        STOP "1"
                    ! sits in middle
                    ELSE
                        DO m = mlo, mhi
                            tempdist = SQRT(tempr**2-(posx(i)-(m*dx-dx/2))**2)
                            tempylo = posy(i) - tempdist
                            tempyhi = posy(i) + tempdist
                            nlo = CEILING((tempylo-yloplate+dy/2)/dy)
                            nhi = FLOOR((tempyhi-yloplate+dy/2)/dy)
                            ! touch lower y boundary
                            IF ((nlo.LT.1).AND.(nhi.LT.Nsolidy)) THEN
                                DO n = 1,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = Nsolidy-(ABS(nlo)),Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch upper y boundary
                            ELSE IF ((nlo.GT.0).AND.(nhi.GT.Nsolidy)) THEN
                                DO n = nlo,Nsolidy
                                    list(m,n) = list(m,n) + 1
                                END DO
                                DO n = 1, (nhi-Nsolidy)
                                    list(m,n) = list(m,n) + 1
                                END DO
                            ! touch both y boundary, WRONG
                            ELSE IF ((nlo.LT.1).AND.(nhi.GT.Nsolidy)) THEN
                                PRINT*, "error, reduce hcone or increase slope"
                                STOP "6"
                            ! sits in the middle
                            ELSE
                                DO n=nlo,nhi
                                    list(m,n) = list(m,n) + 1
                                END DO
                            END IF
                        END DO
                    END IF
                END IF
            END DO
            ! now fill output
            volume = 3.1415926d0/3.0d0*(hcone**3-hplate**3)/slope**2
            output(:,:) = list/volume
        ! fix
        ELSE
            !!!!!! not finished
        END IF
        
        DEALLOCATE(list)
        ! DEALLOCATE(posSolidx)
        ! DEALLOCATE(posSolidy)

    END SUBROUTINE

    ! ** Subroutine to read in time, and then output dens
    ! ** Suitiable for mpi
    ! ** This is the onefile version
    SUBROUTINE calc_dens_mpi3d_onefile(filename,input,T,Nliquid,Nsolidx,Nsolidy,time,output)
        !!!! the filename
        CHARACTER(LEN=*), INTENT(IN) :: filename
        ! input: simulation boundaries (6),
        !        slope (Z/X), height of the cone, height of solid plate,
        !        xlo of solid plate, xhi of solid plate,
        !        ylo of solid plate, yhi of solid plate,
        !        boundary condition in x-direction (0-f, 1-p),
        !        boundary condition in y-direction (0-f, 1-p)
        !        hcut
        REAL(REAL64), DIMENSION(16), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
        ! number of liquid particles
        INTEGER, INTENT(IN) :: Nliquid
        ! number of solid particles
        INTEGER, INTENT(IN) :: Nsolidx, Nsolidy
        !!!! time: init_timestep -- initial timestep of the simulation
        !!!!       begin_timestep -- first timestep to process
        !!!!       end_timestep -- final timestep to process
        !!!!       inter_timestep -- interval between timesteps
        INTEGER, INTENT(IN), DIMENSION(4) :: time
        REAL(REAL64), DIMENSION(T,Nsolidx,Nsolidy), INTENT(OUT) :: output
        INTEGER :: i, j
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens
        INTEGER :: skip_lines
        REAL(REAL64) :: xlen, ylen, zlen
        REAL(REAL64) :: tempx, tempy, tempz, trash1, trash2
        REAL(REAL64), DIMENSION(Nliquid) :: posx, posy, posz

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
            posy = 0
            posz = 0
            DO j=1,Nliquid
                READ(1,*) trash1, trash2, tempx, tempy, tempz
                posx(j) = tempx * xlen + input(1)
                posy(j) = tempy * ylen + input(3)
                posz(j) = tempz * zlen + input(5)
            END DO
            ALLOCATE(dens(Nsolidx,Nsolidy))
            dens = 0
            CALL calc_dens3d_onefile(posx,posy,posz,input,Nliquid,Nsolidx,Nsolidy,dens)
            output(i,:,:) = dens
            DEALLOCATE(dens)
        END DO

        CLOSE(1)
    END SUBROUTINE


END MODULE