! This module is to calculate surface posiiton for 
! the configuration similar to liquid in a mug
! Using cylindrical coordinate
! The cylinder must be centered around z-axis
! Note: the volume calculation is an approximation 
! bases on the cutoff is much smaller than the radius
! of the cylinder, for example 3.5<<18
MODULE my_calc_surf_cylin1

    USE kinds

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE count_lines(filename,n_lines)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(OUT) :: n_lines
        INTEGER :: io

        PRINT*, "COUNT LINES STARTED"
        PRINT*, filename
        OPEN(1,file=filename, iostat=io, status='old')
        IF (io/=0) STOP 'Cannot open file!'

        n_lines = 0
        DO
            READ(1,*,iostat=io)
            IF (io/=0) EXIT
            n_lines = n_lines + 1
        END DO
        CLOSE(1)
        PRINT*, "FILE HAS ", n_lines, " LINES "
    END SUBROUTINE

    !!!!!!!! subroutine calculate density
    SUBROUTINE calc_dens(filename,input,n,output)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        !!!! everything is in reduced unit
        !!!! input: simulation box boundaries (6), interested area boundaries (6)
        !!!!        cutoff, radius of cylinder
        REAL(REAL64), DIMENSION(14), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: n ! number of atoms
        !!!! output: density of each particle [index, (x,y,z,dens)]
        REAL(REAL64), DIMENSION(n,4), INTENT(OUT) :: output
        INTEGER :: n_lines
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: pos_of_atoms
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posx, posy, posz
        REAL(REAL64) :: xlen, ylen, zlen ! length of dimension
        REAL(REAL64) :: dx, dy, dz, s
        REAL(REAL64) :: tempx, tempy, tempz, trash1, trash2
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: head
        INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: cell_xyz
        INTEGER, DIMENSION(3) :: n_cells
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: list
        REAL(REAL64), DIMENSION(3) :: sij
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: PH ! place holder to not repeat counting neighbours
        INTEGER :: ix, iy, iz
        INTEGER :: i, j, l0, l, r0, r, p0, p, t
        REAL(REAL64) :: volume
        REAL(REAL64) :: ccdist ! distance between the particle and the z-axis
        INTEGER :: cell_xyz_N


        !CALL count_lines(filename, n_lines)
        n_lines = n + 9

        ! PRINT*, "OPEN FILE"
        OPEN(1, file=filename)

        ALLOCATE(pos_of_atoms(n_lines-9,3))
        ALLOCATE(posx(n_lines-9))
        ALLOCATE(posy(n_lines-9))
        ALLOCATE(posz(n_lines-9))

        ! header has 9 lines
        DO i=1,9
            READ(1,*)
        END DO
        ! read data
        DO i=1,n_lines-9
            READ(1,*) trash1, trash2, tempx, tempy, tempz
            pos_of_atoms(i,1) = tempx
            pos_of_atoms(i,2) = tempy
            pos_of_atoms(i,3) = tempz
        END DO

        CLOSE(1)

        ! Rescale positions
        xlen = input(2) - input(1)
        ylen = input(4) - input(3)
        zlen = input(6) - input(5)

        posx(:) = pos_of_atoms(:,1) * xlen + input(1)
        posy(:) = pos_of_atoms(:,2) * ylen + input(3)
        posz(:) = pos_of_atoms(:,3) * zlen + input(5)

        n_cells(1) = FLOOR((input(8)-input(7))/input(13))
        n_cells(2) = FLOOR((input(10)-input(9))/input(13))
        n_cells(3) = FLOOR((input(12)-input(11))/input(13))

        dx = (input(8)-input(7))/n_cells(1)
        dy = (input(10)-input(9))/n_cells(2)
        dz = (input(12)-input(11))/n_cells(3)


        cell_xyz_N = FLOOR(n/4.0d0)
        ALLOCATE(head(n_cells(1),n_cells(2),n_cells(3)))
        ALLOCATE(cell_xyz(cell_xyz_N,n_cells(1),n_cells(2),n_cells(3)))
        ALLOCATE(list(n_lines-9))
        head = 0
        cell_xyz = 0
        list = 0
        ! A coarse allocation of atoms
        DO i=1,n_lines-9
            ix = FLOOR((posx(i)-input(7))/dx) + 1
            iy = FLOOR((posy(i)-input(9))/dy) + 1
            iz = FLOOR((posz(i)-input(11))/dz) + 1
            IF (ix.GT.n_cells(1)) THEN
                ix = n_cells(1)
            ELSE IF (ix.LT.1) THEN
                ix = 1
            END IF
            IF (iy.GT.n_cells(2)) THEN
                iy = n_cells(2)
            ELSE IF (iy.LT.1) THEN
                iy = 1
            END IF
            IF (iz.GT.n_cells(3)) THEN
                iz = n_cells(3)
            ELSE IF (iz.LT.1) THEN
                iz = 1
            END IF
            head(ix,iy,iz) = head(ix,iy,iz) + 1
            cell_xyz(head(ix,iy,iz),ix,iy,iz) = i
        END DO
        ALLOCATE(PH(n,n))
        PH = 0
        ! Count neighbours of atoms with in cutoff radius
        ! Different treatment for 0-fix boundary or 1-periodic boundary
        ! We only consider f f f boundaries
        DO i=1,n_lines-9
            ix = FLOOR((posx(i)-input(7))/dx) + 1
            iy = FLOOR((posy(i)-input(9))/dy) + 1
            iz = FLOOR((posz(i)-input(11))/dz) + 1
            IF (ix.GT.n_cells(1)) THEN
                ix = n_cells(1)
            ELSE IF (ix.LT.1) THEN
                ix = 1
            END IF
            IF (iy.GT.n_cells(2)) THEN
                iy = n_cells(2)
            ELSE IF (iy.LT.1) THEN
                iy = 1
            END IF
            IF (iz.GT.n_cells(3)) THEN
                iz = n_cells(3)
            ELSE IF (iz.LT.1) THEN
                iz = 1
            END IF
            DO l0=1,3
                l=ix+l0-2
                IF (l.EQ.0) THEN
                    CYCLE
                ELSEIF (l .EQ. (n_cells(1) + 1)) THEN
                    CYCLE
                END IF
                DO r0=1,3
                    r=iy+r0-2
                    IF (r .EQ. 0) THEN
                        CYCLE
                    ELSEIF (r .EQ. (n_cells(2) + 1)) THEN
                        CYCLE
                    END IF
                    DO p0=1,3
                        p=iz+p0-2
                        IF (p .EQ. 0) THEN
                            CYCLE
                        ELSE IF (p .EQ. (n_cells(3) + 1)) THEN
                            CYCLE
                        END IF
                        IF (head(l,r,p) .GT. 0) THEN
                            DO t=1,head(l,r,p)
                                j = cell_xyz(t,l,r,p)
                                IF ((j.GT.i) .AND. (PH(i,j).EQ.0)) THEN
                                    sij(1)=posx(i)-posx(j)
                                    sij(2)=posy(i)-posy(j)
                                    sij(3)=posz(i)-posz(j)
                                    PH(i,j) = 1
                                    PH(j,i) = 1
                                    ! IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                    !     sij(1) = (input(8)-input(7))-ABS(sij(1))
                                    ! END IF
                                    ! IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                    !     sij(2) = (input(10)-input(9))-ABS(sij(2))
                                    ! END IF
                                    ! IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                    !     sij(3) = (input(12)-input(11))-ABS(sij(3))
                                    ! END IF
                                    s=SQRT(sij(1)**2+sij(2)**2+sij(3)**2)
                                    IF (s.LT.input(13)) THEN
                                        list(i) = list(i) + 1
                                        list(j) = list(j) + 1
                                    END IF
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            END DO
        END DO
        
        
        !!!! fill output
        DO i=1,n_lines-9
            output(i,1) = posx(i)
            output(i,2) = posy(i)
            output(i,3) = posz(i)
            ccdist = SQRT(posx(i)**2+posy(i)**2)
            IF (((ccdist+input(13)).GT.input(14)).AND.(ccdist.LT.input(14))) THEN
                volume = REAL(1)/REAL(3)*3.1415926*(REAL(2)*input(13)**3+REAL(3)*input(13)**2*(input(14)-ccdist)&
                -(input(14)-ccdist)**3)
            ELSE IF (ccdist.GT.input(14)) THEN
                volume = 100000.0d0
                ! PRINT*, "OUT SIDE OF RADIUS", ccdist, input(14)
            ELSE
                volume = REAL(4)/REAL(3)*3.1415926*input(13)**3
            END IF
            ! volume = REAL(4)/REAL(3)*3.1415926*input(4)**3
            output(i,4) = REAL(list(i)+1)/volume
        END DO
    END SUBROUTINE

    !!!! Subroutine that calculates surface position
    !!!! This is only used for a compact object with surfaces
    !!!! This calculates for multiple timesteps
    SUBROUTINE calc_surf(path,input,T,xcolu,ycolu,time,output)
        !!!! path is not full filename, but filenames without timestep
        CHARACTER(LEN=*), INTENT(IN) :: path
        !!!! everything is in reduced unit
        !!!! input: simulation box boundaries (6), interested area boundaries (6)
        !!!!        cutoff, radius of cylinder
        !!!!        n, rate, (0-bot, 1-top)
        REAL(REAL64), DIMENSION(17), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
        INTEGER, INTENT(IN) :: xcolu ! number of x elements
        INTEGER, INTENT(IN) :: ycolu ! number of y elements
        INTEGER, INTENT(IN), DIMENSION(:) :: time
        !!!! output: surface position [timesteps, xcolu, ycolu]
        !!!!         top or bot determined by input(10) (0-bot, 1-top)
        REAL(REAL64), DIMENSION(T,xcolu,ycolu), INTENT(OUT) :: output
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens
        INTEGER :: n
        REAL(REAL64) :: rate, den_inter, z_surf_max, z_surf_min
        REAL(REAL64), DIMENSION(14) :: input2
        CHARACTER(LEN=300) :: filename 
        INTEGER :: i, j, k, s, p, m
        CHARACTER(LEN=10) :: cTemp
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens_liq, dens_vapour
        INTEGER :: indx, ix, iy
        REAL(REAL64) :: dx, dy
        INTEGER, DIMENSION(xcolu,ycolu) :: colu_xy
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: cell_xy
        ! REAL(REAL64) :: begin_time, end_time

        input2(1) = input(1) ! xlo
        input2(2) = input(2) ! xhi
        input2(3) = input(3) ! ylo
        input2(4) = input(4) ! yhi
        input2(5) = input(5) ! zlo
        input2(6) = input(6) ! zhi
        input2(7) = input(7) ! boundary xlo
        input2(8) = input(8) ! boundary xhi
        input2(9) = input(9) ! boundary ylo
        input2(10) = input(10) ! boundary yhi
        input2(11) = input(11) ! boundary zlo
        input2(12) = input(12) ! boundary zhi
        input2(13) = input(13) ! cutoff
        input2(14) = input(14) ! radius of cylinder
        n = INT(input(15))
        rate = input(16)
        output = 0

        DO i=1,T
            ! CALL CPU_TIME(begin_time)
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            ALLOCATE(dens(n,4))
            dens = 0
            CALL calc_dens(TRIM(ADJUSTL(filename)),input2,n,dens)
            den_inter = MAXVAL(dens(:,4))*rate
            indx = COUNT(dens(:,4).GT.den_inter)
            ALLOCATE(dens_liq(indx,4))
            ALLOCATE(dens_vapour(n-indx,4))
            k = 0
            s = 0
            DO j=1,n
                IF (dens(j,4).GT.den_inter) THEN
                    k = k + 1
                    dens_liq(k,:) = dens(j,:)
                ELSE
                    s = s + 1
                    dens_vapour(s,:) = dens(j,:)
                END IF
            END DO
            dx = (input(8)-input(7))/REAL(xcolu)
            dy = (input(10)-input(9))/REAL(ycolu)
            colu_xy = 0
            ALLOCATE(cell_xy(n,xcolu,ycolu))
            cell_xy = 0
            DO j=1,indx
                ix = INT(FLOOR((dens_liq(j,1)-input(7))/dx)) + 1
                iy = INT(FLOOR((dens_liq(j,2)-input(9))/dy)) + 1
                IF (ix.GT.xcolu) THEN
                    ix = xcolu
                ELSE IF (ix.LT.1) THEN
                    ix = 1
                END IF
                IF (iy.GT.ycolu) THEN
                    iy = ycolu
                ELSE IF (iy.LT.1) THEN
                    iy = 1
                END IF
                colu_xy(ix,iy) = colu_xy(ix,iy) + 1
                cell_xy(INT(colu_xy(ix,iy)),ix,iy) = j
            END DO
            k = 0
            IF (input(17).EQ.1) THEN ! output top_surf
                DO j=1,xcolu
                    DO s=1,ycolu
                        z_surf_max = -10000.0
                        z_surf_min = 10000.0
                        IF (colu_xy(j,s).GT.0.0) THEN
                            DO p=1,colu_xy(j,s)
                                m = cell_xy(p,j,s)
                                IF (dens_liq(m,3).GT.z_surf_max) THEN
                                    z_surf_max = dens_liq(m,3)
                                    output(i,j,s) = dens_liq(m,3)
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            ELSE IF (input(17).EQ.0) THEN ! output bot_surf
                DO j=1,xcolu
                    DO s=1,ycolu
                        z_surf_max = -10000.0
                        z_surf_min = 10000.0
                        IF (colu_xy(j,s).GT.0.0) THEN
                            DO p=1,colu_xy(j,s)
                                m = cell_xy(p,j,s)
                                IF (dens_liq(m,3).LT.z_surf_min) THEN
                                    z_surf_min = dens_liq(m,3)
                                    output(i,j,s) = dens_liq(m,3)
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            END IF
            DEALLOCATE(cell_xy)
            DEALLOCATE(dens_liq)
            DEALLOCATE(dens_vapour)


            DEALLOCATE(dens)
            ! CALL CPU_TIME(end_time)
            ! PRINT*, end_time-begin_time

        END DO
        
    END SUBROUTINE

    !!!! Instead of a rectangular mesh, we use a circular mesh
    SUBROUTINE calc_surf_cir(path,input,T,N_layer,maxblock,time,output)
        !!!! path is not full filename, but filenames without timestep
        CHARACTER(LEN=*), INTENT(IN) :: path
        !!!! everything is in reduced unit
        !!!! input: simulation box boundaries (6), interested area boundaries (6)
        !!!!        cutoff, radius of cylinder 
        !!!!        n, rate, (0-bot, 1-top)
        REAL(REAL64), DIMENSION(17), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
        INTEGER, INTENT(IN) :: N_layer ! number of layers in radial direction, including center
        INTEGER, INTENT(IN) :: maxblock ! maxium number of blocks in outer ring, 8*(N_layer-1)
        INTEGER, INTENT(IN), DIMENSION(:) :: time
        !!!! output: surface position [timesteps, xcolu, ycolu]
        !!!!         top or bot determined by input(10) (0-bot, 1-top)
        REAL(REAL64), DIMENSION(T,N_layer,maxblock), INTENT(OUT) :: output
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens
        INTEGER :: n
        REAL(REAL64) :: rate, den_inter, z_surf_max, z_surf_min
        REAL(REAL64), DIMENSION(14) :: input2
        CHARACTER(LEN=300) :: filename 
        INTEGER :: i, j, k, s, p, m
        CHARACTER(LEN=10) :: cTemp
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens_liq, dens_vapour
        REAL(REAL64) :: theta, dtheta
        INTEGER :: indx, ir, itheta
        INTEGER, DIMENSION(N_layer,maxblock) :: colu_rtheta
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: cell_rtheta
        REAL(REAL64) :: radius, temp_dr, temp_r
        ! REAL(REAL64) :: begin_time, end_time

        input2(1) = input(1) ! xlo
        input2(2) = input(2) ! xhi
        input2(3) = input(3) ! ylo
        input2(4) = input(4) ! yhi
        input2(5) = input(5) ! zlo
        input2(6) = input(6) ! zhi
        input2(7) = input(7) ! boundary xlo
        input2(8) = input(8) ! boundary xhi
        input2(9) = input(9) ! boundary ylo
        input2(10) = input(10) ! boundary yhi
        input2(11) = input(11) ! boundary zlo
        input2(12) = input(12) ! boundary zhi
        input2(13) = input(13) ! cutoff
        input2(14) = input(14) ! radius of cylinder
        radius = input(14)
        n = INT(input(15))
        rate = input(16)
        output = 0

        DO i=1,T
            IF (MOD(i,100).EQ.(0)) THEN
                PRINT*, time(i)
            ENDIF
            ! CALL CPU_TIME(begin_time)
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            ALLOCATE(dens(n,4))
            dens = 0
            CALL calc_dens(TRIM(ADJUSTL(filename)),input2,n,dens)
            den_inter = MAXVAL(dens(:,4))*rate
            indx = COUNT(dens(:,4).GT.den_inter)
            ALLOCATE(dens_liq(indx,4))
            ALLOCATE(dens_vapour(n-indx,4))
            k = 0
            s = 0
            DO j=1,n
                IF (dens(j,4).GT.den_inter) THEN
                    k = k + 1
                    dens_liq(k,:) = dens(j,:)
                ELSE
                    s = s + 1
                    dens_vapour(s,:) = dens(j,:)
                END IF
            END DO
            colu_rtheta = 0
            ALLOCATE(cell_rtheta(n,N_layer,maxblock))
            cell_rtheta = 0
            temp_dr = radius/(2*N_layer-1)
            DO j=1,indx
                temp_r = SQRT(dens_liq(j,1)**2+dens_liq(j,2)**2)
                ir = INT(FLOOR(temp_r/2.0d0/temp_dr-0.5d0)) + 2
                dtheta = 3.1415926/4.0d0/(ir-1.0d0)
                IF (ir==1) THEN
                    itheta = INT(1)
                    colu_rtheta(ir,itheta) = colu_rtheta(ir,itheta) + 1
                    cell_rtheta(INT(colu_rtheta(ir,itheta)),ir,itheta) = j
                ELSE IF ((ir.GT.1).AND.(ir.LE.N_layer)) THEN
                    IF (dens_liq(j,2).GE.0) THEN
                        theta = ACOS(dens_liq(j,1)/temp_r)
                        itheta = INT(FLOOR((theta+dtheta/2.0d0)/dtheta)) + 1
                    ELSE
                        theta = 2.0d0*3.1415926-ACOS(dens_liq(j,1)/temp_r)
                        itheta = INT(FLOOR((theta+dtheta/2.0d0)/dtheta)) + 1
                        IF (itheta.GT.(8.0d0*(ir-1.0d0))) THEN
                            itheta = 1
                        ENDIF
                    ENDIF
                    colu_rtheta(ir,itheta) = colu_rtheta(ir,itheta) + 1
                    cell_rtheta(INT(colu_rtheta(ir,itheta)),ir,itheta) = j
                ENDIF
            END DO
            k = 0
            IF (input(17).EQ.1) THEN ! output top_surf
                ! calculate the center cell
                z_surf_max = -10000.0
                z_surf_min = 10000.0
                IF (colu_rtheta(1,1).GT.0.0) THEN
                    DO p=1,colu_rtheta(1,1)
                        m = cell_rtheta(p,1,1)
                        IF (dens_liq(m,3).GT.z_surf_max) THEN
                            z_surf_max = dens_liq(m,3)
                            output(i,1,1) = dens_liq(m,3)
                        END IF
                    END DO
                END IF
                ! calculate outer layer cells
                DO j=2,N_layer
                    DO s=1,INT(8*(j-1))
                        z_surf_max = -10000.0
                        z_surf_min = 10000.0
                        IF (colu_rtheta(j,s).GT.0.0) THEN
                            DO p=1,colu_rtheta(j,s)
                                m = cell_rtheta(p,j,s)
                                IF (dens_liq(m,3).GT.z_surf_max) THEN
                                    z_surf_max = dens_liq(m,3)
                                    output(i,j,s) = dens_liq(m,3)
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            ELSE IF (input(17).EQ.0) THEN ! output bot_surf
                ! calculate the center cell
                z_surf_max = -10000.0
                z_surf_min = 10000.0
                IF (colu_rtheta(1,1).GT.0.0) THEN
                    DO p=1,colu_rtheta(1,1)
                        m = cell_rtheta(p,1,1)
                        IF (dens_liq(m,3).LT.z_surf_min) THEN
                            z_surf_min = dens_liq(m,3)
                            output(i,1,1) = dens_liq(m,3)
                        END IF
                    END DO
                END IF
                ! calculate outer layer cells
                DO j=2,N_layer
                    DO s=1,INT(8*(j-1))
                        z_surf_max = -10000.0
                        z_surf_min = 10000.0
                        IF (colu_rtheta(j,s).GT.0.0) THEN
                            DO p=1,colu_rtheta(j,s)
                                m = cell_rtheta(p,j,s)
                                IF (dens_liq(m,3).LT.z_surf_min) THEN
                                    z_surf_min = dens_liq(m,3)
                                    output(i,j,s) = dens_liq(m,3)
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            END IF
            DEALLOCATE(cell_rtheta)
            DEALLOCATE(dens_liq)
            DEALLOCATE(dens_vapour)


            DEALLOCATE(dens)
            ! CALL CPU_TIME(end_time)
            ! PRINT*, end_time-begin_time

        END DO



    END SUBROUTINE

END MODULE
