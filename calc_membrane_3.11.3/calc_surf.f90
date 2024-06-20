MODULE my_calc_surf

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
    !!!!!!!! suitable for onefile
    SUBROUTINE calc_dens_onefile(posx,posy,posz,input,n,output)
        !!!! position of atoms
        REAL(REAL64), DIMENSION(n), INTENT(IN) :: posx, posy, posz
        !!!! everything is in reduced unit
        !!!! input: simulation box boundaries (6), interested area boundaries (6)
        !!!!        cutoff, boundary (3) (0-f,1-p)
        REAL(REAL64), DIMENSION(16), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: n ! number of atoms
        !!!! output: density of each particle [index, (x,y,z,dens)]
        REAL(REAL64), DIMENSION(n,4), INTENT(OUT) :: output
        INTEGER :: n_lines
        ! REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: pos_of_atoms
        ! REAL(REAL64) :: xlen, ylen, zlen ! length of dimension
        REAL(REAL64) :: dx, dy, dz, s
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: head
        INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: cell_xyz
        INTEGER, DIMENSION(3) :: n_cells
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: list
        REAL(REAL64), DIMENSION(3) :: sij
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: PH ! place holder to not repeat counting neighbours
        INTEGER :: ix, iy, iz
        INTEGER :: i, j, l0, l, r0, r, p0, p, t
        REAL(REAL64) :: volume
        INTEGER :: cell_xyz_N


        !CALL count_lines(filename, n_lines)
        n_lines = n + 9



        ! check if unit is normalised or not
        ! IF ((MAXVAL(pos_of_atoms(:,1)) .LT. 1) .AND. (MAXVAL(pos_of_atoms(:,2)).LT.1) .AND. (MAXVAL(pos_of_atoms(:,3)).LT.1)) THEN
        !     posx(:) = pos_of_atoms(:,1) * input(1)
        !     posy(:) = pos_of_atoms(:,2) * input(2)
        !     posz(:) = pos_of_atoms(:,3) * input(3)
        ! ELSE
        !     posx(:) = pos_of_atoms(:,1)
        !     posy(:) = pos_of_atoms(:,2)
        !     posz(:) = pos_of_atoms(:,3)
        ! END IF

        ! Rescale positions
        ! xlen = input(2) - input(1)
        ! ylen = input(4) - input(3)
        ! zlen = input(6) - input(5)

        ! posx(:) = posx(:) * xlen + input(1)
        ! posy(:) = posy(:) * ylen + input(3)
        ! posz(:) = posz(:) * zlen + input(5)

        n_cells(1) = FLOOR((input(8)-input(7))/input(13))
        n_cells(2) = FLOOR((input(10)-input(9))/input(13))
        n_cells(3) = FLOOR((input(12)-input(11))/input(13))

        dx = (input(8)-input(7))/n_cells(1)
        dy = (input(10)-input(9))/n_cells(2)
        dz = (input(12)-input(11))/n_cells(3)

        ! ! check nothing is out of simulation box, only for periodic boundaries
        ! DO i=1,n_lines-9
        !     IF (input(14).EQ.1) THEN
        !         IF ((posx(i).LT.input(7)) .OR. (posx(i).GT.input(8))) THEN
        !             posx(i) = posx(i) - sign((input(8)-input(7)),(posx(i)-input(7)))
        !         END IF
        !     END IF
        !     IF (input(15).EQ.1) THEN
        !         IF ((posy(i).LT.input(9)) .OR. (posy(i).GT.input(10))) THEN
        !             posy(i) = posy(i) - sign((input(10)-input(9)),(posy(i)-input(9)))
        !         END IF
        !     END IF
        !     IF (input(16).EQ.1) THEN
        !         IF ((posz(i).LT.input(11)) .OR. (posz(i).GT.input(12))) THEN
        !             posz(i) = posz(i) - sign((input(12)-input(11)),(posz(i)-input(11)))
        !         END IF
        !     END IF
        ! END DO


        cell_xyz_N = FLOOR(n/n_cells(1)/n_cells(2)/1.0d0)
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
        PH = .TRUE.
        ! Count neighbours of atoms with in cutoff radius
        ! Different treatment for 0-fix boundary or 1-periodic boundary
        ! If p p p
        IF ((input(14).EQ.1).AND.(input(15).EQ.1).AND.(input(16).EQ.1)) THEN
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
                        l=n_cells(1)
                    ELSEIF (l .EQ. (n_cells(1) + 1)) THEN
                        l = 1
                    END IF
                    DO r0=1,3
                        r=iy+r0-2
                        IF (r .EQ. 0) THEN
                            r=n_cells(2)
                        ELSEIF (r .EQ. (n_cells(2) + 1)) THEN
                            r=1
                        END IF
                        DO p0=1,3
                            p=iz+p0-2
                            IF (p .EQ. 0) THEN
                                p=n_cells(3)
                            ELSE IF (p .EQ. (n_cells(3) + 1)) THEN
                                p=1
                            END IF
                            IF (head(l,r,p) .GT. 0) THEN
                                DO t=1,head(l,r,p)
                                    j = cell_xyz(t,l,r,p)
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If p p f
        ELSE IF ((input(14).EQ.1).AND.(input(15).EQ.1).AND.(input(16).EQ.0)) THEN
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
                        l=n_cells(1)
                    ELSEIF (l .EQ. (n_cells(1) + 1)) THEN
                        l = 1
                    END IF
                    DO r0=1,3
                        r=iy+r0-2
                        IF (r .EQ. 0) THEN
                            r=n_cells(2)
                        ELSEIF (r .EQ. (n_cells(2) + 1)) THEN
                            r=1
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
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If p f p
        ELSE IF ((input(14).EQ.1).AND.(input(15).EQ.0).AND.(input(16).EQ.1)) THEN
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
                        l=n_cells(1)
                    ELSEIF (l .EQ. (n_cells(1) + 1)) THEN
                        l = 1
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
                                p=n_cells(3)
                            ELSE IF (p .EQ. (n_cells(3) + 1)) THEN
                                p=1
                            END IF
                            IF (head(l,r,p) .GT. 0) THEN
                                DO t=1,head(l,r,p)
                                    j = cell_xyz(t,l,r,p)
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If p f f
        ELSE IF ((input(14).EQ.1).AND.(input(15).EQ.0).AND.(input(16).EQ.0)) THEN
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
                        l=n_cells(1)
                    ELSEIF (l .EQ. (n_cells(1) + 1)) THEN
                        l = 1
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
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If f p p
        ELSE IF ((input(14).EQ.0).AND.(input(15).EQ.1).AND.(input(16).EQ.1)) THEN
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
                            r=n_cells(2)
                        ELSEIF (r .EQ. (n_cells(2) + 1)) THEN
                            r=1
                        END IF
                        DO p0=1,3
                            p=iz+p0-2
                            IF (p .EQ. 0) THEN
                                p=n_cells(3)
                            ELSE IF (p .EQ. (n_cells(3) + 1)) THEN
                                p=1
                            END IF
                            IF (head(l,r,p) .GT. 0) THEN
                                DO t=1,head(l,r,p)
                                    j = cell_xyz(t,l,r,p)
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If f p f
        ELSE IF ((input(14).EQ.0).AND.(input(15).EQ.1).AND.(input(16).EQ.0)) THEN
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
                            r=n_cells(2)
                        ELSEIF (r .EQ. (n_cells(2) + 1)) THEN
                            r=1
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
                                    IF ((j.LT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If f f p
        ELSE IF ((input(14).EQ.0).AND.(input(15).EQ.0).AND.(input(16).EQ.1)) THEN
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
                                p=n_cells(3)
                            ELSE IF (p .EQ. (n_cells(3) + 1)) THEN
                                p=1
                            END IF
                            IF (head(l,r,p) .GT. 0) THEN
                                DO t=1,head(l,r,p)
                                    j = cell_xyz(t,l,r,p)
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        ! If f f f
        ELSE IF ((input(14).EQ.0).AND.(input(15).EQ.0).AND.(input(16).EQ.0)) THEN
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
                                    IF ((j.GT.i) .AND. PH(i,j)) THEN
                                        sij(1)=posx(i)-posx(j)
                                        sij(2)=posy(i)-posy(j)
                                        sij(3)=posz(i)-posz(j)
                                        PH(i,j) = .FALSE.
                                        PH(j,i) = .FALSE.
                                        IF (ABS(sij(1)).GT.0.5*(input(8)-input(7))) THEN
                                            sij(1) = (input(8)-input(7))-ABS(sij(1))
                                        END IF
                                        IF (ABS(sij(2)).GT.0.5*(input(10)-input(9))) THEN
                                            sij(2) = (input(10)-input(9))-ABS(sij(2))
                                        END IF
                                        IF (ABS(sij(3)).GT.0.5*(input(12)-input(11))) THEN
                                            sij(3) = (input(12)-input(11))-ABS(sij(3))
                                        END IF
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
        END IF
        
        !!!! fill output
        DO i=1,n_lines-9
            output(i,1) = posx(i)
            output(i,2) = posy(i)
            output(i,3) = posz(i)
            ! note, input(6) should be 1... change later
            IF ((input(14).EQ.0).AND.(input(15).EQ.1).AND.(input(16).EQ.0)) THEN 
                IF ((posx(i)-input(7)).LT.input(13)) THEN
                    volume = REAL(1)/REAL(3)*3.1415926*(REAL(2)*input(13)**3+REAL(3)*input(13)**2*(posx(i)-input(7))&
                    -(posx(i)-input(7))**3)
                ELSE IF ((posx(i)-input(7)).GT.(input(8)-input(13))) THEN
                    volume = REAL(1)/REAL(3)*3.1415926*(REAL(2)*input(13)**3+REAL(3)*input(13)**2*(input(8)-posx(i))&
                    -(input(8)-posx(i))**3)
                ELSE
                    volume = REAL(4)/REAL(3)*3.1415926*input(13)**3
                END IF
            ! if (x,y,z) are periodic, periodic, fix
            ELSE IF ((input(14).EQ.1).AND.(input(15).EQ.1).AND.(input(16).EQ.0)) THEN
                IF ((posz(i)-input(11)).LT.input(13)) THEN
                    volume = REAL(1)/REAL(3)*3.1415926*(REAL(2)*input(13)**3+REAL(3)*input(13)**2*(posz(i)-input(11))&
                    -(posz(i)-input(11))**3)
                ELSE IF ((posz(i)-input(11)).GT.(input(12)-input(13))) THEN
                    volume = REAL(1)/REAL(3)*3.1415926*(REAL(2)*input(13)**3+REAL(3)*input(13)**2*(input(12)-posz(i))&
                    -(input(12)-posz(i))**3)
                ELSE
                    volume = REAL(4)/REAL(3)*3.1415926*input(13)**3
                END IF
            ELSE
                volume = REAL(4)/REAL(3)*3.1415926*input(13)**3
            END IF
            ! volume = REAL(4)/REAL(3)*3.1415926*input(4)**3
            output(i,4) = REAL(list(i)+1)/volume
        END DO
        DEALLOCATE(head)
        DEALLOCATE(cell_xyz)
        DEALLOCATE(list)
        DEALLOCATE(PH)
    END SUBROUTINE

    !!!! Subroutine that calculates surface position
    !!!! This is only used for a compact object with surfaces
    !!!! This calculates for multiple timesteps, but reads from
    !!!!  one big file
    SUBROUTINE calc_surf_onefile(filename,input,T,n,dens_b,xcolu,ycolu,time,output)
        !!!! filename of the 
        CHARACTER(LEN=*), INTENT(IN) :: filename
        !!!! everything is in reduced unit
        !!!! input: cutoff, boundary (3) (0-f,1-p),
        !!!!        meaningless, rate, (0-bot, 1-top)
        REAL(REAL64), DIMENSION(7), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: T ! number of timesteps
        INTEGER, INTENT(IN) :: n ! number of atoms
        REAL(REAL64), INTENT(IN) :: dens_b ! bulk density
        INTEGER, INTENT(IN) :: xcolu ! number of x elements
        INTEGER, INTENT(IN) :: ycolu ! number of y elements
        !!!! time: init_timestep -- initial timestep of the simulation
        !!!!       begin_timestep -- first timestep to process
        !!!!       end_timestep -- final timestep to process
        !!!!       inter_timestep -- interval between timesteps
        INTEGER, INTENT(IN), DIMENSION(4) :: time
        !!!! output: surface position [timesteps, xcolu, ycolu]
        !!!!         top or bot determined by input(10) (0-bot, 1-top)
        REAL(REAL64), DIMENSION(T,xcolu,ycolu), INTENT(OUT) :: output
        REAL(REAL64) :: xlo, xhi, ylo, yhi, zlo, zhi
        REAL(REAL64) :: xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens
        ! greater than rate, then liquid
        REAL(REAL64) :: rate, den_inter, z_surf_max, z_surf_min
        REAL(REAL64), DIMENSION(16) :: input2
        INTEGER :: i, j, k, s, p, m
        ! CHARACTER(LEN=10) :: cTemp
        INTEGER, DIMENSION(:), ALLOCATABLE :: dens_liq_index
        INTEGER :: indx, ix, iy
        REAL(REAL64) :: dx, dy
        INTEGER, DIMENSION(xcolu,ycolu) :: colu_xy
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: cell_xy
        REAL(REAL64), DIMENSION(n):: posx, posy, posz
        REAL(REAL64) :: trash1, trash2, tempx, tempy, tempz
        REAL(REAL64) :: xlen, ylen, zlen
        INTEGER :: skip_lines
        ! REAL(REAL64) :: begin_time, end_time

        input2(13) = input(1) ! cutoff
        input2(14) = input(2) ! boundary 0-f 1-p
        input2(15) = input(3) ! boundary 
        input2(16) = input(4) ! boundary
        rate = input(6)
        output = 0





        OPEN(1, file=filename)

        !! skip lines
        skip_lines = (time(2)-time(1))/time(4)
        DO i=1,skip_lines
            DO j=1,(9+n)
                READ(1,*)
            END DO
        END DO
        

        DO i=1,T
            ! skip first part of headers
            DO j=1,5
                READ(1,*)
            END DO
            ! read simulation box dimensions
            READ(1,*) xlo, xhi
            READ(1,*) ylo, yhi
            READ(1,*) zlo, zhi
            ! calculation area is the same as simulation box
            xlo_bound = xlo
            xhi_bound = xhi
            ylo_bound = ylo
            yhi_bound = yhi
            zlo_bound = zlo
            zhi_bound = zhi
            ! box length
            xlen = xhi - xlo
            ylen = yhi - ylo
            zlen = zhi - zlo
            ! change input2
            input2(1) = xlo ! xlo
            input2(2) = xhi ! xhi
            input2(3) = ylo ! ylo
            input2(4) = yhi ! yhi
            input2(5) = zlo ! zlo
            input2(6) = zhi ! zhi
            input2(7) = xlo_bound ! boundary xlo
            input2(8) = xhi_bound ! boundary xhi
            input2(9) = ylo_bound ! boundary ylo
            input2(10) = yhi_bound ! boundary yhi
            input2(11) = zlo_bound ! boundary zlo
            input2(12) = zhi_bound ! boundary zhi
            !! allocate position vectors
            posx = 0
            posy = 0
            posz = 0
            DO j=1,n
                READ(1,*) trash1, trash2, tempx, tempy, tempz
                posx(j) = tempx * xlen + xlo
                posy(j) = tempy * ylen + ylo
                posz(j) = tempz * zlen + zlo
            END DO

            ! check nothing is out of simulation box, only for periodic boundaries
            DO j=1,n
                IF (input(2).EQ.1) THEN
                    IF ((posx(j).LT.xlo_bound) .OR. (posx(j).GT.xhi_bound)) THEN
                        posx(j) = posx(j) - sign((xhi_bound-xlo_bound),(posx(j)-xlo_bound))
                    END IF
                END IF
                IF (input(3).EQ.1) THEN
                    IF ((posy(j).LT.ylo_bound) .OR. (posy(j).GT.yhi_bound)) THEN
                        posy(j) = posy(j) - sign((yhi_bound-ylo_bound),(posy(j)-ylo_bound))
                    END IF
                END IF
                IF (input(4).EQ.1) THEN
                    IF ((posz(j).LT.zlo_bound) .OR. (posz(j).GT.zhi_bound)) THEN
                        posz(j) = posz(j) - sign((zhi_bound-zlo_bound),(posz(j)-zlo_bound))
                    END IF
                END IF
            END DO

            ! CALL CPU_TIME(begin_time)
            ALLOCATE(dens(n,4))
            dens = 0.0d0
            CALL calc_dens_onefile(posx,posy,posz,input2,n,dens)
            !den_inter = MAXVAL(dens(:,4))*rate
            ! den_inter = SUM(dens(:,4))/MAX(1,SIZE(dens(:,4)))*rate
            den_inter = REAL(dens_b*rate)
            indx = COUNT(dens(:,4).GT.den_inter)
            ALLOCATE(dens_liq_index(indx))
            dens_liq_index = 0
            ! ALLOCATE(dens_vapour(n-indx,4))
            k = 0
            ! s = 0
            DO j=1,n
                IF (dens(j,4).GT.den_inter) THEN
                    k = k + 1
                    dens_liq_index(k) = j
                ! ELSE
                    ! s = s + 1
                    ! dens_vapour(s,:) = dens(j,:)
                END IF
            END DO
            dx = (xhi_bound-xlo_bound)/REAL(xcolu)
            dy = (yhi_bound-ylo_bound)/REAL(ycolu)
            colu_xy = 0
            ALLOCATE(cell_xy(INT(n/xcolu/ycolu*4.0d0),xcolu,ycolu))
            cell_xy = 0
            DO j=1,indx
                ix = INT(FLOOR((dens(dens_liq_index(j),1)-xlo_bound)/dx)) + 1
                iy = INT(FLOOR((dens(dens_liq_index(j),2)-ylo_bound)/dy)) + 1
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
            IF (input(7).EQ.1) THEN ! output top_surf
                DO j=1,xcolu
                    DO s=1,ycolu
                        z_surf_max = -10000.0
                        z_surf_min = 10000.0
                        IF (colu_xy(j,s).GT.0.0) THEN
                            DO p=1,colu_xy(j,s)
                                m = cell_xy(p,j,s)
                                IF (dens(dens_liq_index(m),3).GT.z_surf_max) THEN
                                    z_surf_max = dens(dens_liq_index(m),3)
                                    output(i,j,s) = dens(dens_liq_index(m),3)
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            ELSE IF (input(7).EQ.0) THEN ! output bot_surf
                DO j=1,xcolu
                    DO s=1,ycolu
                        z_surf_max = -10000.0
                        z_surf_min = 10000.0
                        IF (colu_xy(j,s).GT.0.0) THEN
                            DO p=1,colu_xy(j,s)
                                m = cell_xy(p,j,s)
                                IF (dens(dens_liq_index(m),3).LT.z_surf_min) THEN
                                    z_surf_min = dens(dens_liq_index(m),3)
                                    output(i,j,s) = dens(dens_liq_index(m),3)
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            END IF
            DEALLOCATE(dens)
            DEALLOCATE(cell_xy)
            DEALLOCATE(dens_liq_index)
            ! DEALLOCATE(dens_vapour)

            ! CALL CPU_TIME(end_time)
            ! PRINT*, end_time-begin_time

        END DO


        CLOSE(1)
        
    END SUBROUTINE

    !!!! Subroutine output liquid and vapour particle positions
    SUBROUTINE vapour_liquid(path,input,n,T,time,output)
        CHARACTER(LEN=*), INTENT(IN) :: path
        !!!! everything is in reduced unit
        !!!! input: simulation box boundaries (6), interested area boundaries (6)
        !!!!        cutoff, boundary (3) (0-f,1-p), rate
        REAL(REAL64), DIMENSION(17), INTENT(IN) :: input
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: T ! number of timesteps
        INTEGER, INTENT(IN), DIMENSION(:) :: time
        !!!! output: surface position [timesteps, N, (type,x,y,z,dens)]
        !!!!         type 1-vapour 2-liquid
        REAL(REAL64), DIMENSION(T,n,5), INTENT(OUT) :: output
        REAL(REAL64), DIMENSION(:,:), ALLOCATABLE :: dens
        REAL(REAL64) :: rate, den_inter
        REAL(REAL64), DIMENSION(16) :: input2
        CHARACTER(LEN=300) :: filename 
        INTEGER :: i, j
        CHARACTER(LEN=10) :: cTemp

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
        input2(14) = input(14) ! boundary 0-f 1-p
        input2(15) = input(15) ! boundary 
        input2(16) = input(16) ! boundary
        rate = input(17)
        output = 0

        DO i=1,T
            ! CALL CPU_TIME(begin_time)
            WRITE(cTemp,'(i10)') time(i)
            filename = path//TRIM(ADJUSTL(cTemp))//'.dat'
            ALLOCATE(dens(n,4))
            dens = 0
            CALL calc_dens(TRIM(ADJUSTL(filename)),input2,n,dens)
            den_inter = MAXVAL(dens(:,4))*rate
            DO j=1,n
                IF (dens(j,4).GT.den_inter) THEN
                    output(i,j,1) = 2
                    output(i,j,2) = dens(j,1)
                    output(i,j,3) = dens(j,2)
                    output(i,j,4) = dens(j,3)
                    output(i,j,5) = dens(j,4)
                ELSE
                    output(i,j,1) = 1
                    output(i,j,2) = dens(j,1)
                    output(i,j,3) = dens(j,2)
                    output(i,j,4) = dens(j,3)
                    output(i,j,5) = dens(j,4)
                END IF
            END DO
            DEALLOCATE(dens)
            ! CALL CPU_TIME(end_time)
            ! PRINT*, end_time-begin_time
        END DO
        
    END SUBROUTINE

END MODULE
