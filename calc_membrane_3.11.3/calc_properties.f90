MODULE properties

    USE kinds

    IMPLICIT NONE

    CONTAINS


    !!!! This subroutine calculates the orientational order parameter
    !!!! ASSUMING PERIODIC BOUNDARY CONDITION IN ALL THREE DIRECTIONS
    SUBROUTINE order_parameter(filename,nhead,ntail,T,nmol,time,output)
        !! for meaning of variables, see calc_surf.calc_midplane_tail_onefile
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: nhead, ntail, T, nmol
        !!!! time: init_timestep -- initial timestep of the simulation
        !!!!       begin_timestep -- first timestep to process
        !!!!       end_timestep -- final timestep to process
        !!!!       inter_timestep -- interval between timesteps
        INTEGER, INTENT(IN), DIMENSION(4) :: time
        REAL(REAL64), DIMENSION(T), INTENT(OUT) :: output
        REAL(REAL64) :: xlo, xhi, ylo, yhi, zlo, zhi
        REAL(REAL64) :: xlen, ylen, zlen
        REAL(REAL64) :: xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound
        REAL(REAL64), DIMENSION(nmol,3) :: lipid_axis !! unit vector along the axis of lipids
        REAL(REAL64), DIMENSION(nhead,3) :: head_pos !! position of head beads
        REAL(REAL64), DIMENSION(ntail,3) :: tail_pos !! position of tail beads
        REAL(REAL64) :: diffx, diffy, diffz !! differnce between first head and last tail bead in x y z
        REAL(REAL64) :: temp_len, temp_prod
        REAL(REAL64), DIMENSION(nmol) :: lipid_order !! order of each lipid
        INTEGER :: i, j, k, n, skip_lines
        REAL(REAL64) :: trash1, trash2, trash3, trash4, trash5, trash6, tempx, tempy, tempz

        OPEN(1, file=filename)

        n = nmol*(nhead+ntail)

        !! skip lines
        skip_lines = (time(2)-time(1))/time(4)
        DO i=1,skip_lines
            DO j=1,(9+n)
                READ(1,*)
            END DO
        END DO

        !! actual calculation
        DO i=1,T
            ! skip first part of headers
            DO j=1,5
                READ(1,*)
            END DO
            ! read simulation box dimensions
            READ(1,*) xlo, xhi
            READ(1,*) ylo, yhi
            READ(1,*) zlo, zhi
            ! skip second part of headers
            READ(1,*)
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
            !! read beads positions and calculate axis vector
            DO j=1,nmol
                DO k=1,nhead
                    READ(1,*) trash1, trash2, trash3, tempx, tempy, tempz, trash5,&
                    & trash6, trash4
                    head_pos(k,1) = tempx * xlen + xlo
                    head_pos(k,2) = tempy * ylen + ylo
                    head_pos(k,3) = tempz * zlen + zlo
                END DO
                DO k=1,ntail
                    READ(1,*) trash1, trash2, trash3, tempx, tempy, tempz, trash5,&
                    & trash6, trash4
                    tail_pos(k,1) = tempx * xlen + xlo
                    tail_pos(k,2) = tempy * ylen + ylo
                    tail_pos(k,3) = tempz * zlen + zlo
                END DO
                diffx = head_pos(1,1) - tail_pos(ntail,1)
                diffy = head_pos(1,2) - tail_pos(ntail,2)
                diffz = head_pos(1,3) - tail_pos(ntail,3)
                IF (ABS(diffx).GT.(xlen/2)) THEN
                    lipid_axis(j,1) = (xlen - ABS(diffx))*SIGN(1.0d0,-diffx)
                ELSE
                    lipid_axis(j,1) = diffx
                END IF
                IF (ABS(diffy).GT.(ylen/2)) THEN
                    lipid_axis(j,2) = (ylen - ABS(diffy))*SIGN(1.0d0,-diffy)
                ELSE
                    lipid_axis(j,2) = diffy
                END IF
                IF (ABS(diffz).GT.(zlen/2)) THEN
                    lipid_axis(j,3) = (zlen - ABS(diffz))*SIGN(1.0d0,-diffz)
                ELSE
                    lipid_axis(j,3) = diffz
                END IF
                temp_len = SQRT(lipid_axis(j,1)**2+lipid_axis(j,2)**2+lipid_axis(j,3)**2)
                lipid_axis(j,:) = lipid_axis(j,:)/temp_len
                temp_prod = lipid_axis(j,3)
                lipid_order(j) = 3.0d0*(temp_prod)**2 - 1
            END DO
            output(i) = 0.5d0*SUM(lipid_order)/nmol
        END DO

        CLOSE(1)

    END SUBROUTINE

    !!!! This subroutine calculates the density profile across a membrane
    !!!! It actually calculates the 
    SUBROUTINE density_profile(filename,nhead,ntail,thickness,boundary,&
        &T,nmol,xcolu,ycolu,time,deltaz,z1,Nz,dim2,dim3,output)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: nhead ! number of head beads
        INTEGER, INTENT(IN) :: ntail ! number of tail beads
        REAL(REAL64), INTENT(IN) :: thickness ! initial guess of thickness
        INTEGER, DIMENSION(2), INTENT(IN) :: boundary ! 1: periodic, 0: fixed
        INTEGER, INTENT(IN) :: T ! number of timesteps
        INTEGER, INTENT(IN) :: nmol ! number of lipids
        INTEGER, INTENT(IN) :: xcolu ! number of x elements
        INTEGER, INTENT(IN) :: ycolu ! number of y elements
        !!!! time: init_timestep -- initial timestep of the simulation
        !!!!       begin_timestep -- first timestep to process
        !!!!       end_timestep -- final timestep to process
        !!!!       inter_timestep -- interval between timesteps
        INTEGER, INTENT(IN), DIMENSION(4) :: time
        REAL(REAL64), INTENT(IN) :: deltaz ! width of the bin is 2*deltaz, for number density
        REAL(REAL64), INTENT(IN) :: z1 ! range of sampling, full range is 2*z1, for number density
        INTEGER, INTENT(IN) :: Nz ! number of data points, in total 2*Nz+1 points, for number density
        INTEGER, INTENT(IN) :: dim2 ! second dimension of output, dim2 = 2*nhead + 2*ntail + 1
        INTEGER, INTENT(IN) :: dim3 ! third dimension of output, dim3 = 2*Nz + 1
        !!!! output: surface position [timesteps, xcolu, ycolu]
        !!!!         top or bot determined by input(10) (0-bot, 1-top)
        REAL(REAL64), DIMENSION(T,dim2,dim3), INTENT(OUT) :: output
        REAL(REAL64) :: xlo, xhi, ylo, yhi, zlo, zhi
        REAL(REAL64) :: xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound
        INTEGER :: i, j, k, s, m, l
        INTEGER :: n ! number of atoms
        INTEGER :: ix, iy
        REAL(REAL64) :: dx, dy, thickness_reduced
        REAL(REAL64) :: pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062
        INTEGER, DIMENSION(xcolu,ycolu) :: colu_xy
        REAL(REAL64), DIMENSION(:,:,:), ALLOCATABLE :: cell_xy
        REAL(REAL64), DIMENSION(nmol,nhead) :: posx_head, posy_head, posz_head
        REAL(REAL64), DIMENSION(nmol,ntail) :: posx_tail, posy_tail, posz_tail
        REAL(REAL64), DIMENSION(xcolu,ycolu) :: meanz_tail ! average height of tail beads at one timestep
        REAL(REAL64) :: trash1, trash2, trash3, tempx, tempy, tempz, trash5, trash6, trash4
        REAL(REAL64) :: xlen, ylen, zlen
        INTEGER :: skip_lines, counter
        REAL(REAL64), DIMENSION(:), ALLOCATABLE :: posz_temp, posz_ttemp
        INTEGER, DIMENSION(nmol) :: mol_direc ! direction of lipid, 1: upward, 0: downward
        REAL(REAL64) :: diffz ! difference between first head and last tail in z direction
        INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: pointer_head_top, pointer_head_bot ! point to which lipid
        INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: pointer_tail_top, pointer_tail_bot
        INTEGER, DIMENSION(nhead,xcolu,ycolu) :: head_top_count, head_bot_count
        INTEGER, DIMENSION(ntail,xcolu,ycolu) :: tail_top_count, tail_bot_count
        REAL(REAL64) :: rela_posz, meanz, current_z
        REAL(REAL64) :: deltaz_reduce, z1_reduce, dz_reduce
        INTEGER :: iz, izdown, izup, temp_atom_N

        OPEN(1, file=filename)

        n = nmol*(nhead+ntail)

        !! skip lines
        skip_lines = (time(2)-time(1))/time(4)
        DO i=1,skip_lines
            DO j=1,(9+n)
                READ(1,*)
            END DO
        END DO

        !! actual calculation
        !! for each timestep
        output = 0.0d0
        DO i=1,T
            ! skip first part of headers
            DO j=1,5
                READ(1,*)
            END DO
            ! read simulation box dimensions
            READ(1,*) xlo, xhi
            READ(1,*) ylo, yhi
            READ(1,*) zlo, zhi
            ! skip second part of headers
            READ(1,*)
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
            ! reduced deltaz and z1
            deltaz_reduce = deltaz/zlen
            z1_reduce = z1/zlen
            dz_reduce = z1_reduce/Nz
            !! allocate position vectors
            posx_head = 0
            posy_head = 0
            posz_head = 0
            posx_tail = 0
            posy_tail = 0
            posz_tail = 0
            !! calculate reduced thickness
            thickness_reduced = thickness/zlen
            IF (thickness_reduced.GE.0.5) THEN
                PRINT*, "membrane thickness greater than simulation cell! Reduce thickness"
            END IF
            DO j=1,nmol
                DO k=1,nhead
                    READ(1,*) trash1, trash2, trash3, tempx, tempy, tempz, trash5,&
                    & trash6, trash4
                    posx_head(j,k) = tempx * xlen + xlo
                    posy_head(j,k) = tempy * ylen + ylo
                    posz_head(j,k) = tempz
                END DO
                DO k=1,ntail
                    READ(1,*) trash1, trash2, trash3, tempx, tempy, tempz, trash5,&
                    & trash6, trash4
                    posx_tail(j,k) = tempx * xlen + xlo
                    posy_tail(j,k) = tempy * ylen + ylo
                    posz_tail(j,k) = tempz
                END DO
                diffz = zlen * (posz_head(j,1)-posz_tail(j,ntail))
                IF (ABS(diffz).GT.(zlen/2)) THEN
                    diffz = (zlen - ABS(diffz))*SIGN(1.0d0,-diffz)
                END IF
                IF (diffz.GT.0) THEN
                    mol_direc(j) = 1 !! lipid pointing upward
                ELSE
                    mol_direc(j) = 0 !! lipid pointing downward
                END IF
            END DO

            ! check nothing is out of simulation box, only for periodic boundaries
            IF (boundary(1).EQ.1) THEN
                DO j=1,nmol
                    DO k=1,nhead
                        IF ((posx_head(j,k).LT.xlo_bound) .OR. (posx_head(j,k).GT.xhi_bound)) THEN
                            posx_head(j,k) = posx_head(j,k) - sign((xhi_bound-xlo_bound),(posx_head(j,k)-xlo_bound))
                        END IF
                    END DO
                    DO k=1,ntail
                        IF ((posx_tail(j,k).LT.xlo_bound) .OR. (posx_tail(j,k).GT.xhi_bound)) THEN
                            posx_tail(j,k) = posx_tail(j,k) - sign((xhi_bound-xlo_bound),(posx_tail(j,k)-xlo_bound))
                        END IF
                    END DO
                END DO
            END IF
            IF (boundary(2).EQ.1) THEN
                DO j=1,nmol
                    DO k=1,nhead
                        IF ((posy_head(j,k).LT.ylo_bound) .OR. (posy_head(j,k).GT.yhi_bound)) THEN
                            posy_head(j,k) = posy_head(j,k) - sign((yhi_bound-ylo_bound),(posy_head(j,k)-ylo_bound))
                        END IF
                    END DO
                    DO k=1,ntail
                        IF ((posy_tail(j,k).LT.ylo_bound) .OR. (posy_tail(j,k).GT.yhi_bound)) THEN
                            posy_tail(j,k) = posy_tail(j,k) - sign((yhi_bound-ylo_bound),(posy_tail(j,k)-ylo_bound))
                        END IF
                    END DO
                END DO
            END IF

            meanz_tail = 0
            dx = (xhi_bound-xlo_bound)/REAL(xcolu)
            dy = (yhi_bound-ylo_bound)/REAL(ycolu)
            colu_xy = 0
            ALLOCATE(cell_xy(INT(n/xcolu/ycolu*4.0d0),xcolu,ycolu))
            cell_xy = 0
            DO j=1,nmol
                DO k=1,ntail
                    ix = INT(FLOOR((posx_tail(j,k)-xlo_bound)/dx)) + 1
                    iy = INT(FLOOR((posy_tail(j,k)-ylo_bound)/dy)) + 1
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
                    cell_xy(INT(colu_xy(ix,iy)),ix,iy) = posz_tail(j,k)
                END DO
            END DO
            DO j=1,xcolu
                DO k=1,ycolu
                    ALLOCATE(posz_temp(colu_xy(j,k)))
                    posz_temp(:) = cell_xy(1:colu_xy(j,k),j,k)
                    posz_temp = posz_temp*2.0d0*pi-pi
                    meanz_tail(j,k) = (ATAN2(SUM(SIN(posz_temp))/colu_xy(j,k),&
                                           &SUM(COS(posz_temp))/colu_xy(j,k))&
                                    &+pi)/2.0d0/pi
                    posz_temp = 0
                    counter = 0
                    IF ((meanz_tail(j,k)-thickness_reduced).LT.0.0d0) THEN
                        DO s=1,colu_xy(j,k)
                            IF ((cell_xy(s,j,k).GT.(meanz_tail(j,k)+thickness_reduced))&
                            &.AND.(cell_xy(s,j,k).LT.(1.0d0-(thickness_reduced-meanz_tail(j,k))))) THEN
                                colu_xy(j,k) = colu_xy(j,k) - 1
                            ELSE
                                counter = counter + 1
                                posz_temp(counter) = cell_xy(s,j,k)
                            END IF
                        END DO
                    ELSEIF ((meanz_tail(j,k)+thickness_reduced).GE.1.0d0) THEN
                        DO s=1,colu_xy(j,k)
                            IF ((cell_xy(s,j,k).GT.(thickness_reduced-(1.0d0-meanz_tail(j,k))))&
                            &.AND.(cell_xy(s,j,k).LT.(meanz_tail(j,k)-thickness_reduced))) THEN
                                colu_xy(j,k) = colu_xy(j,k) - 1
                            ELSE
                                counter = counter + 1
                                posz_temp(counter) = cell_xy(s,j,k)
                            END IF
                        END DO
                    ELSE
                        DO s=1,colu_xy(j,k)
                            IF ((cell_xy(s,j,k).GE.(meanz_tail(j,k)-thickness_reduced))&
                            &.AND.(cell_xy(s,j,k).LE.(meanz_tail(j,k)+thickness_reduced))) THEN
                                counter = counter + 1
                                posz_temp(counter) = cell_xy(s,j,k)
                            ELSE
                                colu_xy(j,k) = colu_xy(j,k) - 1
                            END IF
                        END DO
                    END IF
                    ALLOCATE(posz_ttemp(colu_xy(j,k)))
                    posz_ttemp(:) = posz_temp(1:colu_xy(j,k))
                    DEALLOCATE(posz_temp)
                    posz_ttemp = posz_ttemp*2.0d0*pi-pi
                    meanz_tail(j,k) = (ATAN2(SUM(SIN(posz_ttemp))/colu_xy(j,k),&
                                           &SUM(COS(posz_ttemp))/colu_xy(j,k))&
                                    &+pi)/2.0d0/pi
                    DEALLOCATE(posz_ttemp)
                END DO
            END DO
            DEALLOCATE(cell_xy)
            !! put top/bottom head/tail beads into bins
            ALLOCATE(pointer_head_top(nhead,INT(FLOOR(REAL(nmol/xcolu))),xcolu,ycolu))
            ALLOCATE(pointer_head_bot(nhead,INT(FLOOR(REAL(nmol/xcolu))),xcolu,ycolu))
            ALLOCATE(pointer_tail_top(ntail,INT(FLOOR(REAL(nmol/xcolu))),xcolu,ycolu))
            ALLOCATE(pointer_tail_bot(ntail,INT(FLOOR(REAL(nmol/xcolu))),xcolu,ycolu))
            head_top_count = 0
            head_bot_count = 0
            tail_top_count = 0
            tail_bot_count = 0
            DO j=1,nmol
                DO k=1,nhead
                    ix = INT(FLOOR((posx_head(j,k)-xlo_bound)/dx)) + 1
                    iy = INT(FLOOR((posy_head(j,k)-ylo_bound)/dy)) + 1
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
                    IF (mol_direc(j).EQ.1) THEN
                        head_top_count(k,ix,iy) = head_top_count(k,ix,iy) + 1
                        pointer_head_top(k,head_top_count(k,ix,iy),ix,iy) = j
                    ELSEIF (mol_direc(j).EQ.0) THEN
                        head_bot_count(k,ix,iy) = head_bot_count(k,ix,iy) + 1
                        pointer_head_bot(k,head_bot_count(k,ix,iy),ix,iy) = j
                    END IF
                END DO
                DO k=1,ntail
                    ix = INT(FLOOR((posx_tail(j,k)-xlo_bound)/dx)) + 1
                    iy = INT(FLOOR((posy_tail(j,k)-ylo_bound)/dy)) + 1
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
                    IF (mol_direc(j).EQ.1) THEN
                        tail_top_count(k,ix,iy) = tail_top_count(k,ix,iy) + 1
                        pointer_tail_top(k,tail_top_count(k,ix,iy),ix,iy) = j
                    ELSEIF (mol_direc(j).EQ.0) THEN
                        tail_bot_count(k,ix,iy) = tail_bot_count(k,ix,iy) + 1
                        pointer_tail_bot(k,tail_bot_count(k,ix,iy),ix,iy) = j
                    END IF
                END DO
            END DO
            DO j=1,xcolu
                DO k=1,ycolu
                    meanz = meanz_tail(j,k)
                    !! heads
                    DO s=1,nhead
                        !! top heads
                        DO m=1,head_top_count(s,j,k)
                            current_z = posz_head(pointer_head_top(s,m,j,k),s)
                            !! calculate relative position of bead to midplane of the bin
                            IF (meanz.GT.current_z) THEN
                                IF ((meanz-current_z).GT.0.5d0) THEN
                                    rela_posz = current_z + 1.0d0 - meanz
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            ELSE
                                IF ((current_z-meanz).GT.0.5d0) THEN
                                    rela_posz = meanz + 1.0d0 - current_z
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            END IF
                            !! [rela_posz-deltaz,rela_posz+deltaz]
                            IF ((rela_posz-deltaz_reduce).LT.(-z1_reduce)) THEN
                                IF ((rela_posz+deltaz_reduce).GT.(-z1_reduce)) THEN
                                    iz = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                    DO l=1,iz
                                        output(i,s,l) = output(i,s,l) + 1.0d0
                                        output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                    END DO
                                END IF
                            ELSE
                                IF ((rela_posz-deltaz_reduce).LT.(z1_reduce)) THEN
                                    IF ((rela_posz+deltaz_reduce).LT.z1_reduce) THEN
                                        izdown = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        izup = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                        DO l=izdown, izup
                                            output(i,s,l) = output(i,s,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    ELSE
                                        iz = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        DO l=iz,INT(2*Nz+1)
                                            output(i,s,l) = output(i,s,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    END IF
                                END IF
                            END IF
                        END DO
                        !! bottom heads
                        temp_atom_N = INT(nhead+s)
                        DO m=1,head_bot_count(s,j,k)
                            current_z = posz_head(pointer_head_bot(s,m,j,k),s)
                            !! calculate relative position of bead to midplane of the bin
                            IF (meanz.GT.current_z) THEN
                                IF ((meanz-current_z).GT.0.5d0) THEN
                                    rela_posz = current_z + 1.0d0 - meanz
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            ELSE
                                IF ((current_z-meanz).GT.0.5d0) THEN
                                    rela_posz = meanz + 1.0d0 - current_z
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            END IF
                            !! [rela_posz-deltaz,rela_posz+deltaz]
                            IF ((rela_posz-deltaz_reduce).LT.(-z1_reduce)) THEN
                                IF ((rela_posz+deltaz_reduce).GT.(-z1_reduce)) THEN
                                    iz = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                    DO l=1,iz
                                        output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                        output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                    END DO
                                END IF
                            ELSE
                                IF ((rela_posz-deltaz_reduce).LT.(z1_reduce)) THEN
                                    IF ((rela_posz+deltaz_reduce).LT.z1_reduce) THEN
                                        izdown = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        izup = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                        DO l=izdown, izup
                                            output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    ELSE
                                        iz = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        DO l=iz,INT(2*Nz+1)
                                            output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    END IF
                                END IF
                            END IF
                        END DO
                    END DO
                    !! tails
                    DO s=1,ntail
                        !! top tails
                        temp_atom_N = INT(2*nhead+s)
                        DO m=1,tail_top_count(s,j,k)
                            current_z = posz_tail(pointer_tail_top(s,m,j,k),s)
                            !! calculate relative position of bead to midplane of the bin
                            IF (meanz.GT.current_z) THEN
                                IF ((meanz-current_z).GT.0.5d0) THEN
                                    rela_posz = current_z + 1.0d0 - meanz
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            ELSE
                                IF ((current_z-meanz).GT.0.5d0) THEN
                                    rela_posz = meanz + 1.0d0 - current_z
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            END IF
                            !! [rela_posz-deltaz,rela_posz+deltaz]
                            IF ((rela_posz-deltaz_reduce).LT.(-z1_reduce)) THEN
                                IF ((rela_posz+deltaz_reduce).GT.(-z1_reduce)) THEN
                                    iz = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                    DO l=1,iz
                                        output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                        output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                    END DO
                                END IF
                            ELSE
                                IF ((rela_posz-deltaz_reduce).LT.(z1_reduce)) THEN
                                    IF ((rela_posz+deltaz_reduce).LT.z1_reduce) THEN
                                        izdown = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        izup = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                        DO l=izdown, izup
                                            output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    ELSE
                                        iz = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        DO l=iz,INT(2*Nz+1)
                                            output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    END IF
                                END IF
                            END IF
                        END DO
                        !! bottom tails
                        temp_atom_N = INT(2*nhead+ntail+s)
                        DO m=1,tail_bot_count(s,j,k)
                            current_z = posz_tail(pointer_tail_bot(s,m,j,k),s)
                            !! calculate relative position of bead to midplane of the bin
                            IF (meanz.GT.current_z) THEN
                                IF ((meanz-current_z).GT.0.5d0) THEN
                                    rela_posz = current_z + 1.0d0 - meanz
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            ELSE
                                IF ((current_z-meanz).GT.0.5d0) THEN
                                    rela_posz = meanz + 1.0d0 - current_z
                                ELSE
                                    rela_posz = current_z - meanz
                                END IF
                            END IF
                            !! [rela_posz-deltaz,rela_posz+deltaz]
                            IF ((rela_posz-deltaz_reduce).LT.(-z1_reduce)) THEN
                                IF ((rela_posz+deltaz_reduce).GT.(-z1_reduce)) THEN
                                    iz = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                    DO l=1,iz
                                        output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                        output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                    END DO
                                END IF
                            ELSE
                                IF ((rela_posz-deltaz_reduce).LT.(z1_reduce)) THEN
                                    IF ((rela_posz+deltaz_reduce).LT.z1_reduce) THEN
                                        izdown = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        izup = INT(CEILING((rela_posz+deltaz_reduce+z1_reduce)/dz_reduce))
                                        DO l=izdown, izup
                                            output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    ELSE
                                        iz = INT(CEILING((rela_posz-deltaz_reduce+z1_reduce)/dz_reduce)+1)
                                        DO l=iz,INT(2*Nz+1)
                                            output(i,temp_atom_N,l) = output(i,temp_atom_N,l) + 1.0d0
                                            output(i,dim2,l) = output(i,dim2,l) + 1.0d0
                                        END DO
                                    END IF
                                END IF
                            END IF
                        END DO
                    END DO
                END DO
            END DO
            output(i,:,:) = output(i,:,:)/xlen/ylen/deltaz/2
            DEALLOCATE(pointer_head_top)
            DEALLOCATE(pointer_head_bot)
            DEALLOCATE(pointer_tail_top)
            DEALLOCATE(pointer_tail_bot)

            ! CALL CPU_TIME(end_time)
            ! PRINT*, end_time-begin_time

        END DO




        CLOSE(1)

    END SUBROUTINE

    !!!! This subroutine calculates probability distribution of displacement
    !!!! 
    SUBROUTINE diffusion_distribution(filename,nhead,ntail,T,ngap,nmol,time,nbins,sstop,output)
        !! for meaning of variables, see calc_surf.calc_midplane_tail_onefile
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(IN) :: nhead, ntail, nmol
        !! note, here T = (end_timestep-begin_timestep)/inter_timestep/(ngap+1)
        INTEGER, INTENT(IN) :: T ! size of output
        INTEGER, INTENT(IN) :: ngap ! gap between timesteps
        !!!! time: init_timestep -- initial timestep of the simulation
        !!!!       begin_timestep -- first timestep to process
        !!!!       end_timestep -- final timestep to process
        !!!!       inter_timestep -- interval between timesteps
        INTEGER, INTENT(IN) :: nbins !! number of bins for probability distribution
        REAL(REAL64), INTENT(IN) :: sstop !! maximum displacement considered
        INTEGER, INTENT(IN), DIMENSION(4) :: time
        REAL(REAL64), DIMENSION(T,nbins), INTENT(OUT) :: output
        REAL(REAL64), DIMENSION(nmol,2) :: lipid_pos_prev !! previous lipid position, only x and y
        REAL(REAL64), DIMENSION(nmol,2) :: lipid_pos_now !! current lipid position, only x and y
        REAL(REAL64), DIMENSION(nhead,2) :: head_pos !! position of head beads, only x and y
        REAL(REAL64), DIMENSION(ntail,2) :: tail_pos !! position of tail beads, only x and y
        REAL(REAL64) :: displacement2 !! square of lateral displacement
        REAL(REAL64) :: dss
        INTEGER :: i, j, k, n, skip_lines, is
        REAL(REAL64) :: trash1, trash2, trash3, trash4, trash5, trash6, tempx, tempy

        OPEN(1, file=filename)

        n = nmol*(nhead+ntail)
        dss = sstop/nbins
        output = 0

        !! skip lines
        skip_lines = (time(2)-time(1))/time(4)
        DO i=1,skip_lines
            DO j=1,(9+n)
                READ(1,*)
            END DO
        END DO

        !! read the position of first time step
        ! skip header
        DO j=1,9
            READ(1,*)
        END DO
        ! calculate lipid position
        DO j=1,nmol
            DO k=1,nhead
                READ(1,*) trash1, trash2, trash3, trash4, trash5,trash6, tempx, tempy
                head_pos(k,1) = tempx
                head_pos(k,2) = tempy
            END DO
            DO k=1,ntail
                READ(1,*) trash1, trash2, trash3, trash4, trash5, trash6, tempx, tempy
                tail_pos(k,1) = tempx
                tail_pos(k,2) = tempy
            END DO
            lipid_pos_prev(j,1) = (SUM(head_pos(:,1))+SUM(tail_pos(:,1)))/(nhead+ntail)
            lipid_pos_prev(j,2) = (SUM(head_pos(:,2))+SUM(tail_pos(:,2)))/(nhead+ntail)
        END DO


        !! actual calculation begin
        DO i=1,T
            ! skip ngap timesteps
            DO j=1,ngap*(n+9)
                READ(1,*)
            END DO
            ! now read the timestep of interest
            ! skip first part of headers
            DO j=1,9
                READ(1,*)
            END DO
            ! calculate lipid position and calculate displacement
            DO j=1,nmol
                DO k=1,nhead
                    READ(1,*) trash1, trash2, trash3, trash4, trash5,trash6, tempx, tempy
                    head_pos(k,1) = tempx
                    head_pos(k,2) = tempy
                END DO
                DO k=1,ntail
                    READ(1,*) trash1, trash2, trash3, trash4, trash5, trash6, tempx, tempy
                    tail_pos(k,1) = tempx
                    tail_pos(k,2) = tempy
                END DO
                lipid_pos_now(j,1) = (SUM(head_pos(:,1))+SUM(tail_pos(:,1)))/(nhead+ntail)
                lipid_pos_now(j,2) = (SUM(head_pos(:,2))+SUM(tail_pos(:,2)))/(nhead+ntail)
                displacement2 = (lipid_pos_prev(j,1)-lipid_pos_now(j,1))**2 &
                                  +(lipid_pos_prev(j,2)-lipid_pos_now(j,2))**2
                ! update previous lipid positions
                lipid_pos_prev(j,1) = lipid_pos_now(j,1)
                lipid_pos_prev(j,2) = lipid_pos_now(j,2)
                ! put displacement into bins
                IF (displacement2>sstop) THEN
                    PRINT*,"timestep",i,"lipid",j, "displacement2 bigger than sstop"
                ELSE
                    is = INT(FLOOR(displacement2/dss)) + 1
                    output(i,is) = output(i,is) + 1
                END IF
            END DO
            output(i,:) = output(i,:)/nmol
        END DO

        CLOSE(1)
    END SUBROUTINE


END MODULE
