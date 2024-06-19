MODULE readBins

!USE kinds

IMPLICIT NONE

CONTAINS

  SUBROUTINE read_bins(bins, time_marks, nbin_x, nbin_y, nbin_z, filename)
    INTEGER, PARAMETER :: REAL32 = SELECTED_REAL_KIND(6, 37)
    CHARACTER(len=*), INTENT(IN) :: filename ! name of file to read from
    INTEGER, INTENT(IN) :: time_marks ! number of timesteps
    INTEGER, INTENT(IN) :: nbin_x, nbin_y, nbin_z ! dimension of bins
    REAL(REAL32), INTENT(OUT), DIMENSION(time_marks,nbin_x,nbin_y,nbin_z) :: bins ! the bins
    INTEGER :: i, k
    INTEGER :: block ! block size, number of lines for each timestep
    INTEGER :: tempZ, tempY, tempX ! tells you which x, y, z coordinate you are in
    INTEGER :: first
    REAL(REAL32) :: second, third, fourth, fifth, sixth

    block = nbin_x*nbin_y*nbin_z + 1
    !time_marks = int((n_lines-3)/(block))
    block = block - 1
    OPEN(1, file = filename)
    DO i=1,3
      READ(1,*)
    END DO
    DO k = 1, time_marks
      READ(1,*)
      DO i = 1, block
        tempZ = INT(MOD(i-1,nbin_z) + 1)
        tempY = INT(MOD(INT((i-tempZ)/nbin_z),nbin_y) + 1)
        tempX = INT((i-tempZ-nbin_y*(tempY-1))/nbin_z/nbin_y + 1)
        READ(1,*) first, second, third, fourth, fifth, sixth
        bins(k,tempX,tempY,tempZ) = fifth
      END DO
    END DO


    CLOSE(1)

  END SUBROUTINE

END MODULE
