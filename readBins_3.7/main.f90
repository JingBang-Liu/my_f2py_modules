PROGRAM main

!USE kinds
USE readBins

IMPLICIT NONE

INTEGER, PARAMETER :: REAL32 = SELECTED_REAL_KIND(6, 37)
INTEGER :: nbin_x, nbin_y, nbin_z
INTEGER :: timemarks
CHARACTER(len=100) :: filename
REAL(REAL32), DIMENSION(:,:,:,:), ALLOCATABLE :: bins

nbin_x = 153
nbin_y = 4
nbin_z = 8
timemarks = 500
filename = "../../MD/MD_results/rupture_1/1/bins_dens.dat"
ALLOCATE(bins(timemarks,nbin_x,nbin_y,nbin_z))

CALL read_bins(bins,timemarks,nbin_x,nbin_y,nbin_z,filename)
PRINT*, "THE SIZE OF BINS", SHAPE(bins) 
PRINT*, bins(100,20,2,3)

END PROGRAM