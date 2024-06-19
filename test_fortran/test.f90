program main
      implicit none

      integer :: i, n
      REAL, DIMENSION(:),ALLOCATABLE :: a, b

      n = 10
      ALLOCATE(a(n))
      ALLOCATE(b(n))
      a = 0.0d0
      do i = 2,n
          !print*, "test do loop", i
          a(i) = REAL(i)
      end do

      !print*, a/3.0d0

      !print*, "test floor", ceiling(-3.5d0)

      !print*, mod(-7,3)

      b = 0
      b = a
      print*, b
      a = a*3.0d0
      print*, b

      stop "12"
      



end program
