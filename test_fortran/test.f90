program main

    use kinds

    implicit none

    integer :: i, n
    REAL(REAL64), DIMENSION(:),ALLOCATABLE :: a, b
    CHARACTER(LEN=100) :: filename
    REAL(REAL64), DIMENSION(1) :: xlo, xhi, ylo, yhi, zlo, zhi
    REAL(REAL64), DIMENSION(6) :: input

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
    
    filename = "equ.dat"
    OPEN(1, file=filename)

    DO i=1,5
        READ(1,*)
    END DO
    READ(1,*) input(1), input(2)
    READ(1,*) input(3), input(4)
    READ(1,*) input(5), input(6)
    PRINT*, input
    CLOSE(1)

    stop "12"
end program
