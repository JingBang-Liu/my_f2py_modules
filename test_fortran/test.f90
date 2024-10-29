program main

    use kinds

    implicit none

    integer :: i, n
    REAL(REAL64), DIMENSION(:),ALLOCATABLE :: a, b
    CHARACTER(LEN=100) :: filename
    REAL(REAL64), DIMENSION(1) :: xlo, xhi, ylo, yhi, zlo, zhi
    REAL(REAL64), DIMENSION(6) :: input
    REAL(REAL64) :: pi = 3.14159265359
    REAL(REAL64), DIMENSION(2) :: x_sin
    REAL(REAL64), DIMENSION(2,3) :: temp

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

    DO i=1,2
        PRINT*, "HAHA"
    ENDDO

    x_sin(1) = pi*0.5d0
    x_sin(2) = pi
    PRINT*, SIN(x_sin)*SIN(x_sin)
    PRINT*, ATAN2(SIN(x_sin),COS(x_sin))
    PRINT*, MOD(-0.2d0,1.0d0)

    temp(1,1) = 1
    temp(1,2) = 2
    temp(1,3) = 3
    temp(2,1) = 4
    temp(2,2) = 5
    temp(2,3) = 6
    PRINT*, SUM(temp,1), "lala"

    stop "12"

end program
