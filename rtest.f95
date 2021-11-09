program rom

double precision :: low, up
double precision :: intsum, ans
double precision :: h, arr(1:2,1:10)
integer :: i
logical :: run


low = 2
up = 4
tolerance = 1E-5
ans = -1.0D0

run = .true.


h = up - low

arr(1,1) =(h/2) * (calcSpline(low) + calcSpline(up))

!print *, "Test case: ", arr(1,i-1) - arr(2,i)


i = 2
do while (i <= 10 .and. run )
    print *, "Tol check : ", arr(1,i-1) - arr(2,i)
    intsum = 0.0D0
    do k = 1, 2**(i-2)
        intsum = intsum + calcSpline(low + (k-0.5D0)*h)
        !print *, "Sum: ", intsum
    end do
    
    arr(2,1) = (0.5D0) * ( arr(1,1) + h * intsum)
    !print *, "new root: ", arr(2,1)
    do j = 2, i
        arr(2,j) = arr(2,j-1) + ( (arr(2,j-1) - arr(1,j-1)) / (4**(j-1) -1) )
    end do

    print *, arr(2,1:i)


    if (abs(arr(1,i-1) - arr(2,i) ) < tolerance) run = .false.

    do l=1,10
        arr(1,l) = arr(2,l)
    end do


    h = h / 2
    i = i + 1

end do

contains

double precision function calcSpline(x)

    double precision, intent (in) :: x

    calcSpline = x**2

end function calcSpline


end program rom
