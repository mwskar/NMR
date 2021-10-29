program rom

double precision :: low, up, calcSpline
double precision :: intsum, ans
double precision :: h, arr(0:10,0:10)
integer :: i

low = 2
up = 4
tolerance = 1E-2
ans = -1.0D0

i = 1

h = low - up
arr(0,0) = (h/2) * (calcSpline(low) + calcSpline(up))

do while ( (arr(1,i-1) - arr(2,i) ) >= tolerance )
    intsum = 0.0D0
    do k = 1, 2**(i-2)
        intsum = intsum + calcSpline(low + (k-0.5D0)*h)
    end do
    
    arr(2,0) = (1/2) * ( arr(i-1,i-1) + h * intsum)

    do j = 1, i
        arr(2,j) = arr(2,j-1) + ( (arr(2,j-1) - arr(1,j-1)) / (4**(j-1) -1) )
    end do
    i = i +1

    do l=0,10
        arr(1,l) = arr(2,l)
    end do

end do

contains

double precision function calcSpline(x)

    double precision, intent (in) :: x

    calcSpline = x**2

end function calcSpline


end program rom