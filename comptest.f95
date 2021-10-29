program composite


double precision:: low,up, tol
double precision :: h, prev, cur, ends, evens, odds, ihold
integer :: run, counter


counter = 0

low = 2.0D0
up = 4.0D0
tol = 1E-5

cur = 1.0D0
prev = 0.0D0
ends = calcSpline(low) + calcSpline(up)
run = 3


!print *, "Diff: ", cur - prev, tol

do while (cur - prev >= tol)
  !print *, "Run"
  prev = cur
  evens = 0.0D0
  odds = 0.0D0
  h = (up - low) / run
  
  do i = 1, run - 1
    ihold = i
    if (mod(i,2)==0) then
    !print *, "Evens: ", evens, " plus ", calcSpline(low + (ihold*h)), " @ ", low + (ihold*h)
      evens = evens + calcspline(low + (ihold*h))
      
    else
    !print *, "Odds: ", odds, " plus ", calcSpline(low + (ihold*h)), " @ ", low + (ihold*h)
      odds = odds + calcSpline(low + (ihold*h))
    end if
  end do
  
  cur = (h* (ends + (2 * evens) + (4* odds)) )/3.0D0
  !print *, "Cur | ", cur, " | Prev | ", prev, " | Diff | ", cur - prev
  run = run + 2
  !counter = counter + 1
end do


print *, cur

contains 


double precision function calcSpline(x)
double precision, intent (in) :: x

calcSpline = x**2
end function calcSpline


end program composite