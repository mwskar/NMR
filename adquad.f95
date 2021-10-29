program tester


double precision :: tol, PI

PI = 4.D0*datan(1.D0)

tol = 1E-5

print*, adquad(2.0D0,40.D0,10 *  tol)

!print*, "Hi"


contains




double precision function ncInt(low,up)

double precision, intent (in) :: low,up
double precision :: h

h = (up - low ) / 2
ncInt = (h/3) * (calcSpline(low) + (4* calcSpline(low + h)) + calcSpline(up))

end function ncInt



double precision function calcSpline(x)
double precision, intent (in) :: x

calcSpline = (x**2) + x + 5
end function calcSpline


double precision recursive function adquad(low, up, tol) result (ans)
double precision, intent (in) :: low,up,tol
double precision :: mid, cur, halves

!print *, "Run"
mid =( up + low ) / 2
cur = ncInt(low,up)
halves = ncInt(low,mid) + ncInt(mid,up)

!print *, cur , halves

if (cur - halves < tol) then
ans= halves
else
ans = adquad(low, mid, tol/2) + adquad(mid, up, tol/2)
end if


end function adquad



end program tester




