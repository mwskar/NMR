program gtester

double precision :: low,up
double precision :: shift, tolerance, cur, prev, root, weight
character (len = 8) :: fileArr(1:5)
integer :: fileChoice

up = 4.0D0
low = 2.0D0
tolerance = 1E-9

fileArr(1:5) = (/'02pt.dat','04pt.dat','08pt.dat','16pt.dat','32pt.dat'/)

cur = 0.0D0
prev = 1.0D0

fileChoice = 1

do while (abs(cur - prev)>=tolerance) 
  prev = cur
  cur = 0.0D0
  open (unit = 100, file = fileArr(fileChoice), status = 'old')
101    read(100,*, end = 102) root, weight
          print *, root, weight
          shift = ( ( (up-low)*root ) + (low + up) )/ 2.0D0
          print *, "Shift = ", shift
          cur = cur + ( calcSpline(shift) * weight)
       goto 101
102 continue
  close(100)
  cur = cur * ( (up-low)/2.0D0 )
  print *, cur
  fileChoice = fileChoice + 1
end do

contains

double precision function calcSpline(x)
double precision, intent (in) :: x
calcSpline = x**2
end function calcSpline

end program gtester
