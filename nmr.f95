program nmr

!!!!!! Decare variables
parameter (MAXPTS = 10000)

!!! Variables for reading in
character (len = 40) :: inputFile, outputFile
double precision :: baseline, tolerance
integer :: filtertype, filterSize, filterPasses, integrationChoice


!!! variables for insertion sort
double precision :: xpts(0:MAXPTS), ypts(0:MAXPTS)
double precision :: xtmp, ytmp
integer :: count


!!! variables for spline
double precision :: bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS), x(0:MAXPTS)

double precision :: hld


!!! Vars for integration and output
double precision :: lowIntPts(1:10), upIntPts(1:10),intAreas(1:10)
integer :: peaks


!!!!!!!Read in instructions and store values

open (unit = 15, file='nmr.in', status = 'old')
    read (15, *) inputFile
    read (15, *) baseline 
    read (15, *) tolerance 
    read (15, *) filtertype 
    read (15, *) filterSize 
    read (15, *) filterPasses 
    read (15, *) integrationChoice 
    read (15, *) outputFile 
close(15)

!!!!! Insertion sort
count = -1
open (unit = 16, file = inputFile, status = 'old')

100 read (16, *, end = 200) xtmp, ytmp
count = count + 1
    do i = 0, count
        if (xtmp < xpts(i)) then
            do j = count, i, -1
                xpts(j + 1) = xpts(j)
                ypts(j + 1) = ypts(j)
            enddo

            xpts(i) = xtmp
            ypts(i) = ytmp
            goto 100
        end if
    enddo

    xpts(count) = xtmp
    ypts(count) = ytmp


goto 100 
200 continue
close(16)

!print *, xpts(0:count)
call TMshift
!print *, xpts(0:count)

if (filtertype /= 0 .and. filterSize /= 0) then
    if (filtertype == 1) then
        do l = 1, filterPasses
            print *, "pass ", l
            call boxFilter
        end do
    else if (filtertype == 2) then
    end if
else
    print *, "Filter is off"
end if



call buildSpline(xpts, ypts, bvals, cvals, dvals, count)


!do i =0, 2
!  hld = xpts(i)
!  print*, "Xpts i = ", hld
!  print*, calcSpline(hld)
!end do


open (unit=17, file=outputfile, status='unknown')
  hld = xpts(0)
  do while (hld <= xpts(count))
    write (17,*) hld, calcSpline(hld)
    hld = hld + (abs(xpts(count) - xpts(0)) / 10000.0D0)
  end do
close(17)

call findIntPoints
do i =1, peaks
  print *, "Peak at: ", lowIntPts(i) , upIntPts(i)
end do



call integrate


do i=1, peaks
  print *, intAreas(i)
end do




contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!! TMS Shift !!!!!!!!!!!!!!!!

subroutine TMShift
        logical :: run
        integer :: p , ac, l
        double precision :: arry(0:500),arrx(0:500), highY, shift
        
        highY = 0
        run = .true.
        ac = 500
        p = count
        do while (ypts(p) < baseline)      
                p = p - 1
        enddo
        
         do while (run)
                if (ypts(p) > baseline) then
                        arry(ac) = ypts(p)
                        arrx(ac) = xpts(p)
                        ac = ac - 1
                        !print *, "Ac ", ac
                else
                        run = .false.
                end if
                p = p - 1
         enddo


        do l = ac, 500, 1
                !print *, "Arry(ac)", arry(l), l
                if (arry(l) > highY) then
                        !print *, "High ", highY
                        highY = arry(l)
                        shift = arrx(l)
                endif                
        enddo
        
        do l = 0, count,1
                xpts(l) = xpts(l) - shift
        end do
        

end subroutine TMShift





!!!!!!!!!!!!! Box Car !!!!!!!!!!!!!!!!!!!

subroutine boxFilter

    !integer, intent (in) :: filterSize, count
    !double precision, intent (inout) :: ypts(0:MAXPTS)
    double precision :: ysum, temp(0:MAXPTS)
    integer :: offset, i, j, l, k, lowBound, upBound

    temp(0:) = ypts(0:)
    offset = filterSize / 2
    
    ysum = 0.0D0

    iloop: do i = 0, count, 1
        ysum = 0
        lowBound = i - offset
        upBound = i + offset
        
        if (lowBound < 0) then
            do k = lowBound , -1, 1
                ysum = ysum + ypts(count + k + 1)
                lowBound = lowBound + 1
            end do
        end if
        
        if (upBound > count) then
           
            do l = upBound - count - 1, 0, -1
                ysum = ysum + ypts(l)
                upBound = upBound - 1
            end do 
        end if


        do j = lowBound, upBound, 1
            ysum = ysum + ypts(j)
        end do

        ysum = ysum / real(filterSize, kind = 8)
        
        temp(i) = ysum
    end do iloop

    ypts(0:) = temp(0:)
end subroutine boxFilter





!!!!!!!!!! SG Filter !!!!!!!!!!!!!!

subroutine SGFilter

double precision :: m
integer :: i

m = count - (filterSize - 1)

do i = 2, filterSize
end do


end subroutine SGFilter



!!!!!!!!!! Build Spline !!!!!!!!!!!!!

subroutine buildSpline(xpts, ypts, bvals, cvals, dvals, count)


double precision, intent (out) :: bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)
double precision, intent (in) :: xpts(0:MAXPTS), ypts(0:MAXPTS)
double precision :: alvals(0:MAXPTS), lvals(0:MAXPTS)
double precision :: muvals(0:MAXPTS),zvals(0:MAXPTS)
double precision :: hvals(0:MAXPTS)
integer, intent (in) :: count
integer :: i, j, l

do i = 0, count
    hvals(i) = xpts(i + 1)- xpts(i)
end do

do i = 1, count, 1
    alvals(i) = ((3/hvals(i))*(ypts(i+1)-ypts(i))) - ((3/hvals(i-1))*(ypts(i)-ypts(i-1))) 
end do


lvals(0) = 1
muvals(0) = 0
zvals(0) = 0


do i = 1, count, 1
    lvals(i) = 2*(xpts(i+1)-xpts(i-1))-hvals(i-1)*muvals(i-1)
    muvals(i) = hvals(i) / lvals(i)
    zvals(i) = (alvals(i) - hvals(i-1)*zvals(i-1)) / lvals(i)
end do

lvals(count) = 1
zvals(count) = 0
cvals(count) = 0

do j = count, 0, -1
    cvals(j) = zvals(j)- muvals(j)*cvals(j+1)
    bvals(j)= (ypts(j+1)-ypts(j))/hvals(j) - hvals(j) * (cvals(j+1) + 2*cvals(j))/3
    dvals(j) = (cvals(j+1) - cvals(j)) / (3 * hvals(j))
end do

!do l = 0, count - 1
    !print *, ypts(l), bvals(l), cvals(l), dvals(l)
!end do

end subroutine buildSpline



!!!!!!!!!!!!Find int points !!!!!!!!!!!!!

subroutine findIntPoints
integer :: peakCount, i

peakCount = 1
i = 1
do while (xpts(i) < 0.0D0)
  if (ypts(i-1)>baseline .and. ypts(i)<baseline) then
    print *, "Upper end"
    upIntPts(peakCount) = bisection(xpts(i-1),xpts(i))
    peakCount = peakCount + 1
  else if(ypts(i-1)<baseline .and. ypts(i)>baseline) then
    print *, "Lower end"
    lowIntPts(peakCount) = bisection(xpts(i-1),xpts(i))
  end if
  
  i = i + 1
end do

peaks = peakCount - 1

end subroutine findIntPoints


subroutine integrate
integer :: k

do k = 1, peaks, 1
  print *, "Running"
  if (integrationChoice == 0) then
    intAreas(k) = CNC(lowIntPts(k), upIntPts(k))
  else if (integrationChoice == 1) then
    intAreas(k) = romberg(lowIntPts(k), upIntPts(k))
  else if (integrationChoice == 2) then
    intAreas(k) = adquad(lowIntPts(k), upIntPts(k), tolerance)
  else
  end if
end do

end subroutine integrate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Functions !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!! Calc Spline !!!!!!!!!!!!!!!!!

double precision function calcSpline(find)!, xpts, ypts, bvals, cvals, dvals, count)


!double precision, intent (in) :: xpts(0:MAXPTS),ypts(0:MAXPTS), bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)
double precision, intent (in) :: find
double precision :: ans
integer :: eqchoose, i


!print *, "Calc spline for : ", find

eqchoose = -1
do i = 0, count
    if (find >= xpts(i)) eqchoose = eqchoose + 1
end do

!print *, "Using equation: ", eqchoose

ans = ypts(eqchoose) + bvals(eqchoose) * (find - xpts(eqchoose))
      ans = ans + cvals(eqchoose) * ((find - xpts(eqchoose))**2)
      ans = ans + dvals(eqchoose) * ((find - xpts(eqchoose))**3)

calcSpline = ans

end function calcSpline


!!!!!!!!!!!!!! Bisection !!!!!!!!!!!!!!!

double precision function bisection (lowerBound, upperBound)

!double precision, intent (in) :: xpts(0:MAXPTS),ypts(0:MAXPTS)
double precision, intent (in) :: lowerBound, upperBound
double precision :: tLower, tUpper, guess, prevAnswer, answer, difference, lval, uval
integer :: i

tLower = lowerBound
tUpper = upperBound
guess = (upperBound + lowerBound ) / 2
lval = calcSpline(tLower)
uval = calcSpline(tUpper)
prevAnswer = 0
answer = calcSpline(guess)
difference = 1



do while(difference > tolerance)

    if (lval > baseline .and. answer > baseline) then
        tLower = guess
    else if (lval < baseline .and. answer < baseline) then
        tLower = guess
    else
        tUpper = guess
    end if

    guess = (tLower + tUpper ) / 2.0D0

    lval = calcSpline(tLower)
    uval = calcSpline(tUpper)

    answer = calcSpline(guess)
    difference = abs(prevAnswer - answer)
    prevAnswer = answer
end do
    
    bisection = guess

end function bisection




!!!!!!!!!! Simpson Integration !!!!!!!!!!!

double precision function simpson(low, up)
double precision, intent (in):: up, low
double precision :: h

h = (up - low) / 2

simpson = (h/3.0D0) * ((calcSpline(low) - baseline) + (4.0D0 * (calcSpline(low+h)-baseline) ) + (calcSpline(up)-baseline))

end function simpson


!!!!!!!! Compotite Newton Cotes !!!!!!!!!!!

double precision function CNC(low, up)

double precision, intent (in):: low,up
double precision :: h, prev, cur, ends, evens, odds, ihold, tol
integer :: run, counter, i

tol = tolerance
counter = 0

cur = 1.0D0
prev = 0.0D0
ends = calcSpline(low) + calcSpline(up) - 2*baseline
run = 3


!print *, "Diff: ", cur - prev, tol

do while (cur - prev >= tol)
  print *, "Calc"
  prev = cur
  evens = 0.0D0
  odds = 0.0D0
  h = (up - low) / run
  
  do i = 1, run - 1
    ihold = i
    if (mod(i,2)==0) then
    !print *, "Evens: ", evens, " plus ", calcSpline(low + (ihold*h)), " @ ", low + (ihold*h)
      evens = evens + ( calcspline(low + (ihold*h)) - baseline )
      
    else
    !print *, "Odds: ", odds, " plus ", calcSpline(low + (ihold*h)), " @ ", low + (ihold*h)
      odds = odds + ( calcSpline(low + (ihold*h)) - baseline)
    end if
  end do
  
  cur = (h* (ends + (2 * evens) + (4* odds)) )/3.0D0
  !print *, "Cur | ", cur, " | Prev | ", prev, " | Diff | ", cur - prev
  run = run + 2
  !counter = counter + 1
end do

CNC = cur
end function CNC


!!!!!!!!! Romberg !!!!!!!!!!!!!!!!

double precision function romberg(low, up)


double precision, intent (in) :: low, up
double precision :: intsum, ans
double precision :: h, arr(1:2,1:10)
integer :: i
logical :: run


ans = -1.0D0

run = .true.


h = up - low

arr(1,1) =(h/2) * (calcSpline(low) + calcSpline(up))

i = 2
do while (i <= 10 .and. run )
    !print *, "Tol check : ", arr(1,i-1) - arr(2,i)
    intsum = 0.0D0
    do k = 1, 2**(i-2)
        intsum = intsum + ( calcSpline(low + (k-0.5D0)*h) - baseline)
        !print *, "Sum: ", intsum
    end do
    
    arr(2,1) = (0.5D0) * ( arr(1,1) + h * intsum)
    !print *, "new root: ", arr(2,1)
    do j = 2, i
        arr(2,j) = arr(2,j-1) + ( (arr(2,j-1) - arr(1,j-1)) / (4**(j-1) -1) )
    end do

    !print *, arr(2,1:i)


    if (abs(arr(1,i-1) - arr(2,i) ) < tolerance) run = .false.

    do l=1,10
        arr(1,l) = arr(2,l)
    end do

    romberg = arr(1,i)

    h = h / 2
    i = i + 1

end do

end function romberg





!!!!!!!!! Adaptive Quadurature !!!!!!!!!!!!!!

double precision recursive function adquad(low, up, tol) result (ans)
double precision, intent (in) :: low,up,tol
double precision :: mid, cur, halves

!print *, "Run"
mid =( up + low ) / 2
cur = simpson(low,up)
halves = simpson(low,mid) + simpson(mid,up)

!print *, cur , halves

if (cur - halves < tol) then
ans= halves
else
ans = adquad(low, mid, tol/2) + adquad(mid, up, tol/2)
end if


end function adquad


end program nmr
