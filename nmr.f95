program nmr

!!!!!! Decare variables
parameter (MAXPTS = 10000)

!!! Variables for reading in
character (len = 40) :: inputFile, outputFile
double precision :: baseline, tolerance
integer :: filtertype, filterSize, filterPasses, integrationChoice


!!! Variables for shift
double precision :: plotShift


!!! variables for insertion sort
double precision :: xpts(0:MAXPTS), ypts(0:MAXPTS)
double precision :: xtmp, ytmp
integer :: count


!!! variables for spline
double precision :: bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS), x(0:MAXPTS)

double precision :: hld


!!! Variables for integration, peaks, hydorgens, and output
double precision :: lowIntPts(1:10), upIntPts(1:10),intAreas(1:10), peakTop(1:10), peakPos(1:10)
integer :: peaks, hydroArr(1:10)



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


!!!! Find the 'TMS' location and shift the data set

call TMshift


!!! Run filter option for the number of times input

if (filtertype /= 0 .and. filterSize /= 0) then
    if (filtertype == 1) then
        do l = 1, filterPasses
            call boxFilter
        end do
    else if (filtertype == 2) then
      do l = 1, filterpasses
        call SGFilter
      end do
    end if
else
    print *, "Filter is off"
end if



!!! Uses the data to build a natural cubic spline

call buildSpline(xpts, ypts, bvals, cvals, dvals, count)

!!! Searches through the data set for the begining and end of peaks

call findIntPoints

!!! Integrates the peaks with the users choice of
!!! method and tolerance

call integrate

!!! Finds the top loactions and values for each peak

call findTops


!!! Calculates the number of hydrogens in each peak

call calcHydrogens



!!!! Presents information in text output and writes to output file
!open (unit = 16, file = inputFile, status = 'old')

open(unit = 102, file = outputFile, status = 'unknown' )

print*, "Program options"
write(102,*) "Program options"
print*, "================================"
write(102,*) "================================"
print*, "Baseline adjustment  :", baseline
write(102,*) "Baseline adjustment  :", baseline
print*, "Tolerance            :", tolerance
write(102,*) "Tolerance            :", tolerance


if (filtertype == 1) then
  print*,"Boxcar Filtering"
  write(102,*) "Boxcar Filtering"
  print*, "Boxcar Size (Cyclic) :", filterSize
  write(102,*) "Boxcar Size (Cyclic) :", filterSize
  print*, "Boxcar Passes        :", filterPasses
  write(102,*) "Boxcar Passes        :", filterPasses
else
  print*,"SG Filtering"
  write(102,*) "SG Filtering"
  print*, "SG Size (Cyclic) :", filterSize
  write(102,*) "SG Size (Cyclic) :", filterSize
  print*, "SG Passes        :", filterPasses
  write(102,*) "SG Passes        :", filterPasses
end if

print*,
write(102,'(/)') 
print*, "Integration Method"
write(102,*) "Integration Method"
print*, "============================="
write(102,*) "============================="

if (integrationChoice == 0) then
  print*, "Composite Newton Cotes"
  write(102,*) "Composite Newton Cotes"
else if (integrationChoice == 1) then
  print*, "Romberg"
  write(102,*) "Romberg"
else if (integrationChoice == 2) then
  print*, "Adaptive Quadurature"
  write(102,*) "Adaptive Quadurature"
else
  print*, "Gauss"
  write(102,*) "Gauss"
end if

print*,
write(102,'(/)')

print*, "Plot File Data"
write(102,*) "Plot File Data"
print*, "====================="
write(102,*) "====================="
print*, "File: ", outputFile
write(102,*) "File: ", outputFile

print*, "Plot shifted ", plotShift, "ppm for TMS calibration"
write(102,*) "Plot shifted ", plotShift, "ppm for TMS calibration"
print*,
write(102,'(/)')


write(*,"(5X, 7a10, 5X, /)", advance = 'no') "Peak ","Begin","End ","Location","Top","Area","Hydrogens"
write(102,"(5X, 7a10, /)", advance = 'no') "Peak ","Begin","End ","Location","Top","Area","Hydrogens"
write(*, '(5X, 6a10)'), " =======", "========= ", "=============== ", " =========== ", " =============== ", " ========="
write(102, '(5X, 6a10)'), " =======", "========= ", "============ ", " ============= ", " =============== ", " ========="

do i=1,peaks, 1
  print*, i, lowIntPts(i), upIntPts(i), peakPos(i),peakTop(i),intAreas(i),hydroArr(i)
  write(102, *) i, lowIntPts(i), upIntPts(i), peakPos(i),peakTop(i),intAreas(i),hydroArr(i)
end do

close(102)

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TMS Shift !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Finds the right most peak and shifts the
!!! data so it occurs at position 0

subroutine TMShift
        logical :: run
        integer :: p , ac, l
        double precision :: arry(0:500),arrx(0:500), highY, shift
        
        highY = 0
        run = .true.
        ac = 500
        p = count

        !! Determmines where the right most 
        !! data point crosses the beaseline
        do while (ypts(p) < baseline)      
                p = p - 1
        enddo
        

        !! Sotores data points untill they fall below baseline
         do while (run)
                if (ypts(p) > baseline) then
                        arry(ac) = ypts(p)
                        arrx(ac) = xpts(p)
                        ac = ac - 1
                else
                        run = .false.
                end if
                p = p - 1
         enddo

        !! Searches for the greatest value
        do l = ac, 500, 1
                !print *, "Arry(ac)", arry(l), l
                if (arry(l) > highY) then
                        !print *, "High ", highY
                        highY = arry(l)
                        shift = arrx(l)
                endif                
        enddo
        

        !! Shifts data set by the x-value of the highest
        !! y-value 
        do l = 0, count,1
                xpts(l) = xpts(l) - shift
        end do
        
        !! Records how big the shift was
        plotShift = shift

end subroutine TMShift





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Box Car !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Smoothes the data by avergaing data points of an input width

subroutine boxFilter
    double precision :: ysum, temp(0:MAXPTS)
    integer :: offset, i, j, l, k, lowBound, upBound

    !! Figures out how many points to have on either side
    offset = (filterSize - 1)/ 2
    
    ysum = 0.0D0


    !! Runs so that every data point is in "the center"
    iloop: do i = 0, count, 1
        ysum = 0.0D0
        lowBound = i - (offset)
        upBound = i + (offset)
        
        !! Calculates the values of the numbers
        !! that need to wrap around the bottom
        !! of the data set
        if (lowBound < 0) then
            do k = lowBound , -1, 1
                ysum = ysum + ypts(count + k + 1)
                lowBound = lowBound + 1
            end do
        end if
        

        !! Calculates the values of the numbers
        !! that need to wrap aroud the top of 
        !! the data set
        if (upBound > count) then
           
            do l = upBound - count - 1, 0, -1
                ysum = ysum + ypts(l)
                upBound = upBound - 1
            end do 
        end if


        !! Calulaets the value of the numebrs
        !! that do not need to wrap around the data set
        do j = lowBound, upBound, 1
            ysum = ysum + ypts(j)
        end do

        !! Finishes the average
        ysum = ysum / real(filterSize, kind = 8)
        

        !! Updates the data point
        ypts(i) = ysum
    end do iloop

end subroutine boxFilter





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SG Filter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Filters the data set by using a filter of input Size
!! and applying weights to vlalues to achieve a more accurate
!! result

subroutine SGFilter

    double precision :: five(1:5), eleven(1:11), seventeen(1:17), weights(1:17), wfive,weleven,wseventeen,div
    double precision :: ysum, temp(0:MAXPTS)
    integer :: offset, i, j, l, k, lowBound, upBound, wcount

    !! Weight values
    five(1:5) = (/-3,12,17,12,-3/)
    wfive = 35
    eleven(1:11) = (/-36,9,44,69,84,89,84,69,44,9,-36/)
    weleven = 429
    seventeen(1:17) = (/-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21/)
    wseventen = 323
    
    offset = (filterSize - 1)/ 2
    wcount = 1    
    ysum = 0.0D0
    
    !! Figures out what weight set to use
    if (filterSize == 5) then
      weights(1:5) = five(1:5)
      div = wfive
    else if (filterSize == 11) then
      weights(1:11) = eleven(1:11)
      div = weleven
    else
      weights(1:17) = seventeen(1:17)
      div = wseventeen
    end if


    iloop: do i = 0, count, 1
        ysum = 0
        wcount = 1
        lowBound = i - offset
        upBound = i + offset

        !! Calculates the values of the numbers
        !! that need to wrap around the bottom
        !! of the data set
        if (lowBound < 0) then
            do k = lowBound , -1, 1
                ysum = ysum + ( ypts(count + k + 1) * weights(wcount) )
                lowBound = lowBound + 1
                wcount = wcount + 1
            end do
        end if
        

        !! Calculates the values of the numbers
        !! that need to wrap around the top
        !! of the data set
        if (upBound > count) then
            do l = upBound - count - 1, 0, -1
                ysum = ysum + ( ypts(l) * weights(wcount) )
                upBound = upBound - 1
                wcount = wcount + 1
            end do 
            wcount = 1
        end if

        !! Calculates the values of the numbers
        !! that do not need to wrap around the
        !! data set
        do j = lowBound, upBound, 1
            ysum = ysum + ( ypts(j) * weights(wcount) )
            wcount = wcount + 1
        end do


        !! Divides to finish the average
        ysum = ysum / div
        
        !! Updates the value at that point
        ypts(i) = ysum
    end do iloop

end subroutine SGFilter



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Build Spline !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Takes in the data set and constructs the weights for a natural cubic spline
!! These values will then be used later to find the values of points

subroutine buildSpline(xpts, ypts, bvals, cvals, dvals, count)


double precision, intent (out) :: bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)
double precision, intent (in) :: xpts(0:MAXPTS), ypts(0:MAXPTS)
double precision :: alvals(0:MAXPTS), lvals(0:MAXPTS)
double precision :: muvals(0:MAXPTS),zvals(0:MAXPTS)
double precision :: hvals(0:MAXPTS)
integer, intent (in) :: count
integer :: i, j, l



!! The next few portions of code
!! calculate intermediary variables
!! used to calulate the values to
!! build the actual spline

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


!! Builds arrays with the values used to calculate 
do j = count, 0, -1
    cvals(j) = zvals(j)- muvals(j)*cvals(j+1)
    bvals(j)= (ypts(j+1)-ypts(j))/hvals(j) - hvals(j) * (cvals(j+1) + 2*cvals(j))/3
    dvals(j) = (cvals(j+1) - cvals(j)) / (3 * hvals(j))
end do

end subroutine buildSpline



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Find int points !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Seaches through the data and finds pairs of points on oppposite sides of the
!! baseline and finds the intersection by using a bisection algorithm
!! Then it stores the points as either the begining or end of a peak

subroutine findIntPoints
integer :: peakCount, i

peakCount = 1
i = 1
do while (xpts(i) < 0.0D0) !! Run while TMS is not reached
  if (ypts(i-1)>baseline .and. ypts(i)<baseline) then !! IF end of peak
    upIntPts(peakCount) = bisection(xpts(i-1),xpts(i))
    peakCount = peakCount + 1
  else if(ypts(i-1)<baseline .and. ypts(i)>baseline) then !! IF begining of peak
    lowIntPts(peakCount) = bisection(xpts(i-1),xpts(i))
  end if
  
  i = i + 1
end do

peaks = peakCount - 1

end subroutine findIntPoints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Integrate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Runs the correct integration based off of user choice
!! Uses the previously stored peak begin and end positions
!! to integrate all areas
!! Stores areas in an array

subroutine integrate
integer :: k

do k = 1, peaks, 1
  if (integrationChoice == 0) then
    intAreas(k) = CNC(lowIntPts(k), upIntPts(k))
  else if (integrationChoice == 1) then
    intAreas(k) = romberg(lowIntPts(k), upIntPts(k))
  else if (integrationChoice == 2) then
    intAreas(k) = adquad(lowIntPts(k), upIntPts(k), tolerance)
  else
    intAreas(k) = Gauss(lowIntPts(k),upIntPts(k))
  end if
end do

end subroutine integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!! Find Tops !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Uses the data on the begining and end of peaks
!! to caluclate the height at their midpoint
!! (adjusts for baseline)
!! stores both the location of the midpoitn and its
!! adjusted height

!! Midpoint = (upper bound + lower bound) / 2

subroutine findTops
integer :: i
  do i=1, peaks, 1
    peakTop(i) = calcSpline( ( lowIntPts(i) + upIntPts(i) ) / 2.0D0)
    peakPos(i) = ( lowIntPts(i) + upIntPts(i) ) / 2.0D0
  end do
end subroutine findTops


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calc Hydrogens !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Caluclates the number of hydrogens in an area
!! by determining the the smallest area and dividing
!! all the areas by it
!! stores hydrogen values in an array

subroutine calcHydrogens
  integer :: i
  double precision :: minArea
  minArea = intAreas(1)
  
  !! Determines the smallest area
  do i=2, peaks
    if (intAreas(i) < minArea) then
      minArea = intAreas(i)
    end if
  end do

  !! Divides all areas by the smallest area
  do i=1, peaks, 1
    hydroArr(i) = IDINT(intAreas(i) / minArea)
  end do

end subroutine calcHydrogens


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Functions !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calc Spline !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Determines which of the previously calculated weights for a spline to 
!! use to calulate the y-value at an x-location

!! Takes the x-location as input
!! Returns the y-value


double precision function calcSpline(find)!, xpts, ypts, bvals, cvals, dvals, count)

double precision, intent (in) :: find
double precision :: ans
integer :: eqchoose, i

!! Figures out what 'equation' (set of weights) to use
eqchoose = -1
do i = 0, count
    if (find >= xpts(i)) eqchoose = eqchoose + 1
end do

!! Calulates the value at that lcation using the correct weights
ans = ypts(eqchoose) + bvals(eqchoose) * (find - xpts(eqchoose))
      ans = ans + cvals(eqchoose) * ((find - xpts(eqchoose))**2)
      ans = ans + dvals(eqchoose) * ((find - xpts(eqchoose))**3)

calcSpline = ans

end function calcSpline


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bisection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Determines the location of a 'zero' with the baseline given
!! two points: one above and one below the baseline
!! is run to user set tolerance

!! takes in the lower and upper values of the peak
!! returns the x-location of the 'zero'

double precision function bisection (lowerBound, upperBound)

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


!! Runs untill tolerance is met
do while(difference > tolerance)

    !! replaces the bound which is 
    !! on the same side of the baseline as
    !! the 'guess' value with the 'guess'
    if (lval > baseline .and. answer > baseline) then
        tLower = guess
    else if (lval < baseline .and. answer < baseline) then
        tLower = guess
    else
        tUpper = guess
    end if

    !! Cals new guess for next run
    guess = (tLower + tUpper ) / 2.0D0

    !! Sets new upper and lower values for next run
    lval = calcSpline(tLower)
    uval = calcSpline(tUpper)

    !! Holds onto answer and calculates difference
    !! between current and prevoid answers
    answer = calcSpline(guess)
    difference = abs(prevAnswer - answer)
    prevAnswer = answer
end do
    
  bisection = guess

end function bisection




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simpson Integration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! A basic integration using Simpsons' method
!! takes in the lower and upper bound of the peak
!! returns the area
!! adjusts for baseline

double precision function simpson(low, up)
double precision, intent (in):: up, low
double precision :: h

h = (up - low) / 2.0D0

simpson = (h/3.0D0) * ((calcSpline(low) - baseline) + (4.0D0 * (calcSpline(low+h)-baseline) ) + (calcSpline(up)-baseline))

end function simpson


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compotite Newton Cotes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Uses the Composite Newton Cotes method to integrate the peak
!! by addinng more subsections to calcualte untill tolerance is met

!! takes in the lower and upper bound for inetegration
!! returns the area

double precision function CNC(low, up)

double precision, intent (in):: low,up
double precision :: h, prev, cur, ends, evens, odds, ihold, tol
integer :: run, counter, i

tol = tolerance
counter = 0

cur = 0.0D0
prev = 1.0D0
ends = calcSpline(low) + calcSpline(up) - 2*baseline
run = 2


!! As the loop increments, more sub sections are added to integrate
do while (abs(cur - prev) >= tol)
  prev = cur
  evens = 0.0D0
  odds = 0.0D0
  h = (up - low) / run
  
  !! Calcualtes the values of the subsections
  do i = 1, run - 1
    ihold = i
    !! Add to evens
    if (mod(i,2)==0) then
      evens = evens + ( calcspline(low + (ihold*h)) - baseline )
    !! Add to odds
    else
      odds = odds + ( calcSpline(low + (ihold*h)) - baseline)
    end if
  end do
  
  !! Calcuates the Composite Newton Cotes integration
  cur = (h* (ends + (2.0D0 * evens) + (4.0D0* odds)) )/3.0D0
  run = run + 2
end do

CNC = cur
end function CNC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Romberg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Integrates using Romberg integration
!! Builds an expanding list of values that determines
!! an accurate integral

!! takes in the lower and upper bounds
!! retuns the area

double precision function romberg(low, up)


double precision, intent (in) :: low, up
double precision :: intsum, ans
double precision :: h, arr(1:2,1:10)
integer :: i, k
logical :: run

ans = -1.0D0
run = .true.
h = up - low

!! Calcs the starting point
arr(1,1) =(h/2) * (calcSpline(low) + calcSpline(up) - 2*baseline)

i = 2
do while (i <= 10 .and. run )
    intsum = 0.0D0
    
    !! Calcs the 'left most' value for the row
    do k = 1, 2**(i-2)
        intsum = intsum + ( calcSpline(low + (k-0.5D0)*h) - baseline)
    end do
    
    !! Recursively calcs the rest of the row
    arr(2,1) = (0.5D0) * ( arr(1,1) + h * intsum)
    do j = 2, i
        arr(2,j) = arr(2,j-1) + ( (arr(2,j-1) - arr(1,j-1)) / (4**(j-1) -1) )
    end do

    !! Checks tolerance
    if (abs(arr(1,i-1) - arr(2,i) ) < tolerance) run = .false.

    !! upadates storage arays
    do l=1,10
        arr(1,l) = arr(2,l)
    end do

    romberg = arr(1,i)

    h = h / 2
    i = i + 1
end do

end function romberg





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Adaptive Quadurature !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! Uses the Adaptive Quadurature Method of integration
!! This method relys on recursively calling itelf using continuously smaller
!! bounds

!! takes in the lower and upper bounds as well as the tolerance
!! retuns the area of the current bounds

double precision recursive function adquad(low, up, tol) result (ans)
double precision, intent (in) :: low,up,tol
double precision :: mid, cur, halves

!! Calcs the areas using Simpson's itegration method
mid =( up + low ) / 2
cur = simpson(low,up)
halves = simpson(low,mid) + simpson(mid,up)

!! If tolerance is met
if (cur - halves < tol) then !! return the value betwene these bounds
  ans= halves
else !! Recursive call with smaller bounds
  ans = adquad(low, mid, tol/2) + adquad(mid, up, tol/2)
end if


end function adquad



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gauss Quad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Integrates using Gaussian Quadurature which
!! calulates by transforming the function onto the [-1,1]
!! domain and using weights and new root (x) values

!! takes lower and upper bounds
!! returns area

double precision function Gauss(low,up)
  double precision, intent (in) :: low,up
  double precision :: shift, cur, prev, root, weight
  character (len = 9) :: fileArr(1:9)
  integer :: fileChoice

  !! Array of files holding the appropritae weights and roots
  fileArr(1:) = (/'002pt.dat','004pt.dat','008pt.dat','016pt.dat','032pt.dat','064pt.dat','128pt.dat', '256pt.dat','512pt.dat'/)
  
  cur = 0.0D0
  prev = 1.0D0
  
  fileChoice = 1
  
  !! Check tolerance
  do while (abs(cur - prev)>=tolerance .and. fileChoice < 10) 
    prev = cur
    cur = 0.0D0

    !! Reading file acts as a de-facto loop
    open (unit = 100, file = fileArr(fileChoice), status = 'old') ! Open current file
  101    read(100,*, end = 102) root, weight ! Read in the roots and weights
            !! Calcuates the apropriate shift and applys the new root
            !! then sums the values of the number of appropritae
            !! calcuations 
            shift = ( ( (up-low)*root ) + (low + up) )/ 2.0D0
            cur = cur + ( ( calcSpline(shift) - baseline )* weight)
         goto 101
  102 continue
    close(100)
    cur = cur * ( (up-low)/2.0D0 )
    fileChoice = fileChoice + 1
  end do

  Gauss = cur

end function Gauss


end program nmr
