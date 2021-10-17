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
double precision :: bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)




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


if (filtertype /= 0 .and. filterSize /= 0) then
    if (filtertype == 1) then
        do l = 1, filterPasses
            print *, "pass ", l
            call boxFilter(filterSize, xpts, ypts, count)
        end do
    else if (filtertype == 2) then
    end if
else
    print *, "Filter is off"
end if


print *, ypts(0)

call buildSpline(xpts, ypts, bvals, cvals, dvals, count)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine test(ypts)
    double precision:: ypts(0:MAXPTS)
    print *, ypts(0) 
end subroutine test


!!!!!!!!!!!!! Box Car !!!!!!!!!!!!!!!!!!!

subroutine boxFilter(filterSize, xpts, ypts, count)

    integer, intent (in) :: filterSize, count
    double precision, intent (inout) :: xpts(0:MAXPTS), ypts(0:MAXPTS)
    double precision :: ysum
    integer :: offset

    offset = filterSize / 2
    do i = offset + 1, (count - offset) + 1, 1
        ysum = 0
        do j = i - offset, i + offset
            ysum = ysum + ypts(j)
        end do
        ysum = ysum / filterSize
        print *, "New Y: ", ysum, " | for ", xpts(i)
        ypts(i) = ysum
    end do

end subroutine boxFilter

!!!!!!!!!! SG Filter !!!!!!!!!!!!!!

subroutine SGFilter(filterSize, count, ypts, m)

integer, intent (in) :: count
double precision, intent (inout) :: ypts(:)
double precision, intent(out) :: m

m = count - (filterSize - 1)




end subroutine SGFilter


!!!!!!!!!! Build Spline !!!!!!!!!!!!!

subroutine buildSpline(xpts, ypts, bvals, cvals, dvals, count)


double precision, intent (out) :: bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)
double precision, intent (in) :: xpts(0:MAXPTS), ypts(0:MAXPTS)
double precision :: alvals(0:MAXPTS), lvals(0:MAXPTS)
double precision :: muvals(0:MAXPTS),zvals(0:MAXPTS)
double precision :: hvals(0:MAXPTS)

integer :: i, count

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

do l = 0, count - 1
    print *, ypts(l), bvals(l), cvals(l)
end do

end subroutine buildSpline





!!!!!!!!!! Calc Spline !!!!!!!!!!!!!!!!!

double precision function calcSpline(find, xpts, ypts, bvals, cvals, dvals, count)


double precision, intent (in) :: xpts(0:MAXPTS),ypts(0:MAXPTS), bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)
double precision, intent (in) :: find
double precision :: ans
integer :: eqchoose, count

eqchoose = -1
do i = 0, count
    if (find >= xpts(i)) eqchoose = eqchoose + 1
end do

ans = ypts(eqchoose) + bvals(eqchoose) * (find -xpts(eqchoose))
      ans = ans + cvals(eqchoose) * ((find -xpts(eqchoose))**2)
      ans = ans + dvals(eqchoose) * ((find -xpts(eqchoose))**3)

calcSpline = ans

end function calcSpline


!!!!!!!!!!!!!! Bisection !!!!!!!!!!!!!!!

double precision function bisection (lowerBound, upperBound, xpts, ypts, tolerance, baseline)

double precision, intent (in) :: xpts(0:MAXPTS),ypts(0:MAXPTS)
double precision, intent (in) :: lowerBound, upperBound, baseline
double precision :: tLower, tUpper, guess, prevAnswer, answer, difference, calcSpline, lval, uval
integer :: i

tLower = lowerBound
tUpper = upperBound
guess = (upperBound + lowerBound ) / 2
lval = calcSpline(tLower, xpts, ypts, bvals, cvals, dvals, count)
uval = calcSpline(tUpper, xpts, ypts, bvals, cvals, dvals, count)
prevAnswer = 0
answer = calcSpline(guess, xpts, ypts, bvals, cvals, dvals, count)
difference = 0

do while(difference > tolerance)
    if (lval > baseline .and. answer > baseline) then
        tLower = guess
        guess = (tLower + tUpper ) / 2
    else if (lval < baseline .and. answer < baseline) then
        tLower = guess
        guess = (tLower + tUpper ) / 2
    else
        tUpper = guess
        guess = (tLower + tUpper ) / 2
    end if

    lval = calcSpline(tLower, xpts, ypts, bvals, cvals, dvals, count)
    uval = calcSpline(tUpper, xpts, ypts, bvals, cvals, dvals, count)
    
    answer = calcSpline(guess, xpts, ypts, bvals, cvals, dvals, count)
    difference = abs(prevAnswer - answer)
    prevAnswer = answer
end do
    
    bisection = answer

end function bisection

end program nmr