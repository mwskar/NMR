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

double precision :: hld


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


call TMshift


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




contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine TMShift
        logical :: run
        integer :: p , ac
        double precision :: arry(15),arrx(15), highY, shift
        
        highY = 0
        run = .true.
        ac = 15
        p = count
        do while (ypts(p) < baseline)      
                p = p - 1
        enddo
        
         do while (run)
                if (ypts(p) > baseline) then
                        arry(ac) = ypts(p)
                        arrx(ac) = xpts(p)
                        ac = ac - 1
                        print *, "Ac ", ac
                else
                        run = .false.
                end if
                p = p - 1
         enddo


        do l = ac, 15, 1
                print *, "Arry(ac)", arry(l), l
                if (arry(l) > highY) then
                        print *, "High ", highY
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

        ysum = ysum / filterSize
        
        temp(i) = ysum
    end do iloop

    ypts(0:) = temp(0:)
end subroutine boxFilter





!!!!!!!!!! SG Filter !!!!!!!!!!!!!!

subroutine SGFilter

double precision :: m

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

!do l = 0, count - 1
!    print *, ypts(l), bvals(l), cvals(l)
!end do

end subroutine buildSpline



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Functions !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!! Calc Spline !!!!!!!!!!!!!!!!!

double precision function calcSpline(find)!, xpts, ypts, bvals, cvals, dvals, count)


!double precision, intent (in) :: xpts(0:MAXPTS),ypts(0:MAXPTS), bvals(0:MAXPTS), cvals(0:MAXPTS), dvals(0:MAXPTS)
double precision, intent (in) :: find
double precision :: ans
integer :: eqchoose, count

eqchoose = -1
do i = 0, count
    if (find >= xpts(i)) eqchoose = eqchoose + 1
end do

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
        guess = (tLower + tUpper ) / 2
    else if (lval < baseline .and. answer < baseline) then
        tLower = guess
        guess = (tLower + tUpper ) / 2
    else
        tUpper = guess
        guess = (tLower + tUpper ) / 2
    end if

    lval = calcSpline(tLower)
    uval = calcSpline(tUpper)

    answer = calcSpline(guess)
    difference = abs(prevAnswer - answer)
    prevAnswer = answer
end do
    
    bisection = guess

end function bisection



function ncInt(low, up)
double precision, intent (in):: up, low
double precision :: h

h = (up - low) / 2

ncInt = (h/3) * (calcSpline(low) + (4 * calcSpline(low+h) ) + calcSpline(up))

end function ncInt

function romberg(low,up)
double precision, intent (in) :: low, up
double precision :: intsum, ans
double precision :: h, arr(0:10,0:10)
integer :: i

ans = -1.0D0

i = 1

h = low - up
arr(0,0) = (h/2) * (calcSpline(low) + calcSpline(up))

do while ( (arr(1,i-1) - arr(2,i) ) < tolerance )
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

romberg = ans

end function romberg


end program nmr