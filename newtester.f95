program fuckshit

parameter (MAXPTS = 10000)


!! Basic vars
double precision :: xpts(0:MAXPTS), ypts(0:MAXPTS), xtmp, ytmp
integer :: count, range, LWORK

integer, allocatable :: IPIV(:)

!! Complex vars
double complex, allocatable :: Zmat(:,:), Gmat(:,:)
double complex, allocatable :: Yvec(:), Cvec(:), Chold(:), WORK(:)
double complex :: sum

call readIn()
print*, "Count :", count
range = count + 1


open(unit = 16, file = "ordered.dat", status = 'unknown')
    do i = 0, count
        write (16, *) xpts(i), ypts(i)
    end do
close(16)




print*, xpts(0), ypts(0)

allocate(Yvec(1:count+1))


do i = 1, range, 1
    Yvec(i) = dcmplx(ypts(i - 1), 0.0D0)
end do

!print*, "Y vec ", Yvec(range)


allocate(Zmat(1:range, 1:range))
call buildZmat()

allocate(Gmat(1:range, 1:range))
call buildGmat()


allocate(Cvec(1:range))
allocate(Chold(1:range))
allocate(IPIV(1:range))
allocate(WORK(1:range))


print*, "Full Zmat"
!do i = 1, range, 1
!    print*, Zmat(i,1:)
!end do




print*,""
print*,"--------------------------"
print*,""


!! Applying filters

print*, "Yvec: "
!print*, Yvec(1:)


!! Calc F Coeffs

do i = 1, range, 1
    sum = (0.0D0,0.0D0)
    do j = 1, range, 1
        sum = sum + ( Zmat(i,j) * Yvec(j) )
    end do
    Cvec(i) = sum
end do

print*, "First C vec: ", Cvec(2)



print*, "Full Gmat"
!do i = 1, range, 1
!    print*, Gmat(i,1:)
!end do



!! Apply filter

do i = 1, range, 1
    sum = (0.0D0,0.0D0)
    do j = 1, range, 1
        sum = sum + ( Gmat(i,j) * Cvec(j) )
    end do
    Chold(i) = sum
end do

Cvec(1:) = Chold(1:)


print*, "Cvec: "
!print*, Cvec(1:)

!! Invert Zmat

print*, "Before invert"
print*, "Range :" , range, count

print *, "Before Zmat: ", Zmat(1,1)

do i = 1, range, 1
    do j = i, range, 1
        Zmat(i,j) = dconjg(Zmat(i,j))
    end do
end do

print*, "after Zmat: ", Zmat(1,1)
print*, "After invert"




!! Reverse filter to get new points
print*, "Reverse Filter"
do i = 1, range, 1
    sum = (0.0D0,0.0D0)
    do j = 1, range, 1
        sum = sum + ( Zmat(i,j) * Cvec(j) )
    end do
    Yvec(i) = sum
end do

print*, "Yvec: "
!print*, Yvec(1:)

print*, "Before Y", ypts(2)
print*,"Replacing ypts"
do i = 0, count
    ypts(i) = Real(Yvec(i+1), kind = 8)
end do
print*, "After Y", ypts(2)


call write()

print*, "Done"


contains


!!!!!!!!!!!!!!!! build Zmat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine buildZmat()

double complex :: DFTBase
double precision :: Pi, dcount
integer :: range, i, j

dcount = count + 1
range = count + 1
Pi = 4.0D0 * DTAN(1.0D0)
DFTBase = dcmplx(dcos( ((2.0D0*Pi)/dcount)), -1.0D0 * dsin((2.0D0*Pi)/dcount))


print*, "dcount :", dcount
print*, "dsqrt :", dsqrt(dcount)
print*, "DFTBase :", DFTBase

do i = 1, range, 1
    do j = 1, range, 1

Zmat(i,j) = ( DFTBase ** ( (i-1)*(j-1) ) ) / dsqrt(dcount)

    end do
end do

print*, "Zmat spot: ", Zmat(1,1)

end subroutine buildZmat


!!!!!!!!!!!!!!!!!!! Build Gmat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine buildGmat()
integer :: i,j, range
double precision :: dcount, e

dcount = count + 1
range = count + 1
e = 2.718281828

print*, "LN2 :", dlog(2.0D0)
do i = 1, range, 1
    do j = 1, range, 1
        if (i == j) then
            Gmat(i,j) = e ** ( (-4.0D0 * dlog(2.0D0) *real(i*j, kind = 8)) / (dcount**(3.0D0/2.0D0))  )
        else
            Gmat(i,j) = 0.0D0
        end if
    end do
end do

print*, "Gmat : ", Gmat(100,100)


end subroutine buildGmat


!!!!!!!!!!!!!!!!!!!! Read In !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readIn()

count = -1
open (unit = 16, file = "test.dat", status = 'old')
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

end subroutine readIn




subroutine write()

open(unit = 16, file = "filt.dat", status = 'unknown')
    do i = 0, count
        write (16, *) xpts(i), ypts(i)
    end do
close(16)



end subroutine write


end program fuckshit