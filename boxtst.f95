program box

parameter (MAXPTS = 10000)

double precision :: xtmp,ytmp, xpts(0:MAXPTS), ypts(0:MAXPTS)
integer :: count, filterSize, runs


filterSize = 11
runs = 3

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


do i=1,3
  call boxFilter
end do


open (unit=17, file="box.dat", status='unknown')

  do i=0, count
    write (17,*) xpts(i), ypts(i)
  end do
close(17)



contains



subroutine boxFilter

    double precision :: ysum, temp(0:MAXPTS)
    integer :: offset, i, j, l, k, lowBound, upBound

    temp(0:) = ypts(0:)
    offset = (filterSize - 1)/ 2
    
    ysum = 0.0D0

    iloop: do i = 0, count, 1
      print *, "Point :  " , i
        ysum = 0
        lowBound = i - offset
        upBound = i + offset
        
        if (lowBound < 0) then
            do k = lowBound , -1, 1
              print *, "  Include below ", count + k + 1, "  value  ", ypts(count + k + 1)
                ysum = ysum + ypts(count + k + 1)
                lowBound = lowBound + 1
            end do
        end if
        
        if (upBound > count) then
            do l = upBound - count - 1, 0, -1
              print *, "  Include above ", l, "  value  ", ypts(l)
                ysum = ysum + ypts(l)
                upBound = upBound - 1
            end do 
        end if


        do j = lowBound, upBound, 1
            print *, "  Include in ", j, "  value  ", ypts(j)
            ysum = ysum + ypts(j)
        end do

        print *, "  Before divide :  ", ysum
        ysum = ysum / real(filterSize, kind = 8)
        print *, "  Sum : ", ysum
        print *, ""
        
        !temp(i) = ysum
        ypts(i) = ysum
    end do iloop

    !ypts(0:) = temp(0:)
end subroutine boxFilter


end program box