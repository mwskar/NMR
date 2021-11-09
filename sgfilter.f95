program filter

parameter (MAXPTS = 10000)

double precision :: xpts(0:MAXPTS), ypts(0:MAXPTS), xtmp, ytmp
integer :: count, filterSize




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
            !print *, xtmp, ytmp
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


call sgFilter



filterSize = 5

contains

subroutine sgFilter

    double precision :: five(1:5), eleven(1:11), seventeen(1:17), weights(1:17), wfive,weleven,wseventeen,div


    double precision :: ysum, temp(0:MAXPTS)
    integer :: offset, i, j, l, k, lowBound, upBound, wcount

    filterSize = 5

    five(1:5) = (/-3,12,17,12,-3/)
    wfive = 35
    eleven(1:11) = (/-36,9,44,69,84,89,84,69,44,9,-36/)
    weleven = 429
    seventeen(1:17) = (/-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21/)
    wseventen = 323
    
    temp(0:) = ypts(0:)
    offset = (filterSize - 1)/ 2
    wcount = 1
    
    print *, "Offset: ", offset, filterSize
    
    ysum = 0.0D0
    
    if (filterSize == 5) then
      weights(1:5) = five(1:5)
      div = wfive
    else if (filterSize == 11) then
      weights(1:11) = eleven(1:11)
    else
      weights(1:17) = seventeen(1:17)
    end if

    iloop: do i = 0, count, 1
      print *, "Point :  " , i
        ysum = 0
        wcount = 1
        lowBound = i - offset
        upBound = i + offset
        print *, "Before first loop"
        if (lowBound < 0) then
            do k = lowBound , -1, 1
              !print *, "Before other print"
              print *, "  Include below ", count + k + 1, "  value  ", ypts(count + k + 1), " weight ", weights(wcount)
                ysum = ysum + ( ypts(count + k + 1) * weights(wcount) )
                lowBound = lowBound + 1
                wcount = wcount + 1
            end do
        end if
        
        if (upBound > count) then
            do l = upBound - count - 1, 0, -1
              print *, "  Include above ", l, "  value  ", ypts(l), " weight ", weights(wcount)
                ysum = ysum + ( ypts(l) * weights(wcount) )
                upBound = upBound - 1
                wcount = wcount + 1
            end do 
            wcount = 1
        end if


        do j = lowBound, upBound, 1
            print *, "  Include in ", j, "  value  ", ypts(j), " weight ", weights(wcount)
            ysum = ysum + ( ypts(j) * weights(wcount) )
            wcount = wcount + 1
        end do

        print *, "  Before divide :  ", ysum
        ysum = ysum / div
        print *, "  Sum : ", ysum
        print *, ""
        
        !temp(i) = ysum
        ypts(i) = ysum
    end do iloop

    !ypts(0:) = temp(0:)

end subroutine sgFilter

end program filter