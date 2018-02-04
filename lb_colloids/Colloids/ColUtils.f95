! fortran subroutines for bottlenecks in the colloid simulation code.
! Created by Joshua Larsen for LB-Colloids as part of a PhD dissertation
! from the University of Arizona Department of Soil Water and Environmental Sciences

subroutine colcolarray(c_arr, colloids, fxlen, fylen, cxlen, cylen, ccenter, collen, f_arr)
    ! subroutine creates the colloid colloid interaction array for DLVO
    ! not faster than numpy!!!!!!
    integer, intent(in)                            :: fxlen
    integer, intent(in)                            :: fylen
    integer, intent(in)                            :: cxlen
    integer, intent(in)                            :: cylen
    integer, intent(in)                            :: ccenter
    integer, intent(in)                            :: collen
    real*8, dimension(5, 5), intent(in)    :: c_arr
    integer, dimension(collen, 2), intent(in)      :: colloids
    real*8, dimension(fylen, fxlen), intent(out)   :: f_arr
    integer                                        :: c = 0, x = 0, y = 0, xi = 0, yi = 0
    integer                                        :: clx = 0, crx = 0, flx = 0, frx = 0
    integer                                        :: xclx = 0, ycly = 0, cty = 0, cby = 0, fty = 0, fby = 0


    ! todo: Test and fuss with the fortran indexing calculation to get proper colloid force set.
    f_arr(:, :) = 0.0  ! set initial fortran array to zero value

    do c=1, collen
        x = colloids(c, 1)
        y = colloids(c, 2)

        if (x.ne.x) then
            continue

        else if (y.ne.y) then
            continue

        else
            xi = x - (center + 1) ! Plus 1 is to facilitate the python to f95 indexing
            yi = y - (center + 1) ! Plus 1 is to facilitate the python to f95 indexing
            xclx = xi + (cxlen - 1) + 1
            ycly = yi + (cylen - 1) + 1

            if (xi.lt.1) then
                clx = -xi + 2
                crx = cxlen
                flx = 1
                frx = cxlen + xi
                ! print *, flx ! this stupid shit, figure out how to remove me!

            else if (xclx.gt.fxlen) then
                clx = 1
                crx = -(xi - fxlen) + cxlen ! Not certain that this is correct yet
                flx = xi
                frx = fxlen

            else
                clx = 1
                crx = cxlen
                flx = xi
                frx = xi + cxlen

            endif

            if (yi.lt.1) then
                cty = -yi + 2
                cby = cylen
                fty = 1
                fby = cylen + yi

            else if (ycly.gt.fylen) then
                cty = 1
                cby = -(yi - fylen) + cylen
                fty = yi
                fby = fylen

            else
                cty = 1
                cby = cylen
                fty = yi
                fby = yi + cylen

            endif

            f_arr(fty:fby, flx:frx) = f_arr(fty:fby, flx:frx) + c_arr(cty:cby, clx:crx)

        endif
    enddo

    !do i=1, fylen
    !    do j = 1, fxlen
    !        if (ISNAN(f_arr(i, j))) then
    !            f_arr(i, j) = 0

    !        else if (f_arr(i, j).gt.1.0) then
    !            f_arr(i, j) = 0

    !        else if (f_arr(i, j).lt.-1.0) then
    !            f_arr(i, j) = 0

    !        else
    !            continue

    !        endif
    !    enddo
    !enddo
end subroutine colcolarray


subroutine add_arrays(arr1, arr2, xlen, ylen, outarr)
    integer, intent(in)                        :: ylen, xlen
    real*8, dimension(ylen, xlen), intent(in)  :: arr1, arr2
    real*8, dimension(ylen, xlen), intent(out) :: outarr
    integer                                    :: i = 0, j = 0

    outarr(:, :) = 0.0

    do j = 1, xlen
        do i = 1, ylen
            outarr(i, j) = arr1(i, j) + arr2(i, j)
        enddo
    enddo

end subroutine add_arrays