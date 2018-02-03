! fortran subroutines for bottlenecks in the colloid simulation code.
! Created by Joshua Larsen for LB-Colloids as part of a PhD dissertation
! from the University of Arizona Department of Soil Water and Environmental Sciences

subroutine colcolarray(c_arr, colloids, fxlen, fylen, cxlen, cylen, ccenter, collen, f_arr)
    ! subroutine creates the colloid colloid interaction array for DLVO
    integer, intent(in)                            :: fxlen, fylen, cxlen, cylen, ccenter, collen
    real*8, dimension(5, 5), intent(in)            :: c_arr
    integer, dimension(collen, 2), intent(in)      :: colloids
    real*8, dimension(fylen, fxlen), intent(out)   :: f_arr
    integer                                        :: c, x, y, xi, yi, clx, crx, flx, frx, xclx, ycly, cty, cby, fty, fby, t
    logical                                        :: verbose

    ! todo: Test and fuss with the fortran indexing calculation to get proper colloid force set.
    f_arr(:, :) = 0.  ! set initial fortran array to zero value
    verbose = .FALSE.

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
                print *, flx ! this stupid shit, figure out how to remove me!

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
end subroutine colcolarray