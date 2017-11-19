! Fortran subroutines to be compiled with f2py for fast and efficient LB2D code
! that can be called from a python wrapper. This is my first attempt with fortran
! Created by Joshua Larsen


subroutine f_rho(f, ylen, xlen, rho)
  ! function calculates density at each node from the LB-distribution function
  integer, intent(in)             				:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in)  	:: f
  real, dimension(ylen, xlen), intent(out)     	:: rho
		
  rho = SUM(f, DIM=1)
end subroutine f_rho
		

subroutine f_u(f, rho, ylen, xlen, uy, ux)
  ! function calculates velocity arrays where up and right are positive
  integer, intent(in)							:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in)	:: f
  real, dimension(ylen, xlen), intent(in)    	:: rho
  real, dimension(ylen, xlen), intent(out)      :: uy, ux
  
  uy = (SUM(f((/2,3,4/),:,:), DIM=1) - SUM(f((/6,7,8/),:,:), DIM=1))/rho
  ! ux = (SUM(f((/1,2,8/),:,:), DIM=1) - SUM(f((/4,5,6/),:,:), DIM=1))/rho
  ux = (SUM(f((/4,5,6/),:,:), DIM=1) - SUM(f((/1,2,8/),:,:), DIM=1))/rho
end subroutine f_u

		
subroutine f_usqr(uy, ux, ylen, xlen, usqr)
  ! function to calculate velocity squared
  integer, intent(in)							:: ylen, xlen
  real, dimension(ylen, xlen), intent(in)  		:: uy, ux
  real, dimension(ylen, xlen), intent(out)	    :: usqr
	
  usqr = uy*uy + ux*ux
end subroutine f_usqr

subroutine f_eu(uy, ux, tau, g, ylen, xlen, eu)
  ! function to apply eigenvectors to LB velocities in calculating eq. distribution
  integer, intent(in)							:: ylen, xlen
  real, intent(in)							:: tau, g
  real, dimension(ylen, xlen), intent(in)				:: uy, ux
  real, dimension(9, ylen, xlen), intent(out)				:: eu
  real, dimension(ylen, xlen)						:: duy
  duy = uy - (tau*g)
  eu(1, :, :) = ux
  eu(2, :, :) = ux + duy
  eu(3, :, :) = duy
  eu(4, :, :) = duy - ux
  eu(5, :, :) = -ux
  eu(6, :, :) = -ux - duy
  eu(7, :, :) = -duy
  eu(8, :, :) = ux - duy
  eu(9, :, :) = ux * 0.
  
end subroutine f_eu

subroutine f_feq(eu, rho, usqr, wi, csqr, ylen, xlen, feq)
  ! function to calculate the LB equilibrium distribution function 
  integer, intent(in)							:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in)	:: eu
  real, dimension(ylen, xlen), intent(in)    	:: rho
  real, dimension(ylen, xlen), intent(in)		:: usqr
  real, intent(in)								:: csqr
  real, dimension(9), intent(in)				:: wi
  real, dimension(9, ylen, xlen), intent(out)	:: feq
  
  feq(1, :, :) = wi(1)*rho*(1 + eu(1, :, :)/csqr + 0.5*(eu(1, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(2, :, :) = wi(2)*rho*(1 + eu(2, :, :)/csqr + 0.5*(eu(2, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(3, :, :) = wi(3)*rho*(1 + eu(3, :, :)/csqr + 0.5*(eu(3, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(4, :, :) = wi(4)*rho*(1 + eu(4, :, :)/csqr + 0.5*(eu(4, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(5, :, :) = wi(5)*rho*(1 + eu(5, :, :)/csqr + 0.5*(eu(5, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(6, :, :) = wi(6)*rho*(1 + eu(6, :, :)/csqr + 0.5*(eu(6, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(7, :, :) = wi(7)*rho*(1 + eu(7, :, :)/csqr + 0.5*(eu(7, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(8, :, :) = wi(8)*rho*(1 + eu(8, :, :)/csqr + 0.5*(eu(8, :, :)/csqr)**2 - usqr/(2*csqr))
  feq(9, :, :) = wi(9)*rho*(1 + eu(9, :, :)/csqr + 0.5*(eu(9, :, :)/csqr)**2 - usqr/(2*csqr))
end subroutine f_feq

subroutine f_collision(f, feq, tau, ylen, xlen, fcol)
  !function to calculate the LB collosion operator
  integer, intent(in)							:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in)	:: f, feq
  real, intent(in)								:: tau
  real, dimension(9, ylen, xlen), intent(out)	:: fcol
  
  fcol = f - ((f - feq)/tau)
end subroutine f_collision

subroutine f_zhohe(f, rho, ylen, xlen, fzhe)
  ! function to calculate the zho-he pressure boundary conditions of LB
  integer, intent(in)							:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in) 	:: f
  real, dimension(ylen, xlen), intent(in)		:: rho
  real, dimension(9, ylen, xlen), intent(out)	:: fzhe
  real, dimension(ylen, xlen)					:: vy_lb, vy_ub, f2, f3, f4, f6, f7, f8
  integer										:: j, k
  
  vy_lb = 1 - (f(9, :, :) + f(1, :, :) + f(5, :, :) + 2*(f(6, :, :) + f(7, :, :) + f(8, :, :))/rho)
  
  vy_ub = -(1 - (f(9, :, :) + f(1, :, :) + f(5, :, :) + 2*(f(2, :, :) + f(3, :, :) + f(4, :, :))/rho))
  
  ! compute the unkown distribution values on the lower boundary
  f2 = (1/6) * rho * vy_lb + (f(5, :, :) - f(1, :, :))/2 + f(6, :, :)
  f3 = (2/3) * rho * vy_lb + f(7, :, :)
  f4 = (1/6) * rho * vy_lb + (f(1, :, :) - f(5, :, :))/2 + f(8, :, :)
  
  ! compute the unkown distribution values on the upper boundary
  f6 = -(1/6) * rho * vy_ub + (f(1, :, :) - f(5, :, :))/2 + f(2, :, :)
  f7 = -(2/3) * rho * vy_ub + f(3, :, :)
  f8 = -(1/6) * rho * vy_ub + (f(5, :, :) - f(1, :, :))/2 + f(4, :, :)
  
  !$OMP PARALLEL DO
  do k=1, xlen
    do j=1, ylen
      
      if (j.eq.ylen) then
        fzhe(2, j, k) = f(2, j, k)
        fzhe(3, j, k) = f(3, j, k)
        fzhe(4, j, k) = f(4, j, k)
        fzhe(6, j, k) = f6(j, k)
        fzhe(7, j, k) = f7(j, k)
        fzhe(8, j, k) = f8(j, k)
        fzhe(1, j, k) = f(1, j, k)
        fzhe(5, j, k) = f(5, j, k)
        fzhe(9, j, k) = f(9, j, k)
      else if (j.eq.1) then
        ! flipped 1 and ylen
        fzhe(2, j, k) = f2(j, k)
        fzhe(3, j, k) = f3(j, k)
        fzhe(4, j, k) = f4(j, k)
        fzhe(6, j, k) = f(6, j, k)
        fzhe(7, j, k) = f(7, j, k)
        fzhe(8, j, k) = f(8, j, k)
        fzhe(1, j, k) = f(1, j, k)
        fzhe(5, j, k) = f(5, j, k)
        fzhe(9, j, k) = f(9, j, k) 
      else
        fzhe(2, j, k) = f(2, j, k)
        fzhe(3, j, k) = f(3, j, k)
        fzhe(4, j, k) = f(4, j, k)
        fzhe(6, j, k) = f(6, j, k)
        fzhe(7, j, k) = f(7, j, k)
        fzhe(8, j, k) = f(8, j, k)
        fzhe(1, j, k) = f(1, j, k)
        fzhe(5, j, k) = f(5, j, k)
        fzhe(9, j, k) = f(9, j, k)
      endif

    enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine f_zhohe

subroutine f_bounceback(f, fcol, image, ylen, xlen, fbounce)
  ! function to apply lattice boltzmann bounceback/no slip boundaries
  integer, intent(in)							:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in)	:: f
  real, dimension(9, ylen, xlen), intent(in) 	:: fcol
  logical, dimension(ylen, xlen), intent(in)	:: image
  real, dimension(9, ylen, xlen), intent(out)	:: fbounce
  integer 										:: j, k
  
  !$OMP PARALLEL DO
  do k=1, xlen
    do j=1, ylen
    
      if (image(j,k)) then
        fbounce(1, j, k) = f(5, j, k)
        fbounce(2, j, k) = f(6, j, k)
        fbounce(3, j, k) = f(7, j, k)
        fbounce(4, j, k) = f(8, j, k)
        fbounce(5, j, k) = f(1, j, k)
        fbounce(6, j, k) = f(2, j, k)
        fbounce(7, j, k) = f(3, j, k)
        fbounce(8, j, k) = f(4, j, k)
        fbounce(9, j, k) = f(9, j, k)
      else
        fbounce(1, j, k) = fcol(1, j, k)
        fbounce(2, j, k) = fcol(2, j, k)
        fbounce(3, j, k) = fcol(3, j, k)
        fbounce(4, j, k) = fcol(4, j, k)
        fbounce(5, j, k) = fcol(5, j, k)
        fbounce(6, j, k) = fcol(6, j, k)
        fbounce(7, j, k) = fcol(7, j, k)
        fbounce(8, j, k) = fcol(8, j, k)
        fbounce(9, j, k) = fcol(9, j, k)
        
      endif
    enddo
  enddo
  !$OMP END PARALLEL DO  

end subroutine f_bounceback

subroutine f_streaming(fcol, ylen, xlen, fstream)
  ! function to complete the LB streaming step 
  integer, intent(in)							:: ylen, xlen
  real, dimension(9, ylen, xlen), intent(in) 	:: fcol
  real, dimension(9, ylen, xlen), intent(out)   :: fstream
  integer										:: j, k, jp, jn, kp, kn
  
  !$OMP PARALLEL DO
  do k=1, xlen
    do j=1, ylen
      
      if (k.eq.1) then
        kp = k + 1
        kn = xlen
      else if (k.eq.xlen) then
        kp = 1
        kn = k - 1
      else
        kp = k + 1
        kn = k - 1
      endif
      
      if (j.eq.1) then
        jp = j + 1
        jn = ylen
      else if (j.eq.ylen) then 
        jp = 1
        jn = j - 1
      else
        jp = j + 1
        jn = j - 1
      endif
      
      fstream(1, j, kp)  = fcol(1, j, k)
      fstream(2, jp, kp) = fcol(2, j, k)
      fstream(3, jp, k)  = fcol(3, j, k)
      fstream(4, jp, kn) = fcol(4, j, k)
      fstream(5, j, kn)  = fcol(5, j, k)
      fstream(6, jn, kn) = fcol(6, j, k)
      fstream(7, jn, k)  = fcol(7, j, k)
      fstream(8, jn, kp) = fcol(8, j, k)
      fstream(9, j, k)   = fcol(9, j, k)
    enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine f_streaming
