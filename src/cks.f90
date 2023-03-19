subroutine coefs(p,phi,q,theta,d,m,cks)
  !----------------------------------------
  ! Calculate the coefficients of 
  ! lambda(z) = [theta(z)/phi(z)](1-z)^{-d}
  ! also works for d = 0.
  !--------------------------------------------------
  ! theta = (/1,theta(1),theta(2),...,theta(q)/)
  ! phi = (/1,-phi(1),...,-phi(p)/)
  !--------------------------------------------------
  implicit none
  integer :: m,p,q
  real (8) :: d
  real (8) :: theta(0:q), phi(0:p)
  real (8), dimension(0:m) :: cks
  real (8), dimension(0:m) :: delta, tau
  integer:: i,j

  cks(:) = 0; cks(0) = 1
	
  ! ARMA(0,0) => lambda(z) = 1
  if(p==0 .and. q==0. .and. d==0) return
  
	! ARFIMA(0,d,0) => lambda(z) = (1-z)^{-d}
  if(p .eq. 0 .and. q.eq. 0) then
     do i = 1, m
        cks(i) = ((i-1+d)/(i))*cks(i-1)
     end do
     return
  end if
	
  ! IF THE PROGRAM GETS HERE IT MEANS THAT EITHER WE HAVE AN
  ! ARMA(p,q), OR AN ARFIMA(p,d,q)
  ! tau(z) = phi(z)*delta(z)
  tau(:) = 0
  if(d == 0) then
     ! ARMA(p,q); tau(z) = phi(z)
     tau(0:p) = phi
  else
     ! ARFIMA(p,d,q)
     delta(:) = 0; delta(0) = 1
     do i = 1, m
        delta(i) = ((i-1-d)/(i))*delta(i-1)
     end do
     if(p == 0) then
        ! ARFIMA(0,d,q); tau(z) = delta(z)
        tau = delta
     else
        !  ARFIMA(p,d,q)_s; tau(z) = phi(z)*delta(z)
        tau(:) = 0
        do i = 0,m
           if (i < p) then
              do j = 0,i
                 tau(i) = tau(i) + delta(i-j)*phi(j)
              end do
           else
              do j = 0,p
                 tau(i) = tau(i) + delta(i-j)*phi(j)
              end do
           end if
        end do
     end if
  end if
  ! lambda(z) = theta(z)/phi(z)*(1-z)^{-d}
  ! THEREFORE, theta(z) = lambda(z)*tau(z) and we have the recursion formula...
  do i = 1,m
     do j = 0,(i-1)
        cks(i)=cks(i)+cks(j)*tau(i-j)
     end do
     if (i <= q) then
        cks(i) = theta(i) - cks(i)
     else
        cks(i) = - cks(i)
     end if
  end do
  return
end subroutine coefs
