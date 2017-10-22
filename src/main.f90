program Euler

! Program to solve Sod's shock tube problem for the compressible Euler Equations
! of gas dynamics in radial co-ordinates.  Solves using Godunov's method with
! Roe's Riemann Solver to calculate the numerical fluxes.  Uses the Entropy fix
! of Harten.  The method can be choosen to be either first or second order in
! space, using a variety of flux limiters.

   implicit none
!---------------------------------------------------------------------------------
   integer :: nodes, count, order, lim_choice, alpha
   double precision, dimension(:, :), allocatable :: u
   double precision :: gamma, r_min, r_max, rho_l, u_l, p_l, rho_r, u_r, p_r
   double precision :: delta_t, delta_r, r_mid, T, output_t, CFL, CFL_step
!---------------------------------------------------------------------------------

   call ReadVariables(nodes, r_min, r_max, gamma, cfl, output_t, order, &
                       lim_choice, rho_l, u_l, p_l, rho_r, u_r, p_r, r_mid, alpha)

   call ValidateVariables(nodes, r_min, r_max, r_mid, output_t, cfl, order, lim_choice)

   T = 0; count = 1
   allocate (u(-1:nodes, 1:3))
   u = 0; delta_r = 0

   call InitialConditions(u, delta_r, gamma, r_min, r_max, rho_l, u_l, p_l, rho_r, u_r, p_r, nodes, r_mid)

   do while (t < output_t)

      if (count <= 5) then
         CFL_step = 0.2*CFL
      else
         CFL_step = CFL
      endif
      count = count + 1

      call AdaptiveTimeStep(delta_t, u, delta_r, nodes, gamma, CFL_step)

      if (T + delta_t > output_T) then
         delta_t = output_T - T
      endif

      T = T + delta_t
      print *, T, delta_t

      call RiemannSolver(u, delta_r, nodes, gamma, delta_t, order, lim_choice)

      call ODESolver(u, delta_r, gamma, delta_t, nodes, alpha)

   end do

   call WriteSolution(u, delta_r, gamma, nodes)

end program Euler

subroutine InitialConditions(u, delta_r, gamma, r_min, r_max, rho_l, u_l, p_l, rho_r, u_r, p_r, nodes, r_mid)

   implicit none
!---------------------------------------------------------------------------------
! GLOBAL VARIABLES
   integer, intent(IN) :: nodes
   double precision, intent(INOUT), dimension(-1:nodes, 1:3) :: u
   double precision, intent(INOUT) :: delta_r
   double precision, intent(IN) :: gamma, r_min, r_max, r_mid
   double precision, intent(IN) :: rho_l, u_l, p_l, rho_r, u_r, p_r
!---------------------------------------------------------------------------------
! LOCAL VARIABLES
   integer :: i
!---------------------------------------------------------------------------------
   delta_r = (r_max + r_min)/nodes

   do i = -1, nodes
      if ((0.5*delta_r + i*delta_r) < r_mid) then
         u(i, 1) = rho_l
         u(i, 2) = rho_l*u_l
         u(i, 3) = rho_l*(p_l/((gamma - 1.0)*rho_l) + (u_l**2.0)/2.0)
      else
         u(i, 1) = rho_r
         u(i, 2) = rho_r*u_r
         u(i, 3) = rho_r*(p_r/((gamma - 1.0)*rho_r) + (u_r**2.0)/2.0)
      endif
   end do

   return

end subroutine InitialConditions

subroutine RiemannSolver(u, delta_r, nodes, gamma, delta_t, order, lim_choice)

   implicit none
!---------------------------------------------------------------------------------
! GLOBAL VARIABLES
   integer, intent(IN) :: nodes, order, lim_choice
   double precision, intent(INOUT), dimension(-1:nodes, 1:3) :: u
   double precision, intent(IN) :: gamma, delta_t, delta_r
!---------------------------------------------------------------------------------
! LOCAL VARIABLES
   double precision, dimension(0:nodes, 1:3, 1:3) :: K_tilde
   double precision, dimension(0:nodes, 1:3) :: lambda, F
   double precision, dimension(-1:nodes + 1, 1:3) :: alpha_tilde
   double precision, dimension(1:3) :: delta_u, limit
   double precision, external :: limiter
   double precision :: rho_l, u_l, E_l, p_l, H_l, rho_r, u_r, E_r, p_r, H_r
   double precision :: u_tilde, H_tilde, a_tilde, theta
   double precision :: a_l, a_r, rho_star_l, rho_star_r, u_star, p_star
   double precision :: a_star_l, a_star_r
   double precision, dimension(1:2) :: lambda_l, lambda_r
   double precision, external :: Eos
   integer :: i, j
!---------------------------------------------------------------------------------
   K_tilde = 0; lambda = 0; alpha_tilde = 0; F = 0; delta_u = 0
   u_tilde = 0; H_tilde = 0; rho_l = 0; E_l = 0; u_l = 0; H_l = 0; p_l = 0
   rho_r = 0; E_r = 0; u_r = 0; H_r = 0; p_r = 0; a_tilde = 0
   limit = 0; theta = 0

   do i = 0, nodes
      rho_l = u(i - 1, 1); rho_r = u(i, 1)
      u_l = u(i - 1, 2)/rho_l; u_r = u(i, 2)/rho_r
      E_l = u(i - 1, 3); E_r = u(i, 3)
      p_l = Eos(u(i-1, :), gamma); p_r = Eos(u(i, :), gamma)
      H_l = (E_l + P_l)/rho_l; H_r = (E_r + P_r)/rho_r

      delta_u(1) = u(i, 1) - u(i - 1, 1)
      delta_u(2) = u(i, 2) - u(i - 1, 2)
      delta_u(3) = u(i, 3) - u(i - 1, 3)

      u_tilde = ((rho_l**0.5)*u_l + (rho_r**0.5)*u_r)/((rho_l**0.5) + (rho_r**0.5))
      H_tilde = ((rho_l**0.5)*H_l + (rho_r**0.5)*H_r)/((rho_l**0.5) + (rho_r**0.5))
      a_tilde = ((gamma - 1.0)*(H_tilde - 0.5*(u_tilde**2)))**0.5

      lambda(i, 1) = u_tilde - a_tilde
      lambda(i, 2) = u_tilde
      lambda(i, 3) = u_tilde + a_tilde

      K_tilde(i, 1, 1) = 1.0
      K_tilde(i, 1, 2) = u_tilde - a_tilde
      K_tilde(i, 1, 3) = H_tilde - u_tilde*a_tilde

      K_tilde(i, 2, 1) = 1.0
      K_tilde(i, 2, 2) = u_tilde
      K_tilde(i, 2, 3) = 0.5*(u_tilde**2)

      K_tilde(i, 3, 1) = 1.0
      K_tilde(i, 3, 2) = u_tilde + a_tilde
      K_tilde(i, 3, 3) = H_tilde + u_tilde*a_tilde

      alpha_tilde(i, 2) = ((gamma - 1.0)/(a_tilde**2))*(delta_u(1)*(H_tilde - u_tilde**2)&
          &+ u_tilde*delta_u(2) - delta_u(3))
      alpha_tilde(i, 1) = (1.0/(2.0*a_tilde))*(delta_u(1)*(u_tilde + a_tilde) - delta_u(2) -&
          &a_tilde*alpha_tilde(i, 2))
      alpha_tilde(i, 3) = delta_u(1) - (alpha_tilde(i, 1) + alpha_tilde(i, 2))

! Entropy fix

      a_l = (gamma*p_l/rho_l)**0.5
      lambda_l(1) = u_l - a_l

      u_star = (rho_l*u_l + alpha_tilde(i, 1)*(u_tilde - a_tilde))/(rho_l + alpha_tilde(i, 1))
      rho_star_l = rho_l + alpha_tilde(i, 1)
      p_star = (gamma - 1.0)*(u(i - 1, 3) + alpha_tilde(i, 1)*(H_tilde - u_tilde*a_tilde) - &
          &0.5*rho_star_l*(u_star**2))
      a_star_l = (gamma*p_star/rho_star_l)**0.5
      lambda_r(1) = u_star - a_star_l

      a_r = (gamma*p_r/rho_r)**0.5
      lambda_r(2) = u_r + a_r

      u_star = (rho_r*u_r - alpha_tilde(i, 3)*(u_tilde + a_tilde))/(rho_r - alpha_tilde(i, 3))
      rho_star_r = rho_r - alpha_tilde(i, 3)
      p_star = (gamma - 1.0)*(u(i, 3) - alpha_tilde(i, 3)*(H_tilde + u_tilde*a_tilde) - &
           &0.5*rho_star_r*(u_star**2))
      a_star_r = (gamma*p_star/rho_star_r)**0.5
      lambda_l(2) = u_star + a_star_r

      if ((lambda_l(1) < 0.0) .AND. (lambda_r(1) > 0.0)) then ! check for left transonic rarefraction
         lambda(i, 1) = lambda_l(1)*((lambda_r(1) - lambda(i, 1))/(lambda_r(1) - lambda_l(1)))
      elseif ((lambda_l(2) < 0.0) .AND. (lambda_r(2) > 0.0)) then ! check for right transonic raefraction
         lambda(i, 3) = lambda_r(2)*((lambda(i, 3) - lambda_l(2))/(lambda_r(2) - lambda_l(2)))
      endif

      F(i, 1) = 0.5*(u_tilde*(rho_l + rho_r))
      F(i, 2) = 0.5*(u_tilde*(rho_l*u_l + rho_r*u_r) + p_l + p_r)
      F(i, 3) = 0.5*(u_tilde*(E_l + E_r) + u_l*p_l + u_r*p_r)
   end do

   alpha_tilde(-1, :) = 0
   alpha_tilde(nodes + 1, :) = 0

   do i = 0, nodes
      if (order == 1) then
         limit = 0
      elseif (order == 2) then
         do j = 1, 3 ! loop over waves
            if (lambda(i, j) > 0.0) then
               if (alpha_tilde(i, j) /= 0.0) then
                  theta = alpha_tilde(i - 1, j)/alpha_tilde(i, j)
               else
                  theta = 1.0
               endif
            else
               if (alpha_tilde(i, j) /= 0.0) then
                  theta = alpha_tilde(i + 1, j)/alpha_tilde(i, j)
               else
                  theta = 1.0
               endif
            endif
            limit(j) = limiter(theta, lim_choice)
            theta = 0
         end do
      endif

      do j = 1, 3
         F(i, 1) = F(i, 1) - 0.5*(1.0-limit(j)*(1.0-(delta_t/delta_r)*abs(lambda(i, j))))&
                  &*abs(lambda(i, j))*alpha_tilde(i, j)*K_tilde(i, j, 1)
         F(i, 2) = F(i, 2) - 0.5*(1.0-limit(j)*(1.0-(delta_t/delta_r)*abs(lambda(i, j))))&
                  &*abs(lambda(i, j))*alpha_tilde(i, j)*K_tilde(i, j, 2)
         F(i, 3) = F(i, 3) - 0.5*(1.0-limit(j)*(1.0-(delta_t/delta_r)*abs(lambda(i, j))))&
                  &*abs(lambda(i, j))*alpha_tilde(i, j)*K_tilde(i, j, 3)
      end do

      limit = 0; theta = 0
   end do

   do i = 0, nodes - 1
      u(i, 1) = u(i, 1) - (delta_t/delta_r)*(F(i + 1, 1) - F(i, 1))
      u(i, 2) = u(i, 2) - (delta_t/delta_r)*(F(i + 1, 2) - F(i, 2))
      u(i, 3) = u(i, 3) - (delta_t/delta_r)*(F(i + 1, 3) - F(i, 3))
   end do

   u(-1, 1) = u(0, 1)
   u(-1, 2) = -u(0, 2)
   u(-1, 3) = u(0, 3)

   u(nodes, :) = u(nodes - 1, :)

   return

end subroutine RiemannSolver

subroutine ODESolver(u, delta_r, gamma, delta_t, nodes, alpha)

   implicit none
!--------------------------------------------------------------------------------
   integer, intent(IN) :: nodes, alpha
   double precision, intent(INOUT), dimension(-1:nodes, 1:3) :: u
   double precision, intent(IN) :: gamma, delta_r, delta_t
!--------------------------------------------------------------------------------
   double precision, dimension(1:3) :: uold
   double precision :: p, r
   double precision, external :: Eos
   integer :: i
!--------------------------------------------------------------------------------
   p = 0; r = 0
   do i = 0, nodes - 1

      r = 0.5*delta_r + i*delta_r
      
      p = Eos(u(i, :), gamma)
      
      uold(1) = u(i, 1)
      uold(2) = u(i, 2)
      uold(3) = u(i, 3)

      u(i, 1) = u(i, 1) - dble(alpha - 1)*delta_t*(uold(2)/r)
      u(i, 2) = u(i, 2) - dble(alpha - 1)*delta_t*(((uold(2)**2)/uold(1))/r)
      u(i, 3) = u(i, 3) - dble(alpha - 1)*delta_t*((uold(2)/uold(1))*(uold(3) + p)/r)
   end do

   u(-1, 1) = u(0, 1)
   u(-1, 2) = -u(0, 2)
   u(-1, 3) = u(0, 3)

   u(nodes, :) = u(nodes - 1, :)

   return

end subroutine ODESolver

subroutine AdaptiveTimeStep(delta_t, u, delta_r, nodes, gamma, CFL)

   implicit none
!--------------------------------------------------------------------------------
! GLOBAL VARIABLES
   integer, intent(IN) :: nodes
   double precision, intent(IN), dimension(-1:nodes, 1:3) :: u
   double precision, intent(IN) :: gamma, CFL, delta_r
   double precision, intent(INOUT) :: delta_t
!--------------------------------------------------------------------------------
! LOCAL VARIABLES
   double precision, dimension(-1:nodes) :: P, v
   double precision :: a_l, a_r, max_lambda
   double precision, external :: Eos
   integer :: j
!--------------------------------------------------------------------------------
   max_lambda = 0; P = 0; v = 0

   do j = -1, nodes
      p(j) = Eos(u(j, :), gamma)
      v(j) = u(j, 2)/u(j, 1)
   end do

   do j = 0, nodes - 1
      a_l = SQRT(gamma*(p(j) + p(j - 1))/(u(j, 1) + u(j - 1, 1)))
      a_r = SQRT(gamma*(p(j + 1) + p(j))/(u(j + 1, 1) + u(j, 1)))

      max_lambda = dmax1(abs(0.5*(v(j - 1) + v(j)) + a_l)/delta_r,&
           &abs(0.5*(v(j) + v(j + 1)) + a_r)/delta_r, max_lambda)
      max_lambda = dmax1(abs(0.5*(v(j - 1) + v(j)))/delta_r,&
           &abs(0.5*(v(j) + v(j + 1)))/delta_r, max_lambda)
      max_lambda = dmax1(abs(0.5*(v(j - 1) + v(j)) - a_l)/delta_r,&
           &abs(0.5*(v(j) + v(j + 1)) - a_r)/delta_r, max_lambda)
   end do

   delta_t = CFL/max_lambda

   return

end subroutine AdaptiveTimeStep

double precision function limiter(r, choice)

   implicit none
!---------------------------------------------------------------------------------
! LOCAL VARIABLES
   integer, intent(IN) :: choice
   double precision, intent(IN) :: r
!---------------------------------------------------------------------------------

   if (choice == 1) then
      ! superbee limiter
      limiter = max(0.0, min(2.0*r, 1.0), min(r, 2.0))
   elseif (choice == 2) then
      ! minbee limiter
      limiter = max(0.0, min(r, 1.0))
   elseif (choice == 3) then
      ! Van Leer limiter
      if (r <= 0.0) THEN
         limiter = 0.0
      else
         limiter = 2.0*r/(1.0+r)
      endif
   endif

end function limiter

subroutine FluxFunction(State, Flux, Gamma)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(in)  :: State(3), Gamma
   double precision, intent(out) :: Flux(3)
   double precision, external :: Eos
!-------------------------------------------------------------------------------
   Flux(1) = State(2)
   Flux(2) = State(2)*State(2)/State(1) + Eos(State, Gamma)
   Flux(3) = (State(2)/State(1))*(State(3) + Eos(State, Gamma))

end subroutine FluxFunction

function Eos(State, Gamma) Result(Pressure)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(in) :: State(3), Gamma
   double precision :: Pressure
!-------------------------------------------------------------------------------

! Relationships:
! --------------
! Pressure:         P = (Gamma-1)*(E-0.5*Rho*u^2)
! Enthalpy:         H = (E+P)/Rho = h+0.5*u^2
! Speed:            c = Gamma*P/Rho  = (Gamma-1)*(H-0.5*u^2)
! Internal energy:  e = E/Rho

   Pressure = (Gamma - 1d0)*(State(3) - 0.5d0*State(2)*State(2)/State(1))

   return
end function Eos

subroutine WriteSolution(u, delta_r, gamma, nodes)

   implicit none

!---------------------------------------------------------------------------------
! GLOBALS VARIABLES
   integer, intent(IN) :: nodes
   double precision, intent(IN), dimension(-1:nodes, 1:3) :: u
   double precision, intent(IN) :: delta_r, gamma
   double precision, external :: Eos
!---------------------------------------------------------------------------------
! LOCAL VARIABLES
   character(LEN=10)::G
   parameter(G='(107F20.8)')
   integer :: i
   double precision :: e, p
!---------------------------------------------------------------------------------

   open (unit=20, file='exact.m')
   do i = 0, nodes - 1
      p = Eos(u(i,:),gamma)
      e = P/((gamma - 1.0)*u(i, 1))
      write (20, G) (0.5*delta_r + i*delta_r), u(i, 1), u(i, 2)/u(i, 1), P, e
      p = 0; e = 0
   end do
   close (20)

   return

end subroutine WriteSolution

subroutine ReadVariables(nodes, r_min, r_max, gamma, cfl, output_t, order, &
                          lim_choice, rho_l, u_l, p_l, rho_r, u_r, p_r, r_mid, alpha)
   implicit none

!---------------------------------------------------------------------------------
! GLOBALS VARIABLES
   integer, intent(OUT) :: nodes, order, lim_choice, alpha
   double precision, intent(OUT) :: gamma, r_min, r_max, rho_l, u_l, p_l, rho_r, u_r, p_r
   double precision, intent(OUT) :: r_mid, output_t, CFL
!---------------------------------------------------------------------------------

   open (unit=10, file='variables.data', status='old', form='formatted')

   read (10, *) nodes
   read (10, *) r_min
   read (10, *) r_max
   read (10, *) gamma
   read (10, *) cfl
   read (10, *) output_t
   read (10, *) order
   read (10, *) lim_choice
   read (10, *) rho_l
   read (10, *) u_l
   read (10, *) p_l
   read (10, *) rho_r
   read (10, *) u_r
   read (10, *) p_r
   read (10, *) r_mid
   read (10, *) alpha

   close (10)

end subroutine ReadVariables

subroutine ValidateVariables(nodes, r_min, r_max, r_mid, output_t, cfl, order, lim_choice)
   implicit none

!---------------------------------------------------------------------------------
! GLOBALS VARIABLES
   integer, intent(IN) :: nodes, order, lim_choice
   double precision, intent(IN) :: r_min, r_max, r_mid, output_t, CFL
!---------------------------------------------------------------------------------
! LOCAL VARIABLES
   logical :: run
!---------------------------------------------------------------------------------

   run = .true.
   if (nodes < 2) run = .false.
   if ((r_max - r_min) < 0.0) run = .false.
   if ((r_mid <= r_min) .or. (r_mid >= r_max)) run = .false.
   if (cfl < 0) run = .false.
   if (output_t < 0) run = .false.
   if ((order /= 1) .and. (order /= 2)) run = .false.
   if ((lim_choice /= 1) .and. (lim_choice /= 2) .and. (lim_choice /= 3)) run = .false.

   if (.not. run) stop 'Invalid variables file'

end subroutine ValidateVariables
