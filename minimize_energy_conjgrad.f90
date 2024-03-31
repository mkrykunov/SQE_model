subroutine minimize_energy_conjgrad(xyz, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, TOL, MAX_ITER, Q)
   use SQE_interfaces, except_this_one => minimize_energy_conjgrad

   IMPLICIT NONE

   integer, intent(IN) :: MAX_ITER
   double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:)
   double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:)
   double precision, allocatable, intent(IN) :: kappa_b(:), chi_b(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   double precision, intent(IN) :: TOL
   double precision, allocatable, intent(INOUT) :: Q(:)
   double precision, allocatable :: Hessian_diag(:), gradient(:), split_q(:), Ap(:), r(:), p(:)
   integer :: ij, istep, num_bonds
   double precision :: b_norm, q_norm, grad_norm, energy, rsold, rsnew, pAp, alpha

   if ((size(kappa_b) .NE. size(chi_b)) .OR. (size(kappa_b) .NE. size(bonds, 1))) then
      write(*,*) "Something is wrong with dimensions:", size(kappa_b), size(chi_b), size(bonds, 1)
      return
   end if

   num_bonds = size(kappa_b)

   allocate(Hessian_diag(num_bonds))
   allocate(gradient(num_bonds))
   allocate(split_q(num_bonds))
   allocate(Ap(num_bonds))
   allocate(r(num_bonds))
   allocate(p(num_bonds))

   call sqe_hessian(bonds, kappa_a, kappa_b, Hessian_diag)

   call sqe_gradient_b(bonds, chi_a, chi_b, gradient)

   b_norm = vector_norm(gradient)
   write(*,"(A8, 2X, F10.5)") "b_norm =", b_norm

   call initial_guess(split_q, Hessian_diag, gradient)

   call split_q_to_Q(bonds, split_q, Q)
   call sqe_gradient_Ax(bonds, xyz, kappa_a, kappa_b, sqe_radii, split_q, Q, Ap)

   ! r = b - A * x
   !
   r = gradient - Ap

   ! p = r
   !
   p = r

   ! rsold = r' * r
   !
   rsold = dot_product(r, r)

   do istep = 1, MAX_ITER
      call split_q_to_Q(bonds, p, Q)

      ! Ap = A * p
      !
      call sqe_gradient_Ax(bonds, xyz, kappa_a, kappa_b, sqe_radii, p, Q, Ap)

      ! Calculating p' * Ap;
      !
      pAp = dot_product(p, Ap)

      ! alpha = rsold / (p' * Ap)
      !
      alpha = rsold / pAp

      ! x = x + alpha * p
      !
      split_q = split_q + alpha * p

      ! r = r - alpha * Ap
      !
      r = r - alpha * Ap

      ! rsnew = r' * r
      !
      rsnew = dot_product(r, r)

      ! p = r + (rsnew / rsold) * p
      !
      p = r + (rsnew / rsold) * p

      ! rsold = rsnew
      !
      rsold = rsnew

      if (mod(istep - 1, 5) == 0) then
         call split_q_to_Q(bonds, split_q, Q)
         energy = sqe_energy(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q)

         write(*,"(A5, 2X, I5, 2X, E10.5, 2X, F10.5)") "step:", istep, rsnew / b_norm, energy
      end if

      if (rsnew / b_norm < TOL) exit
   end do

   call split_q_to_Q(bonds, split_q, Q)
   energy = sqe_energy(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q)

   write(*,"(A8, 2X, F10.5)") "energy =", energy

   deallocate(p)
   deallocate(r)
   deallocate(Ap)
   deallocate(split_q)
   deallocate(gradient)
   deallocate(Hessian_diag)

end subroutine minimize_energy_conjgrad


