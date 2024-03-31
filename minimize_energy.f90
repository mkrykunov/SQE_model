subroutine minimize_energy(xyz, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, alpha, TOL, MAX_ITER, Q)
   use SQE_interfaces, except_this_one => minimize_energy

   IMPLICIT NONE

   integer, intent(IN) :: MAX_ITER
   double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:)
   double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:)
   double precision, allocatable, intent(IN) :: kappa_b(:), chi_b(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   double precision, intent(IN) :: alpha, TOL
   double precision, allocatable, intent(INOUT) :: Q(:)
   double precision, allocatable :: Hessian_diag(:), gradient(:), split_q(:)
   integer :: ij, istep, num_bonds
   double precision :: b_norm, q_norm, grad_norm, energy

   if ((size(kappa_b) .NE. size(chi_b)) .OR. (size(kappa_b) .NE. size(bonds, 1))) then
      write(*,*) "Something is wrong with dimensions:", size(kappa_b), size(chi_b), size(bonds, 1)
      return
   end if

   num_bonds = size(kappa_b)

   allocate(Hessian_diag(num_bonds))
   allocate(gradient(num_bonds))
   allocate(split_q(num_bonds))

   call sqe_hessian(bonds, kappa_a, kappa_b, Hessian_diag)

   call sqe_gradient_b(bonds, chi_a, chi_b, gradient)

   b_norm = vector_norm(gradient)
   write(*,"(A8, 2X, F10.5)") "b_norm =", b_norm

   call initial_guess(split_q, Hessian_diag, gradient);

   do istep = 1, MAX_ITER
      call split_q_to_Q(bonds, split_q, Q)

      call sqe_gradient(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q, gradient)

      split_q = split_q - alpha * gradient

      grad_norm = vector_norm(gradient)

      call split_q_to_Q(bonds, split_q, Q)
      energy = sqe_energy(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q)

      if (mod(istep - 1,50) == 0) then
         write(*,"(A5, 2X, I5, 2X, E10.5, 2X, F10.5)") "step:", istep, grad_norm / b_norm, energy
      end if

      if (grad_norm / b_norm < TOL) exit
   end do

   deallocate(split_q)
   deallocate(gradient)
   deallocate(Hessian_diag)

end subroutine minimize_energy


