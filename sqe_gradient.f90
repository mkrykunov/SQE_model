subroutine sqe_gradient(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q, gradient)
   use SQE_interfaces, except_this_one => sqe_gradient

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:), kappa_b(:), chi_b(:)
   double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:), split_q(:), Q(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   double precision, allocatable, intent(INOUT) :: gradient(:)
   integer :: i, j, m, ij
   double precision :: R_i, R_j, R_m, R_im, R_jm, alpha_im, alpha_jm

   do ij = 1, size(split_q)
      i = bonds(ij,1)
      j = bonds(ij,2)

      gradient(ij) = kappa_b(ij) * split_q(ij) + chi_b(ij)
      gradient(ij) = gradient(ij) + kappa_a(i) * Q(i) - kappa_a(j) * Q(j)
      gradient(ij) = gradient(ij) + (chi_a(i) - chi_a(j))

      R_i = sqe_radii(i)
      R_j = sqe_radii(j)

      do m = 1, size(Q)
         R_m = sqe_radii(m)

         if (i.NE.m) then
            R_im = sqrt((xyz(i,1) - xyz(m,1))**2 + (xyz(i,2) - xyz(m,2))**2 + (xyz(i,3) - xyz(m,3))**2)
            alpha_im = sqrt(1.0 / (2 * R_i**2 + 2 * R_m**2))
            gradient(ij) = gradient(ij) + Q(m) * erf(alpha_im * R_im) / R_im
         end if

         if (j.NE.m) then
            R_jm = sqrt((xyz(j,1) - xyz(m,1))**2 + (xyz(j,2) - xyz(m,2))**2 + (xyz(j,3) - xyz(m,3))**2)
            alpha_jm = sqrt(1.0 / (2 * R_j**2 + 2 * R_m**2))
            gradient(ij) = gradient(ij) - Q(m) * erf(alpha_jm * R_jm) / R_jm
         end if

      end do
   end do

end subroutine sqe_gradient
