function sqe_energy(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q)
   use SQE_interfaces, except_this_one => sqe_energy

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:), kappa_b(:), chi_b(:)
   double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:), split_q(:), Q(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   integer :: i, j, ij
   double precision :: sqe_energy, R_i, R_j, R_ij, alpha_ij

   sqe_energy = 0.0

   do i = 1, size(Q)
      sqe_energy = sqe_energy + 0.5 * kappa_a(i) * Q(i)**2 + chi_a(i) * Q(i)

      R_i = sqe_radii(i)

      do j = i + 1, size(Q)
         R_j = sqe_radii(j)

         R_ij = sqrt((xyz(i,1) - xyz(j,1))**2 + (xyz(i,2) - xyz(j,2))**2 + (xyz(i,3) - xyz(j,3))**2)
         alpha_ij = sqrt(1.0 / (2 * (R_i**2 + R_j**2)))
         sqe_energy = sqe_energy + Q(i) * Q(j) * erf(alpha_ij * R_ij) / R_ij
      end do
   end do

   do ij = 1, size(split_q)
      i = bonds(ij,1)
      j = bonds(ij,2)

      sqe_energy = sqe_energy + 0.5 * kappa_b(ij) * split_q(ij)**2 + chi_b(ij) * split_q(ij)
   end do

end function sqe_energy
