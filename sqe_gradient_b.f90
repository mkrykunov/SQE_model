subroutine sqe_gradient_b(bonds, chi_a, chi_b, gradient_b)
   use SQE_interfaces, except_this_one => sqe_gradient_b

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: chi_a(:), chi_b(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   double precision, allocatable, intent(INOUT) :: gradient_b(:)
   integer :: i, j, ij

   do ij = 1, size(gradient_b)
      i = bonds(ij,1)
      j = bonds(ij,2)

      gradient_b(ij) = -chi_b(ij)
      gradient_b(ij) = gradient_b(ij) - (chi_a(i) - chi_a(j))
   end do

end subroutine sqe_gradient_b
