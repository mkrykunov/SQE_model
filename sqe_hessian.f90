subroutine sqe_hessian(bonds, kappa_a, kappa_b, Hessian_diag)
   use SQE_interfaces, except_this_one => sqe_hessian

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: kappa_a(:), kappa_b(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   double precision, allocatable, intent(INOUT) :: Hessian_diag(:)
   integer :: i, j, ij

   do ij = 1, size(Hessian_diag)
      i = bonds(ij,1)
      j = bonds(ij,2)

      Hessian_diag(ij) = kappa_b(ij)
      Hessian_diag(ij) = Hessian_diag(ij) + kappa_a(i) + kappa_a(j)
   end do

end subroutine sqe_hessian
