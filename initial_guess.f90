subroutine initial_guess(x, H, b)
   use SQE_interfaces, except_this_one => initial_guess

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: H(:), b(:)
   double precision, allocatable, intent(INOUT) :: x(:)
   integer :: ij

   do ij = 1, size(x) ! num_bonds
      if (abs(H(ij)) > 1.0e-5) then
         x(ij) = b(ij) / H(ij)
      else
          x(ij) = b(ij)
      end if
   end do

end subroutine initial_guess
