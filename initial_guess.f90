subroutine initial_guess(x, H, b)
   use SQE_interfaces, except_this_one => initial_guess

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: H(:), b(:)
   double precision, allocatable, intent(INOUT) :: x(:)

   where (abs(H) > 1.0e-5)
      x = b / H
   elsewhere
      x = b
   endwhere

end subroutine initial_guess
