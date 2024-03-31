function vector_norm(x)
   use SQE_interfaces, except_this_one => vector_norm

   IMPLICIT NONE

   double precision :: vector_norm
   double precision, allocatable, intent(IN) :: x(:)
   integer :: i
   double precision :: result

   result = 0.0
   do i = 1, size(x)
      result = result + x(i)**2
   end do

   vector_norm = sqrt(result)

end function vector_norm
