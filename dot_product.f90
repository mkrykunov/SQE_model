function dot_product(x1, x2)
   use SQE_interfaces, except_this_one => dot_product

   IMPLICIT NONE

   double precision :: dot_product
   double precision, allocatable, intent(IN) :: x1(:), x2(:)
   integer :: i

   dot_product = 0.0

   if (size(x1) == size(x2)) then
      do i = 1, size(x1)
         dot_product = dot_product + x1(i) * x2(i)
      end do
   end if

end function dot_product
