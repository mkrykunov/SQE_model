subroutine split_q_to_Q(bonds, split_q, Q)
   use SQE_interfaces, except_this_one => split_q_to_Q

   IMPLICIT NONE

   double precision, allocatable, intent(IN) :: split_q(:)
   integer, allocatable, intent(IN) :: bonds(:,:)
   double precision, allocatable, intent(INOUT) :: Q(:)
   integer :: i, ij, iatom_1, iatom_2

   do i = 1, size(Q)
      Q(i) = 0.0

      do ij = 1, size(split_q)
         iatom_1 = bonds(ij,1)
         iatom_2 = bonds(ij,2)

         if (i == iatom_1) then
            Q(i) = Q(i) + split_q(ij)
         else if (i == iatom_2) then
            Q(i) = Q(i) - split_q(ij)
         end if
      end do
   end do

end subroutine split_q_to_Q
