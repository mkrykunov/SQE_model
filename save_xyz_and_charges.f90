subroutine save_xyz_and_charges(file_name, atoms, xyz, Q)
   use SQE_interfaces, except_this_one => save_xyz_and_charges

   IMPLICIT NONE

   character (len=150), intent(IN) :: file_name
   character (len=5), allocatable, intent(IN) :: atoms(:)
   double precision, allocatable, intent(IN) :: xyz(:,:), Q(:)
   double precision :: x, y, z
   integer :: i, num_atoms
   logical :: exist
   character (len=30), parameter :: info = "Format: atom X Y Z charge"

   if ((size(atoms) .NE. size(Q)) .OR. (size(atoms) .NE. size(xyz, 1))) then
      write(*,*) "Something is wrong with dimensions:", size(atoms), size(Q), size(xyz, 1)
      return
   end if

   num_atoms = size(atoms)

   inquire(file=file_name, exist=exist)

   if (exist) then
      open(10, file=file_name, status="old", action="write")
   else
      open(10, file=file_name, status="new", action="write")
   end if
      write(10, "(I0)") num_atoms
      write(10, "(A26)") info

      do i = 1, num_atoms
         x = xyz(i,1) / 1.8897259886 ! to Angstrom
         y = xyz(i,2) / 1.8897259886 ! to Angstrom
         z = xyz(i,3) / 1.8897259886 ! to Angstrom

         write(10, "(A5, 2X, 3F10.5, 2X, F10.6)") atoms(i), x, y, z, Q(i)
      end do
   close(10)

end subroutine save_xyz_and_charges
