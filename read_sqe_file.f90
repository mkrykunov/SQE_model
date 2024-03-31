subroutine read_sqe_file(file_name, atoms, xyz, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii)
   use SQE_interfaces, except_this_one => read_sqe_file

   IMPLICIT NONE

   character (len=150), intent(IN) :: file_name
   character (len=5), allocatable, intent(INOUT) :: atoms(:)
   double precision, allocatable, intent(INOUT) :: xyz(:,:), sqe_radii(:)
   double precision, allocatable, intent(INOUT) :: kappa_a(:), chi_a(:)
   double precision, allocatable, intent(INOUT) :: kappa_b(:), chi_b(:)
   integer, allocatable, intent(INOUT) :: bonds(:,:)
   integer :: num_atoms, num_bonds

   character (len=50) :: Format_1, Format_2
   integer :: i, idx
   character (len=5) atom_str

   Format_1 = "(I6, 2X, A5, 3F10.5, 2X, 3F10.6)"
   Format_2 = "(I6, 2X, 2I6, 2X, 2F10.6)"

   open(10, file=file_name, status="old", action="read")
      read(10, *) num_atoms, num_bonds
      write(*, "(2I6)") num_atoms, num_bonds

      allocate(xyz(num_atoms,3))
      allocate(bonds(num_bonds,2))
      allocate(atoms(num_atoms))

      allocate(sqe_radii(num_atoms))

      allocate(kappa_a(num_atoms))
      allocate(chi_a(num_atoms))
      allocate(kappa_b(num_bonds))
      allocate(chi_b(num_bonds))

      do i = 1, num_atoms
         read(10, *) idx, atom_str, xyz(idx,1:3), kappa_a(idx), chi_a(idx), sqe_radii(idx)
         if (num_atoms < 20) write(*, Format_1) idx, atom_str, xyz(idx,1:3), kappa_a(idx), chi_a(idx), sqe_radii(idx)
         atoms(idx) = atom_str
         xyz(idx,1) = xyz(idx,1) * 1.8897259886 ! Angstrom to Bohr conversion
         xyz(idx,2) = xyz(idx,2) * 1.8897259886 ! Angstrom to Bohr conversion
         xyz(idx,3) = xyz(idx,3) * 1.8897259886 ! Angstrom to Bohr conversion
      end do

      do i = 1, num_bonds
         read(10, *) idx, bonds(idx,1:2), kappa_b(idx), chi_b(idx)
         if (num_bonds < 20) write(*, Format_2) idx, bonds(idx,1:2), kappa_b(idx), chi_b(idx)
      end do
   close(10)

end subroutine read_sqe_file
