! COMPILATION:
!
! gfortran -c SQE_interfaces.f90
!
! gfortran read_sqe_file.f90 sqe.f90 minimize_energy.f90 minimize_energy_conjgrad.f90 sqe_hessian.f90 sqe_gradient_b.f90
! initial_guess.f90 sqe_gradient.f90 sqe_gradient_Ax.f90 split_q_to_Q.f90 vector_norm.f90 dot_product.f90 sqe_energy.f90
! save_xyz_and_charges.f90 -o calc_sqe.exe
!
! EXECUTION (EXAMPLE):
!
! calc_sqe.exe Dat/Naphthalene_SQE.inp Dat/output_charges.xyz
!
program sqe
   use SQE_interfaces

   IMPLICIT NONE

   character (len=5), allocatable :: atoms(:)
   double precision, allocatable :: coords(:,:), sqe_radii(:)
   double precision, allocatable :: kappa_a(:), chi_a(:)
   double precision, allocatable :: kappa_b(:), chi_b(:)
   integer, allocatable :: bonds(:,:)
   character (len=150) :: in_file_name, out_file_name
   integer :: num_atoms, conjgrad_opt
   integer, parameter :: MAX_ITER = 1000
   double precision, parameter :: alpha = 0.02, TOL = 1.e-8
   double precision, allocatable :: Q(:)
   integer :: i, nargs
   character (len=80) :: arg

   in_file_name = "Dat/Naphthalene_SQE.inp"
   out_file_name = "Dat/output_charges.xyz"

   nargs = iargc() 
   do i = 0, nargs 
      call getarg(i, arg) 
      print '(a)', arg
      if (i == 1) in_file_name = arg
      if (i == 2) out_file_name = arg
   end do

   call read_sqe_file(in_file_name, atoms, coords, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii)

   if ((size(atoms) .NE. size(kappa_a)) .OR. (size(atoms) .NE. size(coords, 1))) then
      write(*,*) "Something is wrong with dimensions:", size(atoms), size(kappa_a), size(coords, 1)
      return
   end if

   num_atoms = size(atoms)

   allocate(Q(num_atoms))

   conjgrad_opt = 1

   if (conjgrad_opt == 1) then
      call minimize_energy_conjgrad(coords, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, TOL, MAX_ITER, Q)
   else
      call minimize_energy(coords, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, alpha, TOL, MAX_ITER, Q)
   end if

   call save_xyz_and_charges(out_file_name, atoms, coords, Q)

   deallocate(Q)

   deallocate(kappa_a)
   deallocate(chi_a)
   deallocate(kappa_b)
   deallocate(chi_b)
   deallocate(sqe_radii)

   deallocate(atoms)
   deallocate(bonds)
   deallocate(coords)
end program sqe