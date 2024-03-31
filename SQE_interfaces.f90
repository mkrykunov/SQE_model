module SQE_interfaces
   INTERFACE
      subroutine read_sqe_file(file_name, atoms, cartesian, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii)
         character (len=150), intent(IN) :: file_name
         character (len=5), allocatable, intent(INOUT) :: atoms(:)
         double precision, allocatable, intent(INOUT) :: cartesian(:,:), sqe_radii(:)
         double precision, allocatable, intent(INOUT) :: kappa_a(:), chi_a(:)
         double precision, allocatable, intent(INOUT) :: kappa_b(:), chi_b(:)
         integer, allocatable, intent(INOUT) :: bonds(:,:)
      END subroutine  read_sqe_file
   END INTERFACE

   INTERFACE
      subroutine minimize_energy(xyz, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, alpha, TOL, MAX_ITER, Q)
         integer, intent(IN) :: MAX_ITER
         double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:)
         double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:)
         double precision, allocatable, intent(IN) :: kappa_b(:), chi_b(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, intent(IN) :: alpha, TOL
         double precision, allocatable, intent(INOUT) :: Q(:)
      end subroutine minimize_energy
   END INTERFACE

   INTERFACE
      subroutine minimize_energy_conjgrad(xyz, bonds, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, TOL, MAX_ITER, Q)
         integer, intent(IN) :: MAX_ITER
         double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:)
         double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:)
         double precision, allocatable, intent(IN) :: kappa_b(:), chi_b(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, intent(IN) :: TOL
         double precision, allocatable, intent(INOUT) :: Q(:)
      end subroutine minimize_energy_conjgrad
   END INTERFACE

   INTERFACE
      subroutine sqe_hessian(bonds, kappa_a, kappa_b, Hessian_diag)
         double precision, allocatable, intent(IN) :: kappa_a(:), kappa_b(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, allocatable, intent(INOUT) :: Hessian_diag(:)
      end subroutine sqe_hessian
   END INTERFACE

   INTERFACE
      subroutine sqe_gradient_b(bonds, chi_a, chi_b, gradient_b)
         double precision, allocatable, intent(IN) :: chi_a(:), chi_b(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, allocatable, intent(INOUT) :: gradient_b(:)
      end subroutine sqe_gradient_b
   END INTERFACE

   INTERFACE
      subroutine initial_guess(x, H, b)
         double precision, allocatable, intent(IN) :: H(:), b(:)
         double precision, allocatable, intent(INOUT) :: x(:)
      end subroutine initial_guess
   END INTERFACE

   INTERFACE
      subroutine sqe_gradient(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q, gradient)
         double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:), kappa_b(:), chi_b(:)
         double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:), split_q(:), Q(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, allocatable, intent(INOUT) :: gradient(:)
      end subroutine sqe_gradient
   END INTERFACE

   INTERFACE
      subroutine sqe_gradient_Ax(bonds, xyz, kappa_a, kappa_b, sqe_radii, split_q, Q, gradient)
         double precision, allocatable, intent(IN) :: kappa_a(:), kappa_b(:)
         double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:), split_q(:), Q(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, allocatable, intent(INOUT) :: gradient(:)
      end subroutine sqe_gradient_Ax
   END INTERFACE

   INTERFACE
      subroutine split_q_to_Q(bonds, split_q, Q)
         double precision, allocatable, intent(IN) :: split_q(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
         double precision, allocatable, intent(INOUT) :: Q(:)
      end subroutine split_q_to_Q
   END INTERFACE

   INTERFACE
      function vector_norm(x)
         double precision, allocatable, intent(IN) :: x(:)
      end function vector_norm
   END INTERFACE

   INTERFACE
      function dot_product(x1, x2)
         double precision, allocatable, intent(IN) :: x1(:), x2(:)
      end function dot_product
   END INTERFACE

   INTERFACE
      function sqe_energy(bonds, xyz, kappa_a, chi_a, kappa_b, chi_b, sqe_radii, split_q, Q)
         double precision, allocatable, intent(IN) :: kappa_a(:), chi_a(:), kappa_b(:), chi_b(:)
         double precision, allocatable, intent(IN) :: xyz(:,:), sqe_radii(:), split_q(:), Q(:)
         integer, allocatable, intent(IN) :: bonds(:,:)
      end function sqe_energy
   END INTERFACE

   INTERFACE
      subroutine save_xyz_and_charges(file_name, atoms, xyz, Q)
         character (len=150), intent(IN) :: file_name
         character (len=5), allocatable, intent(IN) :: atoms(:)
         double precision, allocatable, intent(IN) :: xyz(:,:), Q(:)
      end subroutine save_xyz_and_charges
   END INTERFACE

end module