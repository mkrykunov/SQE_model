calc_sqe.exe: sqe.o read_sqe_file.o minimize_energy.o minimize_energy_conjgrad.o sqe_hessian.o sqe_gradient_b.o initial_guess.o sqe_gradient.o sqe_gradient_Ax.o split_q_to_Q.o vector_norm.o dot_product.o sqe_energy.o save_xyz_and_charges.o SQE_interfaces.mod
	gfortran read_sqe_file.o sqe.o minimize_energy.o minimize_energy_conjgrad.o sqe_hessian.o sqe_gradient_b.o initial_guess.o sqe_gradient.o sqe_gradient_Ax.o split_q_to_Q.o vector_norm.o dot_product.o sqe_energy.o save_xyz_and_charges.o -o calc_sqe.exe

sqe.o: sqe.f90 SQE_interfaces.f90 read_sqe_file.f90 minimize_energy.f90 minimize_energy_conjgrad.f90 sqe_hessian.f90 sqe_gradient_b.f90 initial_guess.f90 sqe_gradient.f90 sqe_gradient_Ax.f90 split_q_to_Q.f90 vector_norm.f90 dot_product.f90 sqe_energy.f90 save_xyz_and_charges.f90
	gfortran -c SQE_interfaces.f90
	gfortran -c sqe.f90
	gfortran -c read_sqe_file.f90
	gfortran -c minimize_energy.f90
	gfortran -c minimize_energy_conjgrad.f90
	gfortran -c sqe_hessian.f90
	gfortran -c sqe_gradient_b.f90
	gfortran -c initial_guess.f90
	gfortran -c sqe_gradient.f90
	gfortran -c sqe_gradient_Ax.f90
	gfortran -c split_q_to_Q.f90
	gfortran -c vector_norm.f90
	gfortran -c dot_product.f90
	gfortran -c sqe_energy.f90
	gfortran -c save_xyz_and_charges.f90

clean:
	rm *.o *.mod calc_sqe.exe
