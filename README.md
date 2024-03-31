# SQE_model
SQE charge equilibration model
There are two parts of the code. The first one is in Python. It reads mol2 and csv files and prepares the input file for the Fortran code.
Example:
python mol2_to_sqe.py Dat/Naphthalene_coords.mol2 Dat/sqe_params.csv Dat/Naphthalene_SQE.inp
Then the Fortran code reads inp file and performs the calculations of partial atomic charges based on the SQE split-charge equilibration model.
Example:
calc_sqe.exe Dat/Naphthalene_SQE.inp Dat/output_charges.xyz
The results (partial atomic charges) are in the extended XYZ format in the fifth column.
The parameters of the model are in Dat/sqe_params.csv, and they are not optimized yet.
The code is for the neutral molecules yet and does not take into account the external electric field.
