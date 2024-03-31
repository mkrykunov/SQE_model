import sys
import pandas as pd

# 1.5 times the size of covalent radii in Bohr in accordance to other sqe implementations
#
sqe_radii = {
"H": 1.75744516977 ,
"He": 0.79368491538 ,
"Li": 3.62827389888 ,
"Be": 2.72120542416 ,
"B": 2.38105474614 ,
"C": 2.15428762746 ,
"N": 2.01255817828 ,
"O": 1.87082872911 ,
"F": 1.61571572059 ,
"Ne": 1.64406161043 ,
"Na": 4.70541771261 ,
"Mg": 3.99677046673 ,
"Al": 3.42985267003 ,
"Si": 3.14639377169 ,
"P": 3.03301021234 ,
"S": 2.97631843268 ,
"Cl": 2.89128076317 ,
"Ar": 3.00466432251 ,
"K": 5.75421563651 ,
"Ca": 4.98887661096 ,
"Sc": 4.81880127195 ,
"Ti": 4.5353423736 ,
"V": 4.33692114476 ,
"Cr": 3.94007868706 ,
"Mn": 3.94007868706 ,
"Fe": 3.74165745822 ,
"Co": 3.57158211921 ,
"Ni": 3.51489033954 ,
"Cu": 3.74165745822 ,
"Zn": 3.45819855987 ,
"Ga": 3.45819855987 ,
"Ge": 3.4015067802 ,
"As": 3.37316089036 ,
"Se": 3.4015067802 ,
"Br": 3.4015067802 ,
"Kr": 3.28812322086 ,
"Rb": 6.2360957637 ,
"Sr": 5.52744851782 ,
"Y": 5.38571906865 ,
"Zr": 4.96053072113 ,
"Nb": 4.64872593294 ,
"Mo": 4.36526703459 ,
"Tc": 4.16684580574 ,
"Ru": 4.13849991591 ,
"Rh": 4.02511635657 ,
"Pd": 3.94007868706 ,
"Ag": 4.11015402607 ,
"Cd": 4.08180813624 ,
"In": 4.02511635657 ,
"Sn": 3.94007868706 ,
"Sb": 3.94007868706 ,
"Te": 3.91173279723 ,
"I": 3.94007868706 ,
"Xe": 3.9684245769 ,
"Cs": 6.91639711974 ,
"Ba": 6.09436631452 ,
"La": 5.86759919584 ,
"Ce": 5.78256152634 ,
"Pr": 5.75421563651 ,
"Nd": 5.69752385683 ,
"Pm": 5.64083207716 ,
"Sm": 5.61248618733 ,
"Eu": 5.61248618733 ,
"Gd": 5.55579440766 ,
"Tb": 5.49910262799 ,
"Dy": 5.44241084832 ,
"Ho": 5.44241084832 ,
"Er": 5.35737317882 ,
"Tm": 5.38571906865 ,
"Yb": 5.30068139915 ,
"Lu": 5.30068139915 ,
"Hf": 4.96053072113 ,
"Ta": 4.81880127195 ,
"W": 4.59203415327 ,
"Re": 4.28022936508 ,
"Os": 4.08180813624 ,
"Ir": 3.99677046673 ,
"Pt": 3.85504101756 ,
"Au": 3.85504101756 ,
"Hg": 3.74165745822 ,
"Tl": 4.11015402607 ,
"Pb": 4.13849991591 ,
"Bi": 4.19519169558 ,
"Po": 3.9684245769 ,
"At": 4.25188347525 ,
"Rn": 4.25188347525 ,
"Fr": 7.3699313571 ,
"Ra": 6.26444165354 ,
"Ac": 6.09436631452 ,
"Th": 5.83925330601 ,
"Pa": 5.669177967 ,
"U": 5.55579440766 ,
"Np": 5.38571906865 ,
"Pu": 5.30068139915 ,
"Am": 5.1022601703 ,
"Cm": 4.79045538212 ,
"Bk": 4.5353423736 ,
"Cf": 4.5353423736 ,
"Es": 4.5353423736 ,
"Fm": 4.5353423736 ,
"Md": 4.5353423736 ,
"No": 4.5353423736 ,
"Lr": 4.5353423736 ,
"Rf": 4.5353423736 ,
"Db": 4.5353423736 ,
"Sg": 4.5353423736 ,
"Bh": 4.5353423736 ,
"Hs": 4.5353423736 ,
"Mt": 4.5353423736 ,
"Ds": 4.5353423736 ,
"Rg": 4.5353423736 ,
"Cn": 4.5353423736 ,
"Uut": 4.5353423736 ,
"Fl": 4.5353423736 ,
"Uup": 4.5353423736 ,
"Lv": 4.5353423736 ,
"Uuh": 4.5353423736 
}

def read_mol2(mol2_file):
    finp = open(mol2_file, 'r')
    lines = finp.readlines()
    finp.close()

    num_atoms = 0
    num_bonds = 0

    atoms = list()
    xyz = list()
    bonds = list()

    for k, line in enumerate(lines):
        if "@<TRIPOS>MOLECULE" in line:
            num_atoms = int(lines[k + 2].split()[0])
            num_bonds = int(lines[k + 2].split()[1])
        
        if "@<TRIPOS>ATOM" in line and num_atoms != 0:
            for idx in range(k + 1, k + num_atoms + 1):
                atom_str = lines[idx].split()[1]
                xyz_list = [float(item) for item in lines[idx].split()[2:5]]
                atoms.append(atom_str)
                xyz.append(xyz_list)
            
        if "@<TRIPOS>BOND" in line and num_bonds != 0:
            for idx in range(k + 1, k + num_bonds + 1):
                bonds_list = [int(item) for item in lines[idx].split()[1:3]]
                bonds.append(bonds_list)

    return atoms, xyz, bonds


def read_params(csv_file):
    df = pd.read_csv(csv_file)

    kappa_a = dict()
    chi_a = dict()
    kappa_b = dict()
    chi_b = dict()

    for k, item in enumerate(df['atoms_bonds']):
        if "_" in item:
            bond_ab = item
            bond_ba = item.split('_')[-1] + "_" + item.split('_')[0]
            kappa_b[bond_ab] = df['kappa'][k]
            kappa_b[bond_ba] = df['kappa'][k]
            chi_b[bond_ab] = df['chi'][k]
            if bond_ab == bond_ba:
                chi_b[bond_ab] = 0.0
            else:
                chi_b[bond_ba] = -df['chi'][k]
        else:
            kappa_a[item] = df['kappa'][k]
            chi_a[item] = df['chi'][k]

    return kappa_a, chi_a, kappa_b, chi_b


def write_sqe_input(output_file, atoms, xyz, bonds, kappa_a, chi_a, kappa_b, chi_b):
    fout = open(output_file, "w")
    fout.write("%4d %4d\n" % (len(atoms), len(bonds)))
    for i in range(len(atoms)):
        fout.write("%4d %3s %9.4f %9.4f %9.4f  %11.6f %11.6f %11.6f \n" % (i + 1, atoms[i], xyz[i][0], xyz[i][1], xyz[i][2], 
                                                                           kappa_a[atoms[i]], chi_a[atoms[i]], sqe_radii[atoms[i]]))
    for ij in range(len(bonds)):
        i = bonds[ij][0]
        j = bonds[ij][1]
        fout.write("%4d %4d %4d  %11.6f %11.6f\n" % (ij + 1, i, j, kappa_b[atoms[i - 1] + "_" + atoms[j - 1]],
                                                    chi_b[atoms[i - 1] + "_" + atoms[j - 1]]))
    fout.close()

def mol2_to_sqe(mol2_file, csv_file, sqe_input_file):
    atoms, xyz, bonds = read_mol2(mol2_file)
    kappa_a, chi_a, kappa_b, chi_b = read_params(csv_file)
    write_sqe_input(sqe_input_file, atoms, xyz, bonds, kappa_a, chi_a, kappa_b, chi_b)


# example: python mol2_to_sqe.py Dat/Naphthalene_coords.mol2 Dat/sqe_params.csv Dat/Naphthalene_SQE.inp
#
if __name__ == '__main__':
    print (sys.argv[1], sys.argv[2], sys.argv[3])
    
    mol2_file = sys.argv[1]
    csv_file = sys.argv[2]
    sqe_input_file = sys.argv[3]

    mol2_to_sqe(mol2_file, csv_file, sqe_input_file)