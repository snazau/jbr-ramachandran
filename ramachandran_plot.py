from Bio.PDB import *
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set()


def get_atom(residue, atom_key, is_success):
	try:
		atom = residue[atom_key].get_coord()
		is_success = is_success and True
	except Exception:
		atom = None
		is_success = is_success and False
	return atom, is_success


def calc_dihedral_angle(p1, p2, p3, p4):
	q1 = p2 - p1
	q2 = p3 - p2
	q3 = p4 - p3

	m1 = np.cross(q1, q2)
	m2 = np.cross(q2, q3)

	n1 = m1 / np.sqrt(sum(m1 * m1))
	n2 = m2 / np.sqrt(sum(m2 * m2))

	u1 = n2
	u3 = q2 / np.sqrt(sum(q2 * q2))
	u2 = np.cross(u3, u1)

	cos_theta = np.dot(n1, u1)
	sin_theta = np.dot(n1, u2)

	theta = -math.atan2(sin_theta,cos_theta)
	return np.degrees(theta)


def save_ramachandran_plot(structure, save_dir, save_dpi=150):
	for model in structure:
		for chain in model:
			name, digit, letter = chain.get_full_id()
			print(name, digit, letter)

			phis = []
			psis = []
			for residue_index, next_residue in enumerate(chain):
				if residue_index == 0:
					prev_residue = next_residue
					continue
				if residue_index == 1:
					curr_residue = next_residue
					continue

				is_success = True
				n, is_success = get_atom(curr_residue, "N", is_success)
				ca, is_success = get_atom(curr_residue, "CA", is_success)
				c, is_success = get_atom(curr_residue, "C", is_success)
				cp, is_success = get_atom(prev_residue, "C", is_success)
				nn, is_success = get_atom(next_residue, "N", is_success)

				if is_success:
					phi = calc_dihedral_angle(cp, n, ca, c)
					psi = calc_dihedral_angle(n, ca, c, nn)
					phis.append(phi)
					psis.append(psi)

				prev_residue = curr_residue
				curr_residue = next_residue

			g = sns.jointplot(x=phis, y=psis, linewidth=0.5)
			g.set_axis_labels('phi', 'psi')
			g.ax_joint.grid(True)
			save_path = os.path.join(save_dir, name + "-" + str(digit) + "-" + letter + ".png")
			plt.savefig(save_path, dpi=save_dpi)


if __name__ == "__main__":
	pdb_parser = PDBParser()
	pdb_dir = "./data"

	ramachandran_plots_dir = "./plots"
	if not os.path.exists(ramachandran_plots_dir):
		os.makedirs(ramachandran_plots_dir)

	for filename in os.listdir(pdb_dir):
		if filename.endswith(".pdb"):
			name = os.path.splitext(filename)[0]

			pdb_path = os.path.join(pdb_dir, filename)
			structure = pdb_parser.get_structure(name, pdb_path)

			save_dir = os.path.join(ramachandran_plots_dir, name)
			if not os.path.exists(save_dir):
				os.makedirs(save_dir)

			save_ramachandran_plot(structure, save_dir)
