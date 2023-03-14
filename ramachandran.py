import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PPBuilder, calc_dihedral

# Load the PDB file containing EGFR-EGF data
parser = PDBParser()
structure = parser.get_structure("EGFR-EGF", "egfr-egf.pdb")

# Extract the protein chains from the structure
chains = list(structure.get_chains())

# Initialize a list to store phi/psi angles for each residue
angles = []

# Iterate over each residue in each chain and calculate phi/psi angles
for chain in chains:
    residues = chain.get_residues()
    for residue in residues:
        if residue.get_resname() == "GLY":
            continue  # skip glycine residues
        try:
            phi, psi = calc_dihedral(
                residue["C"].get_vector(),
                residue["N"].get_vector(),
                residue["CA"].get_vector(),
                residue["C"].get_vector(),
            )
            angles.append((phi, psi))
        except:
            continue  # skip residues with missing atoms or bad geometry

# Plot the Ramachandran plot using matplotlib

fig, ax = plt.subplots()
ax.set_xlabel("Phi")
ax.set_ylabel("Psi")
ax.set_title("Ramachandran Plot")

# Plot the data points
x, y = zip(*angles)
ax.scatter(x, y)

# Identify outliers
outliers = []
for angle in angles:
    if angle[0] < -180 or angle[0] > 180 or angle[1] < -180 or angle[1] > 180:
        outliers.append(angle)

# Plot the outliers as red points
if len(outliers) > 0:
    ox, oy = zip(*outliers)
    ax.scatter(ox, oy, color="red")

plt.show()
