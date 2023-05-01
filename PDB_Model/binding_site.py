from Bio.PDB import PDBParser, Selection

# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure("EGFR", "PDB_Model/1nql.pdb")

# Select the chain and residue range for the kinase domain
chain = structure[0]['A']
start_res = 670
end_res = 998

# Define the selection function to extract the residues in the kinase domain


def select_kinase(residue):
    if residue.get_id()[1] >= start_res and residue.get_id()[1] <= end_res:
        return True
    else:
        return False


# Use the Selection module to extract the kinase residues
kinase_residues = Selection.unfold_entities(chain, 'R')
kinase_residues = list(filter(select_kinase, kinase_residues))

# Print the list of kinase residues
print("Kinase Residues:")
for residue in kinase_residues:
    print(residue.get_resname(), residue.get_id()[1])
