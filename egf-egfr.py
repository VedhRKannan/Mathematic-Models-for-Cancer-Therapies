from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import nglview
from Bio import SeqIO
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser


# Download the PDB files for 1EGF and 4HJO
pdb_list = PDBList()
pdb_list.retrieve_pdb_file('1EGF')
pdb_list.retrieve_pdb_file('4HJO')

# Load the PDB files for EGF and EGFR
parser = PDBParser()
egf_structure = parser.get_structure('EGF', 'eg/1egf.cif')
egfr_structure = parser.get_structure('EGFR', 'hj/4hjo.cif')

# Load the UniProt sequence for EGFR
handle = open('P00533.fasta')
seq_record = SeqIO.read(handle, 'fasta')
handle.close()

# Extract the sequence of the kinase domain of EGFR (residues 711-980)
kinase_domain = seq_record.seq[710:979]

# Find the chain that contains the kinase domain of EGFR in the PDB structure
kinase_domain_chain = None
for chain in egfr_structure.get_chains():
    for residue in chain.get_residues():
        if residue.get_resname() == 'TPO':
            if chain.id in residue.get_full_id():
                kinase_domain_chain = chain
                break
    if kinase_domain_chain is not None:
        break

# Print the sequence of the kinase domain in the PDB file
if kinase_domain_chain is not None:
    pdb_seq = ''
    for residue in kinase_domain_chain.get_residues():
        if residue.get_resname() == 'TPO':
            pdb_seq += 'T'
        else:
            pdb_seq += residue.get_resname()
    print(pdb_seq)
else:
    print('Kinase domain chain not found in PDB file')


# pdb_list = PDBList()
# pdb_list.retrieve_pdb_file('1EGF')
# pdb_list.retrieve_pdb_file('4HJO')


# parser = MMCIFParser()

# # Load EGF structure
# egf_structure = parser.get_structure('EGF', 'pdb1egf.ent')
# egf_chain = egf_structure[0]['A']

# # Load EGFR structure and extract active site chain
# egfr_structure = parser.get_structure('EGFR', 'pdb4hjo.ent')
# active_site_chain = None
# for chain in egfr_structure[0]:
#     for residue in chain:
#         if residue.get_resname() == 'ASP' and residue.get_id()[1] == 846:
#             active_site_chain = chain
#             break
#     if active_site_chain:
#         break


# # Add EGF to active site of EGFR
# egfr_model = egfr_structure[0]
# for atom in egf_chain.get_atoms():
#     atom_name = atom.get_name()
#     if is_aa(atom.get_parent()):
#         residue_id = atom.get_parent().get_id()
#         active_site_residue = active_site_chain[residue_id]
#         active_site_residue.add(atom)

# # Calculate distance between EGF and active site residues
# active_site_atoms = list(active_site_chain.get_atoms())
# egf_atoms = list(egf_chain.get_atoms())
# for active_site_atom in active_site_atoms:
#     for egf_atom in egf_atoms:
#         distance = active_site_atom - egf_atom
#         if distance < 5:
#             print(
#                 f'Distance between {active_site_atom.get_full_id()} and {egf_atom.get_full_id()} is {distance:.2f} Å')


# # Generate PDB files for EGF-EGFR complex
# io = MMCIF2Dict('pdb4hjo.ent')
# io['atom_site']['pdbx_PDB_model_num'] = 1
# with open('egf-egfr.pdb', 'w') as f:
#     f.write(io.as_pdb_string())

# # Visualize EGF-EGFR complex with nglview
# view = nglview.show_structure_file('egf-egfr.pdb')
# view.add_representation('licorice', selection='protein')
# view.add_representation('ball+stick', selection='LIG')
# view.add_surface(opacity=0.3, color='white',
#                  wireframe=False, selection='protein')
# view.add_distance(atom_pair=[[f'{active_site_chain.get_id()} and name OD1', f'{egf_chain.get_id()} and name OE1'], [f'{active_site_chain.get_id()} and name OD2', f'{egf_chain.get_id()} and name OE1'], [
#                   f'{active_site_chain.get_id()} and name OD2', f'{egf_chain.get_id()} and name OE2'], [f'{active_site_chain.get_id()} and name OD1', f'{egf_chain.get_id()} and name OE2']], label_visible=True, color='green', radius=0.1)
# view.center()
# view.gui()
