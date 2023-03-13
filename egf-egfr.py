import os
import urllib.request
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqUtils import seq1
from Bio import SeqIO
import nglview


# Load the CIF file for EGF and EGFR
egf_structure = MMCIFParser().get_structure('EGF', 'eg/1egf.cif')
egfr_structure = MMCIFParser().get_structure('EGFR', 'hj/4hjo.cif')


# Download the CIF file for EGF and EGFR
urllib.request.urlretrieve(
    "https://files.rcsb.org/download/1EGF.cif", "1EGF.cif")
urllib.request.urlretrieve(
    "https://files.rcsb.org/download/4HJO.cif", "4HJO.cif")


# Load the UniProt sequence for EGFR
fasta_file = os.path.expanduser('P00533.fasta')
handle = open(fasta_file)
seq_record = SeqIO.read(handle, 'fasta')
handle.close()

# Extract the sequence of the kinase domain of EGFR (residues 711-980)
kinase_domain = seq_record.seq[710:979]

# Find the chain that contains the kinase domain of EGFR in the CIF structure
kinase_domain_chain = None
for chain in egfr_structure.get_chains():
    for residue in chain.get_residues():
        if residue.get_resname() == 'TPO':
            if chain.id in residue.get_full_id():
                kinase_domain_chain = chain
                break
    if kinase_domain_chain is not None:
        break

# Print the sequence of the kinase domain in the CIF file
if kinase_domain_chain is not None:
    cif_seq = ''
    for residue in kinase_domain_chain.get_residues():
        if residue.get_resname() == 'TPO':
            cif_seq += 'T'
        else:
            cif_seq += seq1(residue.get_resname())
    print(cif_seq)
else:
    print('Kinase domain chain not found in CIF file')

# Visualize EGF-EGFR complex with nglview
view = nglview.show_structure_file('4HJO.cif')
view.add_component('1EGF.cif')
view.add_representation('cartoon', selection='protein')
view.add_representation('spacefill', selection='EGF', opacity=0.5)
view.add_representation('surface', selection='EGF',
                        opacity=0.3, color='yellow')
view.add_representation(
    'spacefill', selection='chain A and resnum 708-980', opacity=0.5)
view.add_representation(
    'line', selection='chain A and resnum 708-980', color='green')
view.add_representation(
    'surface', selection='chain A and resnum 708-980', opacity=0.3, color='purple')
view.center()
view.camera = 'perspective'
view.background = 'white'
view.update_representation()
view
