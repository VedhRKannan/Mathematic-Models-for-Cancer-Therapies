from Bio.PDB import *


parser = MMCIFParser()
structure = parser.get_structure("1NQL", "1nql.pdb")


def get_chain_ids(structure):
    chain_ids = []
    for model in structure:
        for chain in model:
            if "EGF" in chain.get_full_id()[3]:
                chain_ids.append(chain.get_id())
            elif "EGFR" in chain.get_full_id()[3]:
                chain_ids.append(chain.get_id())
    return chain_ids


def get_binding_domain(structure):
    chain_ids = get_chain_ids(structure)
    binding_residues = []
    for residue in structure[0][chain_ids[0]]:
        if residue.get_resname() == "EGF":
            for atom in residue:
                if "CA" in atom.get_id():
                    egf_ca_coord = atom.get_coord()
        elif residue.get_resname() == "EGFR":
            for atom in residue:
                if "CA" in atom.get_id():
                    egfr_ca_coord = atom.get_coord()
            if abs(egf_ca_coord - egfr_ca_coord) < 15.0:
                binding_residues.append(residue.get_id()[1])
    return binding_residues


binding_residues = get_binding_domain(structure)
print("Binding domain residues:", binding_residues)
