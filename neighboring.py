from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

def calculate_distance(atom1, atom2):
    return atom1 - atom2

def get_residues_within_distance(pdb_file, ligand_id, distance_threshold):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    model = structure[0]

    ligand_atoms = [atom for atom in model.get_atoms() if atom.get_parent().get_id() == ligand_id]
    residues_within_distance = []

    for residue in model.get_residues():
        if is_aa(residue):
            for atom in residue.get_atoms():
                for ligand_atom in ligand_atoms:
                    distance = calculate_distance(atom.get_coord(), ligand_atom.get_coord())
                    if distance <= distance_threshold:
                        residues_within_distance.append(residue)
                        break

    return residues_within_distance

# Example usage
pdb_file = 'path/to/your/pdb/file.pdb'
ligand_id = 'A'
distance_threshold = 5.0

residues = get_residues_within_distance(pdb_file, ligand_id, distance_threshold)
for residue in residues:
    print(residue.get_id())