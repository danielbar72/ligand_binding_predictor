from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Atom import Atom


def get_ligand_id(pdb_file):
    # Load PDB file
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure("pdb", pdb_file)

    # List of common ligand residue names, excluding water
    ligand_names = {'H_'}  # Add other ligand residue names as needed

    for model in structure:
        for chain in model:
            for residue in chain:
                residue_name = residue.resname.strip()
                if residue.id[0][0:2] == 'H_' and residue_name not in {'HOH', 'WAT'}:
                    if len(residue) > 10:  # Check for large HETATM
                        return (residue.id[0])



def residues_within_distance(pdb_file, ligand_distance):
    # Parse the PDB file
    # Get ligand id
    ligand_id = get_ligand_id(pdb_file)


    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)

    ligand_atoms = []

    # Get ligand atoms
    for res in structure.get_residues():
        if res.id[0] == ligand_id:
            for atom in res:
                ligand_atoms.append(atom)

    # Extract residues within specified distance of ligand, consider only amino acids
    residues_within_distance = set()
    for atom in ligand_atoms:
        for residue in structure.get_residues():
            if not is_aa(residue):
                    continue
            for residue_atom in residue.get_atoms():
                if atom - residue_atom <= ligand_distance:
                    name = residue.resname + str(residue.id[1])
                    residues_within_distance.add(name)

    return residues_within_distance

"""
# Example usage
pdb_file = "ligand_protein_pdb/4yay.pdb"
ligand_distance = 6 # Specify the distance threshold
residues = residues_within_distance(pdb_file, ligand_distance)
print("Residues within distance of ligand:")
for residue in residues:
    print(f"Chain: {residue.get_full_id()[2]}, Residue: {residue.get_resname()}{residue.id[1]}")
"""
