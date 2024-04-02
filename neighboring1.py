from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom

def residues_within_distance(pdb_file, ligand_distance):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)

    # Extract ligand atoms
    ligand_atoms = [atom for atom in structure.get_atoms() if isinstance(atom, Atom) and atom.get_id() == 'HETATM']

    # Extract residues within specified distance of ligand
    residues_within_distance = set()
    for atom in ligand_atoms:
        for residue in structure.get_residues():
            for residue_atom in residue.get_atoms():
                if atom - residue_atom <= ligand_distance:
                    residues_within_distance.add(residue)

    return residues_within_distance

# Example usage
pdb_file = "path/to/your/pdb/file.pdb"
ligand_distance = 5.0  # Specify the distance threshold
residues = residues_within_distance(pdb_file, ligand_distance)
print("Residues within distance of ligand:")
for residue in residues:
    print(f"Chain: {residue.get_full_id()[2]}, Residue: {residue.get_resname()}{residue.id[1]}")
