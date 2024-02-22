from Bio.PDB import PDBParser
from Bio.PDB import Selection

# Load PDB file
pdb_parser = PDBParser()
structure = pdb_parser.get_structure("PHA-L", "1fat.pdb")



for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom)


# Extract all residues with a ligand
ligand_residues = Selection.unfold_entities(structure, "R")

# You can then perform further analysis or prediction on the ligand residues
# For example, identifying ligand binding sites based on certain criteria
# or using machine learning algorithms

# Here's a simple example:
# Assume you want to find ligand residues within a certain distance threshold from a metal ion
metal_ion_coord = (x, y, z)  # Coordinates of the metal ion
distance_threshold = 5.0  # Threshold distance in angstroms

ligand_binding_sites = []
for residue in ligand_residues:
    for atom in residue:
        atom_coord = atom.get_coord()
        distance = ((atom_coord[0] - metal_ion_coord[0])**2 + 
                    (atom_coord[1] - metal_ion_coord[1])**2 + 
                    (atom_coord[2] - metal_ion_coord[2])**2)**0.5
        if distance <= distance_threshold:
            ligand_binding_sites.append(residue)

# Now ligand_binding_sites contains residues within the specified distance from the metal ion
# You can further analyze these residues as needed
