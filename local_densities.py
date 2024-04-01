from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

def calculate_local_density(pdb_file, residue_type):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    # Get all residues of the specified type
    residues = [residue for residue in structure.get_residues() if residue.get_resname() == residue_type]

    # Create a NeighborSearch object
    ns = NeighborSearch(list(structure.get_atoms()))

    # Calculate the local density for each residue
    local_densities = []
    for residue in residues:
        neighbors = ns.search(residue.get_coord(), 10.0)  # Adjust the distance cutoff as needed
        local_density = len(neighbors) / (4/3 * 3.14 * 10.0**3)  # Adjust the volume calculation as needed
        local_densities.append(local_density)

    return local_densities

# Example usage
pdb_file = "path/to/your.pdb"
residue_type = "ALA"
densities = calculate_local_density(pdb_file, residue_type)
print(densities)

# returns a list of local densities for each residue of the specified type in the PDB file.