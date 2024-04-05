from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Polypeptide import is_aa

def calculate_local_density(pdb_file):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    # Get all residues
    residues = structure.get_residues()

    # Create a NeighborSearch object
    ns = NeighborSearch(list(structure.get_atoms()))

    # Calculate the local density for each residue
    local_densities = {}
    
    for residue in residues:
        if not is_aa(residue):
            continue

        c_alfa = None
        for atom in residue:
            if atom.id == 'CA':
                c_alfa = atom
                break

        if c_alfa is None:
            continue
        
        neighbors = ns.search(c_alfa.get_coord(), 10.0)  # Adjust the distance cutoff as needed
        local_density = len(neighbors) / (4/3 * 3.14 * 10.0**3)  # Adjust the volume calculation as needed

        name = residue.resname + str(residue.id[1])
        local_densities[name] = len(neighbors)
        """
        if name not in local_densities:
            local_densities[name] = [local_density]
        else:
                            
            local_densities[name].append(full_residue)
        """

    return local_densities


# calculate_local_density("ligand_protein_pdb/4yay.pdb")
# returns a list of local densities for each residue in the PDB file.
