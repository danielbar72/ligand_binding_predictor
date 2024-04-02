from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def extract_solvent_accessibility(pdb_file):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)

    # Calculate solvent accessibility using DSSP
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    # Extract solvent accessibility for each residue
    solvent_accessibility = {}
    for residue in model.get_residues():
        residue_id = residue.get_id()
        chain_id = residue_id[2]
        residue_number = residue_id[1]
        key = (chain_id, residue_number)
        solvent_accessibility[key] = dssp[key][3]

    return solvent_accessibility

# Usage example
pdb_file = "path/to/your/pdb/file.pdb"
solvent_accessibility = extract_solvent_accessibility(pdb_file)
print(solvent_accessibility)

# returns a dictionary where the keys are tuples (chain_id, residue_number) 
# and the values are the corresponding solvent accessibilities.


