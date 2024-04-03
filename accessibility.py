from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def extract_solvent_accessibility(pdb_file):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_file)

    # Calculate solvent accessibility using DSSP
    model = structure[0]
    dssp = DSSP(model, pdb_file)

    # Extract solvent accessibility for each residue
    solvent_accessibility = {}
    for residue in model.get_residues():
        residue_name = residue.get_resname()
        residue_number = residue.id[1]
        full_identifier = f"{residue_name}{residue_number}"
        solvent_accessibility[full_identifier] = dssp[residue.get_full_id()][3]

    return solvent_accessibility

# Usage example
pdb_file = "path/to/your/pdb/file.pdb"
solvent_accessibility = extract_solvent_accessibility(pdb_file)
print(solvent_accessibility)

#returns dictionary with residue id as key(eg LEU123) and accessibility float as value
