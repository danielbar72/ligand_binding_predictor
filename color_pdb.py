from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Color import Color

def color_residues(pdb_input_file):
    # Load the PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_input_file)
    highlighted_residues = [] 
    # Define the residues to be highlighted
    predicted_sites = predict_binding_sites(pdb_input_file, model)
    for pocket_num, residues in predicted_sites.items():
        if pocket_num == 100:
            for r in residues:
                highlighted_residues.append(r[1])


    # Color the residues
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] in highlighted_residues:
                    residue.color = Color(1, 0, 0)  # Red

    # Save the modified PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save("output.pdb")
