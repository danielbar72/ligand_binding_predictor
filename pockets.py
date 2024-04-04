import os
import subprocess
from Bio.PDB import PDBParser  


def find_pockets(pdb_file):
    # Run fpocket command to generate pocket files
    subprocess.run(['fpocket', '-f', pdb_file])

    # Get the path to the fpocket output directory
    output_dir = os.path.splitext(pdb_file)[0] + '_out'

    # Iterate over the pocket files in the output directory
    for file_name in os.listdir(output_dir):
        if file_name.endswith('.pdb'):
            pocket_num = file_name.split('.')[0].split('_')[-1]
            pocket_residues = set()  # Use a set to avoid duplicate residues
            
            file_name = os.path.join(output_dir, file_name)

            parser = PDBParser()
            structure = parser.get_structure("file", file_name)
            for r in structure.get_residues():
                print(r.id)
                
            """
            # Read the pocket file and extract the residues
            with open(os.path.join(output_dir, file_name), 'r') as pocket_file:
                for line in pocket_file:
                    if line.startswith('ATOM'):
                        residue = line[17:20].strip()
                        residue_id = line[22:26].strip()
                        full_residue = f"{residue}{residue_id}"
                        pocket_residues.add(full_residue)  # Add residues to set
            """
            # Output the pocket number and the residues involved
            #print(f"Residues in pocket {pocket_num}: {', '.join(pocket_residues)}")
