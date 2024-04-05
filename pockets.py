import os
import subprocess

def find_pockets(pdb_file):
    # Run fpocket command to generate pocket files
    subprocess.run(['fpocket', '-f', pdb_file])

    # Get the path to the fpocket output directory
    output_dir = os.path.splitext(pdb_file)[0] + '_out'

    # Get the path to the pockets directory
    pockets_dir = os.path.join(output_dir, 'pockets')

    pocket_residues = {} 

    # Iterate over the pocket files in the pockets directory
    for file_name in os.listdir(pockets_dir):
        if file_name.endswith('.pdb'):
            pocket_num = int(file_name.split('.')[0].split('_')[0][6:])
            
            # Read the pocket file and extract the residues
            with open(os.path.join(pockets_dir, file_name), 'r') as pocket_file:
                for line in pocket_file:
                    if line.startswith('ATOM'):
                        residue = line[17:20].strip()
                        residue_id = line[22:26].strip()
                        full_residue = f"{residue}{residue_id}"
                        if pocket_num not in pocket_residues:
                            pocket_residues[pocket_num] = [full_residue]
                        else:
                            
                            pocket_residues[pocket_num].append(full_residue)

    return pocket_residues

