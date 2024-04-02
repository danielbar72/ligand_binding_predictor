import subprocess
from collections import defaultdict

def run_fpocket(pdb_file):
    # Run fpocket command
    command = f"fpocket -f {pdb_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Check for errors
    if process.returncode != 0:
        raise Exception(f"fpocket failed with error: {stderr.decode()}")

    # Process the output
    pockets = defaultdict(list)
    lines = stdout.decode().split("\n")
    current_pocket_id = None
    for line in lines:
        if line.startswith("Pocket"):
            current_pocket_id = int(line.split(":")[1].strip())
        elif line.startswith("Residue"):
            residue_info = line.split(":")[1].strip()
            residue_chain, residue_id = residue_info.split("_")
            pockets[current_pocket_id].append((residue_chain, residue_id))

    return pockets

# Example usage
pdb_file = "path/to/your/pdb/file.pdb"
pockets = run_fpocket(pdb_file)
print(f"Identified {len(pockets)} surface pockets:")
for pocket_id, residues in pockets.items():
    print(f"Pocket {pocket_id}: Residues {residues}")
