import subprocess

def run_fpocket(pdb_file):
    # Run fpocket command
    command = f"fpocket -f {pdb_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Check for errors
    if process.returncode != 0:
        raise Exception(f"fpocket failed with error: {stderr.decode()}")

    # Process the output
    pockets = []
    lines = stdout.decode().split("\n")
    for line in lines:
        if line.startswith("Pocket"):
            pocket_id = int(line.split(":")[1].strip())
            pockets.append(pocket_id)

    return pockets

# Example usage
pdb_file = "path/to/your/pdb/file.pdb"
pockets = run_fpocket(pdb_file)
print(f"Identified {len(pockets)} surface pockets: {pockets}")

# returns a list of pocket IDs identified by fpocket in the PDB file.