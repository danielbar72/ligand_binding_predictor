import subprocess
from collections import defaultdict

def run_fpocket(pdb_file):
    try:
        # Run fpocket command
        result = subprocess.run(['fpocket', '-f', pdb_file], capture_output=True, text=True, check=True)

        # Process the output
        pockets = defaultdict(list)
        current_pocket_id = None
        for line in result.stdout.split("\n"):
            if line.startswith("Pocket"):
                current_pocket_id = int(line.split(":")[1].strip())
            elif line.startswith("Residue"):
                residue_info = line.split(":")[1].strip()
                residue_chain, residue_id = residue_info.split("_")
                pockets[current_pocket_id].append((residue_chain, residue_id))

        return pockets

    except subprocess.CalledProcessError as e:
        # fpocket command failed, raise custom exception with error message
        raise Exception(f"fpocket failed with error: {e.stderr}")

    except Exception as e:
        # Other exceptions (e.g., subprocess.run error), raise with original message
        raise e
