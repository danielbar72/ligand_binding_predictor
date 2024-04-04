import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from accessibility import extract_solvent_accessibility
from local_densities import calculate_local_density
from pockets import find_pockets
from neighboring_final import residues_within_distance
import subprocess
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

def load_pdb_file(file_path):
    # Code to load and process the PDB file
    # For example, you can use the Biopython library to parse the PDB file
    from Bio.PDB import PDBParser

    parser = PDBParser()
    structure = parser.get_structure("pdb", file_path)

    # Process the structure and extract relevant features
    # For example, you can extract the coordinates of atoms or residues

    return structure

def extract_features(pdb_file):
    # Run fpocket to get pockets
    pockets = find_pockets(pdb_file)

    # Get residues within distance of ligand
    #residues_near_ligand = residues_within_distance(pdb_file, ligand_distance)

    # Calculate solvent accessibility
    #solvent_accessibility = extract_solvent_accessibility(pdb_file)

    # Calculate local density
    #local_density = calculate_local_density(pdb_file)

    df = pd.DataFrame(columns=['name', 'pocket_num'])  # You can specify column names here

    structure = load_pdb_file(pdb_file)

    print(len(list(structure.get_residues())))


    for r in structure.get_residues():
        if not is_aa(r):
            continue

        name = r.resname + str(r.id[1])
        was_in_pocket = False
        for id, residues in pockets.items():
            if name in residues:
                df = df._append({'name': name, 'pocket_num': id}, ignore_index=True)
                was_in_pocket = True
                break
        if not was_in_pocket:
            df = df._append({'name': name, 'pocket_num': 0}, ignore_index=True)

    return df
        


    

# Extract labels from the PDB file
def extract_labels(pdb_file, ligand_distance=5):
    # Get residues within distance of ligand
    residues_near_ligand = residues_within_distance(pdb_file, ligand_distance)

    # Extract labels for each residue
    labels = {}
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.id[1]
                residue_chain = chain.id
                label = 1 if residue in residues_near_ligand else 0
                labels[(residue_chain, residue_id)] = label

    return labels


def prepare_dataset(pdb_files):
    features = []
    labels = []
    for pdb_file in pdb_files:
        # Load PDB file and extract features
        structure = load_pdb_file(pdb_file)
        protein_features = extract_features(pdb_file)
        features.append(protein_features)

        # Extract labels from PDB file
        label = extract_labels(pdb_file)
        labels.append(label)

    # Convert features and labels lists to numpy arrays
    X = np.array(features)
    y = np.array(labels)

    return X, y

def train_model(X, y):
    # Code to train the model using the extracted features and labels
    model = RandomForestClassifier()  # Example: using RandomForestClassifier
    model.fit(X, y)
    return model

def predict_binding_sites(pdb_file, model):
    # Load PDB file and extract features
    structure = load_pdb_file(pdb_file)
    features = extract_features(structure)

    # Predict binding sites using the trained model
    predicted_labels = model.predict(features)

    # Extract residues predicted as binding sites
    predicted_sites = []
    for residue, label in zip(structure, predicted_labels):
        residue_id = residue.id[1]
        residue_chain = residue.parent.id
        if label == 1:
            predicted_sites.append((residue_chain, residue_id))

    # Print predicted sites
    print("Predicted Sites:")
    for i, site in enumerate(predicted_sites, start=1):
        print(f"Site {i} Residues:")
        for residue_chain, residue_id in site:
            print(f"{residue_chain}{residue_id}")

    return predicted_sites

# Main function
def main():

    # Prepare the dataset
    pdb_files = ["ligand_protein_pdb/4yay.pdb"]  # Example: list of known binding site PDB files
    X, y = prepare_dataset(pdb_files)

    # Train the machine learning model
    model = train_model(X, y)

    # Predict the binding sites
    #predicted_sites = predict_binding_sites(pdb_file, model)

    # Print the predicted binding sites
    """
    print("Predicted Binding Sites:")
    for site in predicted_sites:
        print(site)
    """
if __name__ == "__main__":
    main()
