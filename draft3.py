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

def get_resname(residue):
    return residue.resname + str(residue.id[1])

def extract_features(pdb_file):
    # Run fpocket to get pockets
    pockets = find_pockets(pdb_file)

    # Get residues within distance of ligand
    residues_near_ligand = residues_within_distance(pdb_file, 6)

    # Calculate solvent accessibility
    #solvent_accessibility = extract_solvent_accessibility(pdb_file)

    # Calculate local density
    local_densities = calculate_local_density(pdb_file)

    df = pd.DataFrame(columns=['name', 'pocket_number', 'local_density', 'hydrophobicity' ,'is_ligand_binding'])
    structure = load_pdb_file(pdb_file)

    amino_acid_hydrophobicity = {
    'ALA': 1.8,
    'ILE': 4.5,
    'LEU': 3.8,
    'VAL': 4.2,
    'PHE': 2.8,
    'TRP': -0.9,
    'MET': 1.9,
    'PRO': -1.6,
    'TYR': -1.3,
    'CYS': 2.5,
    'ASN': -3.5,
    'GLN': -3.5,
    'SER': -0.8,
    'THR': -0.7,
    'ASP': -3.5,
    'GLU': -3.5,
    'HIS': -3.2,
    'LYS': -3.9,
    'ARG': -4.5,
    'GLY': -0.4
    }

    for r in structure.get_residues():
        if not is_aa(r):
            continue

        name = get_resname(r)
        pocket_number = 0
        local_density = local_densities[name]
        is_ligand_binding = name in residues_near_ligand
        hydrophobicity = amino_acid_hydrophobicity[r.resname]
        for id, residues in pockets.items():
            if name in residues:
                pocket_number = id
                break
        df = df._append({'name': name, 'pocket_number': pocket_number, 'local_density': local_density,'hydrophobicity': hydrophobicity,
                          'is_ligand_binding': is_ligand_binding}, ignore_index=True)

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
    combined_df = pd.DataFrame(columns=['name', 'pocket_number', 'local_density', 'hydrophobicity' ,'is_ligand_binding'])
    for pdb_file in pdb_files:
        df = extract_features(pdb_file)
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    return combined_df

def train_model(X_train, y_train, X_test, y_test):
    rf_classifier = RandomForestClassifier(n_estimators=27)

    rf_classifier.fit(X_train, y_train)


    y_pred = rf_classifier.predict(X_test)

    df_results = pd.DataFrame({'y_test': y_test, 'y_pred': y_pred})


    print(df_results[df_results['y_test'] == 1])

    print("______________")
    print(df_results[df_results['y_pred'] == 1])


    accuracy = accuracy_score(y_test, y_pred)

    print("Accuracy:", accuracy)

    return rf_classifier

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
    pdb_files = ["ligand_protein_pdb/2bqv.pdb"]  # Example: list of known binding site PDB files
    df = prepare_dataset(pdb_files)

    # True/False needs to be as 1 or 0
    df['is_ligand_binding'] = df['is_ligand_binding'].astype(int)
    X = df[['pocket_number', 'local_density', 'hydrophobicity']]
    y = df['is_ligand_binding']

    # Splitting the dataset into training and testing sets (80% train, 20% test)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

    # Train the machine learning model
    model = train_model(X_train, y_train, X_test, y_test)

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
