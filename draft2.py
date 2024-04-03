import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import neighboring_final
from Bio.PDB import PDBParser, is_aa

# Load the PDB file and extract relevant features
def load_pdb_file(file_path):
    # Code to load and process the PDB file
    pass

# Prepare the dataset for training
def prepare_dataset(pdb_file):

    parser = PDBParser()

# Parse the PDB file
    structure = parser.get_structure('1XYZ', pdb_file)

    
    residues_within_ligand = [r.id for r in neighboring_final.residues_within_distance(pdb_file, 6)]


    data = {
        'X': [],
        'Y': [],
        'Z': [],
        'is_ligand_binding': []
    }

    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue):
                    continue
                for atom in residue:
                    data['X'].append(atom.get_coord()[0])
                    data['Y'].append(atom.get_coord()[1])
                    data['Z'].append(atom.get_coord()[2])

                    if residue.id in residues_within_ligand:
                        data['is_ligand_binding'].append(1)
                    else:
                        data['is_ligand_binding'].append(0)


    df = pd.DataFrame(data)
    return df



# Train the machine learning model
def train_model(X_train, y_train, X_test, y_test):

# Initialize the Random Forest Classifier model
    rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
    rf_classifier.fit(X_train, y_train)

# Predict on the testing set
    y_pred = rf_classifier.predict(X_test)

    df_results = pd.DataFrame({'y_train': y_test, 'y_test': y_pred})

# Print the concatenated DataFrame
    print(df_results)
# Evaluate the model
    accuracy = accuracy_score(y_test, y_pred)

    print("Accuracy:", accuracy)

    return rf_classifier

def predict_binding_sites(pdb_file, model):
    # Code to predict the binding sites using the trained model
    pass





# Main function
def main():



    pdb_file = "ligand_protein_pdb/4yay.pdb"

    # Data frame
    df = prepare_dataset(pdb_file)

    X = df[['X', 'Y', 'Z']]
    y = df['is_ligand_binding']

    # Splitting the dataset into training and testing sets (80% train, 20% test)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train the machine learning model
    model = train_model(X_train, y_train, X_test, y_test)


    # Try the ML model on another protein
    df_val = prepare_dataset("ligand_protein_pdb/7ofk.pdb")

    X_val = df_val[['X', 'Y', 'Z']]
    y_val = df_val['is_ligand_binding']

    y_pred = model.predict(X_val)

    accuracy = accuracy_score(y_val, y_pred)

    print("Accuracy:", accuracy)

    '''
    # Predict the binding sites
    predicted_sites = predict_binding_sites(pdb_file, model)

    # Print the predicted binding sites
    print("Predicted Binding Sites:")
    for site in predicted_sites:
        print(site)

'''


if __name__ == "__main__":
    main()