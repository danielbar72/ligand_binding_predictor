import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from local_densities import calculate_local_density
from pockets import find_pockets
from neighboring_final import residues_within_distance
import os, sys, warnings
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

def load_pdb_file(file_path):

    parser = PDBParser()
    structure = parser.get_structure("pdb", file_path)

    return structure

def get_resname(residue):
    return residue.resname + str(residue.id[1])

def extract_features(pdb_file):
    # Run fpocket to get pockets
    pockets = find_pockets(pdb_file)

    # Get residues within distance of ligand
    residues_near_ligand = residues_within_distance(pdb_file, 6)

    # Calculate local density
    local_densities = calculate_local_density(pdb_file)

    df = pd.DataFrame(columns=['name', 'pocket_number', 'local_density', 'hydrophobicity' ,'is_ligand_binding'])
    structure = load_pdb_file(pdb_file)

    # Hydrophobicity of each AA
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

    # Prepare a dataframe
    for r in structure.get_residues():
        if not is_aa(r):
            continue

        name = get_resname(r)
        pocket_number = 0
        local_density = local_densities.get(name, 0)
        is_ligand_binding = name in residues_near_ligand
        hydrophobicity = amino_acid_hydrophobicity[r.resname]

        for id, residues in pockets.items():
            if name in residues:
                pocket_number = id
                break
        df = df._append({'name': name, 'pocket_number': pocket_number, 'local_density': local_density,'hydrophobicity': hydrophobicity,
                          'is_ligand_binding': is_ligand_binding}, ignore_index=True)

    return df


def prepare_dataset(pdb_files_dir):
    combined_df = pd.DataFrame(columns=['name', 'pocket_number', 'local_density', 'hydrophobicity' ,'is_ligand_binding'])
    for pdb_file in os.listdir(pdb_files_dir):
        print("Training model on a file: ", pdb_file)
        if not pdb_file.endswith('.pdb'):
            continue
        df = extract_features(os.path.join(pdb_files_dir,pdb_file))
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    return combined_df

def train_model(X_train, y_train, X_test, y_test):
    rf_classifier = RandomForestClassifier(n_estimators=27)

    rf_classifier.fit(X_train, y_train)

    y_pred = rf_classifier.predict(X_test)

    df_results = pd.DataFrame({'y_test': y_test, 'y_pred': y_pred})

    accuracy = accuracy_score(y_test, y_pred)

    print("Accuracy on the training data: ", accuracy)

    return rf_classifier

def predict_binding_sites(pdb_file, model):

    df = extract_features(pdb_file)

    # Predict binding sites using the trained model
    df['is_ligand_binding'] = df['is_ligand_binding'].astype(int)

    X = df[['pocket_number', 'local_density', 'hydrophobicity']]
    y = df['is_ligand_binding']


    predicted_labels = model.predict(X)

    df['ligand_binding_predicted'] = predicted_labels

    df = df[df['ligand_binding_predicted'] == 1][['name', 'pocket_number']]

    ret_dict = {}
    for _, row in df.iterrows():
        name = row['name']
        pocket_num = row['pocket_number'] if row['pocket_number'] != 0 else 100

        if pocket_num not in ret_dict:
            ret_dict[pocket_num] = [name]
        else:
            ret_dict[pocket_num].append(name)
    return {k: ret_dict[k] for k in sorted(ret_dict.keys())}

# Main function
def main(pdb_input_file):

    # Prepare the dataset
    pdb_files_dir = "train_pdb_data"
    df = prepare_dataset(pdb_files_dir)

    # True/False needs to be as 1 or 0
    df['is_ligand_binding'] = df['is_ligand_binding'].astype(int)
    X = df[['pocket_number', 'local_density', 'hydrophobicity']]
    y = df['is_ligand_binding']

    # Splitting the dataset into training and testing sets (80% train, 20% test)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

    # Train the machine learning model
    model = train_model(X_train, y_train, X_test, y_test)

    
    predicted_sites = predict_binding_sites(pdb_input_file, model)

    commands = ["select:"]
    print("Predicted Binding Sites:")
    for pocket_num, residues in predicted_sites.items():
        if pocket_num == 100:
            print("Residues, that are not part of any pocket: ")
        else:
            print("Binding site pocket number: ", pocket_num)
            print("Residues: ")
        for r in residues:
            print("     ", r)
            commands[0] += r[3:] + ","

    commands[0] = commands[0][:len(commands[0]) - 1]
    commands.append("\ncolor red sel")

    with open("chimera_cmds.cmd", "w") as file:
        for c in commands:
            file.write(c)
    
    

if __name__ == "__main__":

    # Ignore warnings
    warnings.filterwarnings("ignore")
    # Has to be exactly one cmd line argument
    if len(sys.argv) == 2:
        file_path = sys.argv[1]
        if not os.path.exists(file_path):
            print(f"The file '{file_path}' does not exist.")
        else:        
            main(sys.argv[1])
    else:
        print("Please provide exactly one pdb file as an input.")

