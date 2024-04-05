import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from local_densities import calculate_local_density
from pockets import find_pockets
from neighboring_final import residues_within_distance
import os, sys
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBParser, PDBIO
import warnings

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

    # print("Accuracy:", accuracy)

    return rf_classifier

def predict_binding_sites(pdb_file, model):
    # Load PDB file and extract features
    structure = load_pdb_file(pdb_file)
    df = extract_features(pdb_file)

    # Predict binding sites using the trained model
    df['is_ligand_binding'] = df['is_ligand_binding'].astype(int)

    X = df[['pocket_number', 'local_density', 'hydrophobicity']]
    y = df['is_ligand_binding']

    

    predicted_labels = model.predict(X)

    df_results = pd.DataFrame({'y_test': y, 'y_pred': predicted_labels})

    print("Final model")
    print(df_results[df_results['y_test'] == 1])
    print("______________")
    print(df_results[df_results['y_pred'] == 1])

    accuracy = accuracy_score(y, predicted_labels)

    print("Accuracy:", accuracy)
    
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


from Bio.PDB import PDBParser, PDBIO
import matplotlib.colors as mpl_colors


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


# Main function
def main(pdb_input_file):

    # Prepare the dataset
    # pdb_files_dir = "ligand_protein_pdb"
    pdb_files_dir = "train"
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
    print("Predicted Binding Sites:")
    for pocket_num, residues in predicted_sites.items():
        if pocket_num == 100:
            print("Residues, that are not part of any pocket: ")
        else:
            print("Binding site pocket number: ", pocket_num)
            print("Residues: ")
        for r in residues:
            print("     ", r)
            
    output_file = color_residues(pdb_input_file, model)
    

if __name__ == "__main__":

    # Ignore warnings
    warnings.filterwarnings("ignore")
    if len(sys.argv) == 2:  # One argument is the script name itself, so we check if the count is 2
        file_path = sys.argv[1]
        if not os.path.exists(file_path):
            print(f"The file '{file_path}' does not exist.")
        else:        
            main(sys.argv[1])
    else:
        print("Please provide exactly one pdb file as an input.")

