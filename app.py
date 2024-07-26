from flask import Flask, request, jsonify
from flask_cors import CORS
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

app = Flask(__name__)
CORS(app)

# Load the trained model and other objects
best_model = joblib.load('best_model_RandomForest.pkl')
label_encoder_smiles = joblib.load('label_encoder_smiles.pkl')
scaler = joblib.load('scaler.pkl')
imputer = joblib.load('imputer.pkl')
feature_columns = joblib.load('feature_columns.pkl')

# List all available molecular descriptors
descriptor_names = [desc_name[0] for desc_name in Descriptors._descList]
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return calculator.CalcDescriptors(mol)
    else:
        return [np.nan] * len(descriptor_names)

def calculate_morgan_fingerprints(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return np.array(fp)
    else:
        return np.array([np.nan] * n_bits)

def generate_best_fragment(smiles, model, feature_columns):
    fingerprints = calculate_morgan_fingerprints(smiles)
    descriptors = calculate_descriptors(smiles)
    if any(np.isnan(fingerprints)) or any(np.isnan(descriptors)):
        return "Invalid SMILES input"
    
    input_data = pd.DataFrame([list(fingerprints) + list(descriptors)], columns=feature_columns)
    
    for col in feature_columns:
        if col not in input_data.columns:
            input_data[col] = np.nan
    input_data = input_data[feature_columns]
    
    input_data = imputer.transform(input_data)
    input_data = scaler.transform(input_data)
    
    fragment_label = model.predict(input_data)[0]
    fragment_smiles = label_encoder_smiles.inverse_transform([int(fragment_label)])[0]
    
    return fragment_smiles

def cleanup_molecule_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)

def smiles_to_pdb(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, nonBondedThresh=500.0)
    pdb_block = Chem.MolToPDBBlock(mol)
    return pdb_block

@app.route('/get_3d_structure', methods=['POST'])
def get_3d_structure():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, nonBondedThresh=500.0)
    pdb_block = Chem.MolToPDBBlock(mol)

    return jsonify({"pdb": pdb_block})

@app.route('/predict_fragment', methods=['POST'])
def predict_fragment():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    fragment_smiles = generate_best_fragment(smiles, best_model, feature_columns)
    cleaned_fragment_smiles = cleanup_molecule_rdkit(fragment_smiles)
    fragment_pdb = smiles_to_pdb(cleaned_fragment_smiles)

    if fragment_pdb is None:
        return jsonify({"error": "Failed to generate PDB for fragment"}), 500

    return jsonify({"fragment_smiles": cleaned_fragment_smiles, "pdb": fragment_pdb})

if __name__ == '__main__':
    app.run(debug=True, port=5000)
