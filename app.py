from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from frag_predict import generate_best_fragment, cleanup_molecule_rdkit, calculate_properties, get_3d_structure, best_model, feature_columns
import os

app = Flask(__name__)
CORS(app)

@app.route('/get_3d_structure', methods=['POST'])
def get_3d_structure_route():
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

@app.route('/get_2d_structure', methods=['POST'])
def get_2d_structure_route():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    img = Draw.MolToImage(mol)
    img_path = os.path.join(os.getcwd(), 'molecule.png')
    img.save(img_path)

    return send_file(img_path, mimetype='image/png')

@app.route('/predict_fragment', methods=['POST'])
def predict_fragment():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    fragment_smiles = generate_best_fragment(smiles, best_model, feature_columns)
    cleaned_fragment_smiles = cleanup_molecule_rdkit(fragment_smiles)
    
    if not cleaned_fragment_smiles:
        return jsonify({"error": "Failed to generate a valid fragment"}), 500
    
    fragment_pdb = get_3d_structure(cleaned_fragment_smiles)
    properties = calculate_properties(cleaned_fragment_smiles)

    if not fragment_pdb:
        return jsonify({"error": "Failed to generate PDB for fragment"}), 500

    return jsonify({
        "fragment_smiles": cleaned_fragment_smiles,
        "pdb": fragment_pdb,
        "properties": properties
    })

@app.route('/download_pdb', methods=['POST'])
def download_pdb():
    data = request.json
    smiles = data.get('smiles')
    filename = data.get('filename', 'structure.pdb')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, nonBondedThresh=500.0)
    pdb_block = Chem.MolToPDBBlock(mol)
    file_path = os.path.join(os.getcwd(), filename)

    with open(file_path, 'w') as file:
        file.write(pdb_block)

    return send_file(file_path, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True, port=5000)
