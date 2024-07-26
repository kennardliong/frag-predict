from flask import Flask, request, jsonify
from flask_cors import CORS  # Import CORS
from rdkit import Chem
from rdkit.Chem import AllChem

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

@app.route('/get_3d_structure', methods=['POST'])
def get_3d_structure():
    data = request.json
    smiles = data.get('smiles')

    if not smiles:
        return jsonify({"error": "SMILES string is required"}), 400

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return jsonify({"error": "Invalid SMILES string"}), 400

    # Add hydrogens and generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, nonBondedThresh=500.0)

    # Convert to PDB format
    pdb_block = Chem.MolToPDBBlock(mol)

    return jsonify({"pdb": pdb_block})

if __name__ == '__main__':
    app.run(debug=True, port=5000)  # Ensure Flask is running on port 5000
