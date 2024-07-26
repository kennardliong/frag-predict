document.getElementById('inputForm').addEventListener('submit', function(event) {
    event.preventDefault();
    const smiles = document.getElementById('smilesInput').value;
    const protein = document.getElementById('proteinSelector').value;

    if (!protein) {
        alert('Please select a protein.');
        return;
    }

    console.log("SMILES String:", smiles);
    console.log("Selected Protein:", protein);

    document.getElementById('subtitle-inputted').style.display = 'block';
    document.getElementById('subtitle-fragmented').style.display = 'block';

    fetch('http://localhost:5000/get_3d_structure', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ smiles: smiles })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            alert(data.error);
            return;
        }

        const viewer = $3Dmol.createViewer("viewer", {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: 'black'
        });
        
        viewer.addModel(data.pdb, "pdb");
        viewer.setStyle({}, {stick: {colorscheme: 'Jmol'}});
        viewer.zoomTo();
        viewer.render();
    })
    .catch(error => {
        console.error('Error:', error);
    });

    fetch('http://localhost:5000/predict_fragment', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ smiles: smiles })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            alert(data.error);
            return;
        }

        const fragmentSmiles = data.fragment_smiles;
        const fragmentPdb = data.pdb;
        console.log("Best Fragment SMILES:", fragmentSmiles);
        document.getElementById('fragment-output').textContent = `Best Fragment SMILES: ${fragmentSmiles}`;

        const viewerFragmented = $3Dmol.createViewer("viewer-fragmented", {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: 'black'
        });

        viewerFragmented.addModel(fragmentPdb, "pdb");
        viewerFragmented.setStyle({}, {stick: {colorscheme: 'Jmol'}});
        viewerFragmented.zoomTo();
        viewerFragmented.render();
    })
    .catch(error => {
        console.error('Error:', error);
    });
});
