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

    // Display submitted information

    // Show subtitles
    document.getElementById('subtitle-inputted').style.display = 'block';
    document.getElementById('subtitle-fragmented').style.display = 'block';

    // Fetch 3D coordinates from backend
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

        // Initialize viewers
        const viewer = $3Dmol.createViewer("viewer", {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: 'black'
        });
        
        const viewerFragmented = $3Dmol.createViewer("viewer-fragmented", {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: 'black'
        });

        // Add models and render
        viewer.addModel(data.pdb, "pdb");
        viewer.setStyle({}, {stick: {colorscheme: 'Jmol'}});
        viewer.zoomTo();
        viewer.render();

        // You can use the same PDB data for now or adjust to use different data
        viewerFragmented.addModel(data.pdb, "pdb");
        viewerFragmented.setStyle({}, {stick: {colorscheme: 'Jmol'}});
        viewerFragmented.zoomTo();
        viewerFragmented.render();
    })
    .catch(error => {
        console.error('Error:', error);
    });
});
