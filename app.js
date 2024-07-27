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

        const inputDownloadLink = document.getElementById('input-download');
        inputDownloadLink.style.display = 'block';
        inputDownloadLink.onclick = function() {
            fetch('http://localhost:5000/download_pdb', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ smiles: smiles, filename: 'input_structure.pdb' })
            })
            .then(response => {
                if (response.ok) {
                    return response.blob();
                } else {
                    alert('Failed to download PDB file');
                }
            })
            .then(blob => {
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                a.download = 'input_structure.pdb';
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
            })
            .catch(error => {
                console.error('Error:', error);
            });
        };
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
        const properties = data.properties;

        console.log("Best Fragment SMILES:", fragmentSmiles);
        document.getElementById('fragment-properties').innerHTML = `
            <p></p>
            <p><span class="pink-title">Best Fragment SMILES: </span> ${fragmentSmiles}</p>   
            <p><span class="pink-title">Molecular Weight: </span>${properties.molecular_weight.toFixed(2)}</p>
            <p><span class="pink-title">LogP: </span>${properties.log_p.toFixed(2)}</p>
            <p><span class="pink-title">Hydrogen Bond Donors: </span>${properties.hydrogen_bond_donors}</p>
            <p><span class="pink-title">Hydrogen Bond Acceptors: </span>${properties.hydrogen_bond_acceptors}</p>
            <p><span class="pink-title">TPSA: </span> ${properties.tpsa.toFixed(2)}</p>
        `;

        const viewerFragmented = $3Dmol.createViewer("viewer-fragmented", {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: 'black'
        });

        viewerFragmented.addModel(fragmentPdb, "pdb");
        viewerFragmented.setStyle({}, {stick: {colorscheme: 'Jmol'}});
        viewerFragmented.zoomTo();
        viewerFragmented.render();

        const fragmentDownloadLink = document.getElementById('fragment-download');
        fragmentDownloadLink.style.display = 'block';
        fragmentDownloadLink.onclick = function() {
            fetch('http://localhost:5000/download_pdb', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ smiles: fragmentSmiles, filename: 'fragment_structure.pdb' })
            })
            .then(response => {
                if (response.ok) {
                    return response.blob();
                } else {
                    alert('Failed to download PDB file');
                }
            })
            .then(blob => {
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.style.display = 'none';
                a.href = url;
                a.download = 'fragment_structure.pdb';
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
            })
            .catch(error => {
                console.error('Error:', error);
            });
        };
    })
    .catch(error => {
        console.error('Error:', error);
    });
});
