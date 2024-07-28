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

    // Show sections and buttons
    document.getElementById('subtitle-inputted').style.display = 'block';
    document.getElementById('subtitle-fragmented').style.display = 'block';
    document.getElementById('toggle-input-view').style.display = 'block';
    document.getElementById('toggle-fragment-view').style.display = 'block';

    // Fetch 2D image for input
    fetch('http://localhost:5000/get_2d_structure', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ smiles: smiles })
    })
    .then(response => response.blob())
    .then(blob => {
        const url = window.URL.createObjectURL(blob);
        document.getElementById('input-2d').src = url;
    })
    .catch(error => {
        console.error('Error:', error);
    });

    // Fetch 3D structure for input
    fetch('http://localhost:5000/get_3d_structure', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ smiles: smiles })
    })
    .then(response => response.json())
    .then(data => {
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
                }
                throw new Error('Network response was not ok.');
            })
            .then(blob => {
                const url = window.URL.createObjectURL(blob);
                inputDownloadLink.href = url;
            })
            .catch(error => {
                console.error('Error:', error);
            });
        };
    })
    .catch(error => {
        console.error('Error:', error);
    });

    // Predict fragment
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

        // Display 2D image for fragmented molecule
        fetch('http://localhost:5000/get_2d_structure', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ smiles: data.fragment_smiles })
        })
        .then(response => response.blob())
        .then(blob => {
            const url = window.URL.createObjectURL(blob);
            document.getElementById('fragment-2d').src = url;
        })
        .catch(error => {
            console.error('Error:', error);
        });

        // Display 3D structure for fragmented molecule
        const viewerFragmented = $3Dmol.createViewer("viewer-fragmented", {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: 'black'
        });

        viewerFragmented.addModel(data.pdb, "pdb");
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
                body: JSON.stringify({ smiles: data.fragment_smiles, filename: 'fragment_structure.pdb' })
            })
            .then(response => {
                if (response.ok) {
                    return response.blob();
                }
                throw new Error('Network response was not ok.');
            })
            .then(blob => {
                const url = window.URL.createObjectURL(blob);
                fragmentDownloadLink.href = url;
            })
            .catch(error => {
                console.error('Error:', error);
            });
        };

        // Display fragment properties in a formatted manner
        const propertiesDiv = document.getElementById('fragment-properties');
        propertiesDiv.innerHTML = `
            <h5 class="pink-title">Properties:</h5>
            <p></p>
            <p><strong>SMILES:</strong> ${data.fragment_smiles}</p>
            <p><strong>Molecular Weight:</strong> ${data.properties.molecular_weight} Da</p>
            <p><strong>LogP Value:</strong> ${data.properties.log_p}</p>
            <p><strong>Hydrogen Bond Acceptors:</strong> ${data.properties.hydrogen_bond_acceptors}</p>
            <p><strong>Hydrogen Bond Donors:</strong> ${data.properties.hydrogen_bond_donors}</p>
            <p><strong>Topological Polar Surface Area:</strong> ${data.properties.tpsa} Å²</p>
        `;
    })
    .catch(error => {
        console.error('Error:', error);
    });
});

// Toggle functionality
document.getElementById('toggle-input-view').addEventListener('click', function() {
    const viewer = document.getElementById('viewer');
    const img2d = document.getElementById('input-2d');
    const is2d = img2d.style.display === 'block';

    if (is2d) {
        img2d.style.display = 'none';
        viewer.style.display = 'block';
        this.textContent = 'Show 2D View';
    } else {
        img2d.style.display = 'block';
        viewer.style.display = 'none';
        this.textContent = 'Show 3D View';
    }
});

document.getElementById('toggle-fragment-view').addEventListener('click', function() {
    const viewer = document.getElementById('viewer-fragmented');
    const img2d = document.getElementById('fragment-2d');
    const is2d = img2d.style.display === 'block';

    if (is2d) {
        img2d.style.display = 'none';
        viewer.style.display = 'block';
        this.textContent = 'Show 2D View';
    } else {
        img2d.style.display = 'block';
        viewer.style.display = 'none';
        this.textContent = 'Show 3D View';
    }
});

function createSparkle(x, y) {
    const sparkle = document.createElement('div');
    sparkle.classList.add('sparkle');
    sparkle.style.left = `${x}px`;
    sparkle.style.top = `${y}px`;
    document.body.appendChild(sparkle);

    setTimeout(() => {
        sparkle.remove();
    }, 500); // Duration of the sparkle animation
}

// Add event listeners to all pink buttons
document.querySelectorAll('.btn-pink').forEach(button => {
    button.addEventListener('click', function(event) {
        // Calculate the sparkle coordinates relative to the entire document
        const x = event.pageX;
        const y = event.pageY;
        // Create the sparkle at the calculated coordinates
        createSparkle(x, y);
    });
});
