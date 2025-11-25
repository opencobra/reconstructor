function loadMetaboliteInfoDiv(reactionData) {
    const metaboliteInfoDiv = document.getElementById('metaboliteinfo-div');
    if (!metaboliteInfoDiv) {
        console.error('Metabolite info div not found.');
        return;
    }

    let contentContainer = metaboliteInfoDiv.querySelector('.metabolite-panel-content');
    if (!contentContainer) {
        contentContainer = document.createElement('div');
        contentContainer.className = 'panel-content metabolite-panel-content';
        metaboliteInfoDiv.appendChild(contentContainer);
    }
    contentContainer.innerHTML = '';
    fillMetaboliteInfoTab(reactionData, contentContainer);
}

function fillMetaboliteInfoTab(data, targetContainer) {
    const metaboliteInfoContainer =
        targetContainer || document.querySelector('#metaboliteinfo-div .metabolite-panel-content') || document.getElementById('metaboliteinfo-div');
    if (!metaboliteInfoContainer) {
        return;
    }
    metaboliteInfoContainer.innerHTML = '';

    // Add the legend at the top
    const legend = document.createElement('div');
    legend.className = 'color-legend';
    legend.innerHTML = `
        <strong>Color Legend:</strong>
        <div class="legend-content">
            <span style="color: #C8C8C8;">■</span> C
            <span style="color: #FFFFFF; margin-left: 1em;">■</span> H
            <span style="color: #FF0D0D; margin-left: 1em;">■</span> O
            <span style="color: #3050F8; margin-left: 1em;">■</span> N
            <span style="color: #FFFF30; margin-left: 1em;">■</span> S
            <span style="color: #FF00FF; margin-left: 1em;">■</span> Unspecified Stereo
        </div>
    `;
    metaboliteInfoContainer.appendChild(legend);

    // Iterate through metabolite data
    data.metabolite_names.forEach((name, index) => {
        const metaboliteDiv = document.createElement('div');
        metaboliteDiv.classList.add('metabolite');

        const toggleDiv = document.createElement('div');
        toggleDiv.classList.add('metabolite-header');

        const nameElement = document.createElement('h3');
        nameElement.textContent = name;

        const toggleButton = document.createElement('button');
        toggleButton.textContent = 'Show 3D Structure';
        toggleButton.classList.add('toggle-button');
        toggleButton.onclick = function() {
            const structureContainer = this.parentNode.parentNode.querySelector('.structure-container');
            if (structureContainer.style.display === 'none') {
                structureContainer.style.display = 'block';
                this.textContent = 'Hide 3D Structure';

                if (!structureContainer.hasAttribute('data-viewer-initialized')) {
                    structureContainer.style.height = '400px';
                    structureContainer.style.width = '400px';
                    structureContainer.style.position = 'relative';

                    let viewer = new $3Dmol.createViewer(structureContainer, { backgroundColor: 'white' });
                    let molecularData = data.metabolite_mol_file_strings[index];

                    let model = viewer.addModel(molecularData, 'sdf');
                    viewer.setStyle({}, {
                        stick: { radius: 0.15, colorscheme: 'Jmol' },
                        sphere: { scale: 0.25, colorscheme: 'Jmol' }
                    });

                    if (data.stereo_locations_list) {
                        let stereoLocations = data.stereo_locations_list[index];
                        for (const loc of stereoLocations) {
                            viewer.setStyle({ model: model, index: loc }, { stick: { color: 'magenta' } });
                        }
                    }
                    viewer.setClickable({}, true, function(atom) {
                        viewer.addLabel(atom.atom, { position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8 });
                    });
                    viewer.zoomTo();
                    viewer.render();
                    structureContainer.setAttribute('data-viewer-initialized', 'true');
                }
            } else {
                structureContainer.style.display = 'none';
                this.textContent = 'Show 3D Structure';
            }
        };

        toggleDiv.appendChild(nameElement);
        toggleDiv.appendChild(toggleButton);
        metaboliteDiv.appendChild(toggleDiv);

        // **Compact Info Table**
        const infoTable = document.createElement('div');
        infoTable.classList.add('metabolite-info-grid');
        infoTable.innerHTML = `
            <div class="info-item"><strong>Charged Formula:</strong> ${data.metabolite_formulas[index] || 'N/A'}</div>
            <div class="info-item"><strong>SMILES:</strong> ${data.metabolite_smiles[index] || 'N/A'}</div>
            <div class="info-item"><strong>InChI:</strong> ${data.metabolite_inchis[index] || 'N/A'}</div>
            <div class="info-item"><strong>InChI Key:</strong> ${data.metabolite_inchi_keys[index] || 'N/A'}</div>
            <div class="info-item"><strong>Molecular Weight:</strong> ${data.metabolite_mol_weights[index] || 'N/A'} g/mol</div>
        `;
        metaboliteDiv.appendChild(infoTable);

        // Stereo count (if exists)
        if (data.stereo_counts || data.stereo_locations_list) {
            const stereoCount = data.stereo_counts[index];
            const stereoCountElement = document.createElement('p');
            stereoCountElement.textContent = stereoCount > 0
                ? `Number of unspecified Stereo Centers: ${stereoCount} (magenta in 3D viewer)`
                : 'No unspecified stereo centers detected.';
            metaboliteDiv.appendChild(stereoCountElement);
        }

        const structureContainer = document.createElement('div');
        structureContainer.className = 'structure-container';
        structureContainer.style.display = 'none';
        metaboliteDiv.appendChild(structureContainer);

        metaboliteInfoContainer.appendChild(metaboliteDiv);
    });
}



function toggleStructure() {
    const buttons = document.querySelectorAll('.toggle-button');
    buttons.forEach(button => {
        button.addEventListener('click', () => {
            const messageContainer = button.parentNode.nextElementSibling;
            if (messageContainer.style.display === 'none' || messageContainer.style.display === '') {
                messageContainer.style.display = 'block';
                button.textContent = 'Hide Message';
            } else {
                messageContainer.style.display = 'none';
                button.textContent = 'Show Message';
            }
        });
    });
}
