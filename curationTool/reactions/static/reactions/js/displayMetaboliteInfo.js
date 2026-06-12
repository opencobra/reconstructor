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

        // **Compact Info Table** — values truncate (…) and carry a copy button
        // that copies the full, untruncated text.
        const infoTable = document.createElement('div');
        infoTable.classList.add('metabolite-info-grid');
        const molWeight = data.metabolite_mol_weights[index];
        const infoFields = [
            ['Charged Formula', data.metabolite_formulas[index]],
            ['SMILES', data.metabolite_smiles[index]],
            ['InChI', data.metabolite_inchis[index]],
            ['InChI Key', data.metabolite_inchi_keys[index]],
            // Display includes the unit; copy the bare numeric value.
            ['Molecular Weight', (molWeight != null && molWeight !== '') ? `${molWeight} g/mol` : null, molWeight],
        ];
        infoFields.forEach(function (field) {
            infoTable.appendChild(createMetaboliteInfoItem(field[0], field[1], field[2]));
        });
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

/**
 * Build one metabolite info row: a label, a single-line value that truncates
 * with an ellipsis (full text on hover), and a copy button that copies the full
 * untruncated value. `copyValue` overrides what is copied when the displayed
 * text differs from the raw value (e.g. molecular weight shows a unit).
 */
function createMetaboliteInfoItem(label, displayValue, copyValue) {
    const hasValue = displayValue != null && displayValue !== '';
    const shown = hasValue ? String(displayValue) : 'N/A';

    const item = document.createElement('div');
    item.className = 'info-item';

    const labelEl = document.createElement('strong');
    labelEl.className = 'info-label';
    labelEl.textContent = label + ':';
    item.appendChild(labelEl);

    const valueEl = document.createElement('span');
    valueEl.className = 'info-value';
    valueEl.textContent = shown;
    valueEl.title = shown; // hover reveals the full value
    item.appendChild(valueEl);

    if (hasValue) {
        const toCopy = String(copyValue != null ? copyValue : displayValue);
        const copyBtn = document.createElement('button');
        copyBtn.type = 'button';
        copyBtn.className = 'info-copy';
        copyBtn.title = 'Copy ' + label;
        copyBtn.setAttribute('aria-label', 'Copy ' + label);
        copyBtn.innerHTML = '<i class="fas fa-copy" aria-hidden="true"></i>';
        copyBtn.addEventListener('click', function () {
            copyMetaboliteValue(toCopy, copyBtn);
        });
        item.appendChild(copyBtn);
    }
    return item;
}

/** Copy `text` to the clipboard and give brief in-place feedback on the button. */
function copyMetaboliteValue(text, btn) {
    const feedback = function (ok) {
        const icon = btn ? btn.querySelector('i') : null;
        if (icon) {
            const prev = icon.className;
            icon.className = ok ? 'fas fa-check' : 'fas fa-times';
            btn.classList.add(ok ? 'copied' : 'copy-failed');
            setTimeout(function () {
                icon.className = prev;
                btn.classList.remove('copied', 'copy-failed');
            }, 1200);
        }
        if (typeof Notify !== 'undefined') {
            if (ok) Notify.success('Copied to clipboard');
            else Notify.error('Could not copy');
        }
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text)
            .then(function () { feedback(true); })
            .catch(function () { fallbackCopyText(text, feedback); });
    } else {
        fallbackCopyText(text, feedback);
    }
}

function fallbackCopyText(text, feedback) {
    try {
        const ta = document.createElement('textarea');
        ta.value = text;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.focus();
        ta.select();
        const ok = document.execCommand('copy');
        document.body.removeChild(ta);
        if (feedback) feedback(ok);
    } catch (e) {
        if (feedback) feedback(false);
    }
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
