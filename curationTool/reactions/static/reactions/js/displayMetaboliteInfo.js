function loadMetaboliteInfoDiv(reactionData) {
    const metaboliteInfoDiv = document.querySelector('.content-div[name="metaboliteinfo-div"]');
    if (metaboliteInfoDiv) {
        metaboliteInfoDiv.innerHTML = `
        <div class="div-header">Metabolite Information</div>
        `;
        fillMetaboliteInfoTab(reactionData);
    } else {
        console.error('Metabolite info div not found.');
    }
}


function fillMetaboliteInfoTab(data) {
    const metaboliteInfoContainer = document.getElementById('metaboliteinfo-div');

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

                    let config = { backgroundColor: 'white' };
                    var viewer = new $3Dmol.createViewer(structureContainer, config);
                    var molecularData = data.metabolite_mol_file_strings[index];
                    var format = 'sdf';
                    viewer.addModel(molecularData, format);
                    viewer.setClickable({}, true, function(atom, _viewer, _event, _container) {
                        viewer.addLabel(atom.atom, { position: atom, backgroundColor: 'darkgreen', backgroundOpacity: 0.8 });
                    });
                    viewer.setStyle({}, { 'stick': { 'colorscheme': 'greenCarbon' } });
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

        const formulaElement = document.createElement('p');
        formulaElement.textContent = `Charged Formula: ${data.metabolite_formulas[index]}`;
        metaboliteDiv.appendChild(formulaElement);

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
