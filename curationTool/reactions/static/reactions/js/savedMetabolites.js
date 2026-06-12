document.addEventListener('DOMContentLoaded', function () {
    const searchInput = document.getElementById('metaboliteSearchInput');
    if (searchInput) {
        searchInput.addEventListener('input', function () {
            const searchTerm = this.value.toLowerCase();
            document.querySelectorAll('#metabolitesList .item').forEach(item => {
                const nameElem = item.querySelector('.name');
                const abbrElem = item.querySelector('.abbr');
                const nameText = nameElem ? nameElem.textContent.toLowerCase() : '';
                const abbrText = abbrElem ? abbrElem.textContent.toLowerCase() : '';
                if (nameText.includes(searchTerm) || abbrText.includes(searchTerm)) {
                    item.style.display = '';  // Show item
                } else {
                    item.style.display = 'none';  // Hide item
                }
            });
        });
    }

    // Share selected metabolites. This script is reused across pages/modal variants,
    // so optional controls need guards.
    const shareSelectedBtn = document.getElementById('shareSelectedBtn');
    if (shareSelectedBtn) {
        shareSelectedBtn.addEventListener('click', shareSelectedMetabolites);
    }

});
async function shareSelectedMetabolites() {
    const selectedCheckboxes = document.querySelectorAll('.share-checkbox:checked');
    if (selectedCheckboxes.length === 0) {
        Notify.error('Please select at least one metabolite to share.');
        return;
    }

    // Show the warning/confirmation dialog.
    const proceed = await Notify.confirm({
        title: 'Share metabolites',
        message: 'Sharing these metabolites will allow the recipient to make changes that affect your copy. Continue?',
        confirmText: 'Share',
        cancelText: 'Cancel',
        danger: true,
    });
    if (!proceed) {
        return;
    }

    // Prompt for the target username.
    const targetUsername = await Notify.prompt({
        title: 'Share with',
        message: 'Enter the username to share with:',
        placeholder: 'Username',
        confirmText: 'Share',
    });
    if (!targetUsername) {
        Notify.error('Username is required.');
        return;
    }
    const sharingUser = sessionStorage.getItem('userID');
    if (!sharingUser) {
        Notify.error('Please log in to share metabolites.');
        return;
    }
    // Gather selected metabolite IDs.
    const metaboliteIds = Array.from(selectedCheckboxes).map(cb => cb.getAttribute('data-metabolite-id'));
    
    // Make the API call.
    fetch('/share_metabolites/', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrfToken  // assuming csrfToken is defined globally
        },
        body: JSON.stringify({
            metabolite_ids: metaboliteIds,
            target_username: targetUsername,
            sharing_user_id: sharingUser
        })
    })
    .then(response => response.json())
    .then(data => {
        if (data.status === 'success') {
            showToast(data.message);
            // Optionally clear checkboxes.
            selectedCheckboxes.forEach(cb => cb.checked = false);
        } else {
            showToast('Error sharing metabolites: ' + data.message, '#cc0000');
        }
    })
    .catch(error => {
        showToast('Error sharing metabolites: ' + error, '#cc0000');
    });
}


function toggleEditField(button, fieldType) {
    const input = button.parentNode.querySelector('input');
    const isEditing = input.readOnly;
    const generateButton = input.parentNode.querySelector('.generate-abbr-btn');

    input.readOnly = !isEditing;
    input.classList.toggle('readonly-field');
    button.textContent = isEditing ? 'Save' : 'Edit';

    if (isEditing) {
        // Store the original value when entering edit mode
        input.dataset.originalValue = input.value;
        input.focus();

        if (fieldType === 'abbr' && generateButton) {
            generateButton.style.display = 'inline-block';
        }
    } else {
        if (fieldType === 'abbr' && generateButton) {
            generateButton.style.display = 'none';
        }

        // Save changes
        const metaboliteId = input.dataset.metaboliteId;
        const originalValue = input.dataset.originalValue;
        const value = input.value.trim();
        console.log('Saving changes:', { metaboliteId, fieldType, originalValue, value });
        // Determine the correct data structure
        let data;
        if (fieldType === 'name') {
            data = { name: value };
        } else if (fieldType === 'abbr') {
            data = { vmh_abbr: value };
        } else {
            data = { [fieldType]: value };
        }

        // API call to update the metabolite
        fetch(`/update_metabolite/${metaboliteId}/`, {
            method: 'PUT',
            headers: {
                'Content-Type': 'application/json',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify(data)
        })
        .then(response => response.json())
        .then(data => {
            if (data.status === 'error') {
                showToast('Error saving changes: ' + data.message, '#cc0000');
                input.value = originalValue;  // Revert to the original value if there's an error
            } else {
                const fieldLabel = (fieldType === 'name' ? 'name' : (fieldType === 'abbr' ? 'abbreviation' : fieldType) ?? fieldType);
                showToast(`Changes saved successfully for ${fieldLabel}`);
                // Instead of refreshing the modal, update only the changed row:
                updateMetaboliteRow(data.metabolite);
            }
        })        
        .catch(error => {
            showToast('Error saving changes: ' + error, '#cc0000');
            input.value = originalValue;  // Revert on network error
        });
    }
}
function updateMetaboliteRow(updatedMetabolite) {
    // Locate the metabolite row. For example, assume you add a data attribute to the main container:
    const row = document.querySelector(`.content[data-metabolite-id="${updatedMetabolite.id}"]`);
    if (!row) return;

    // Update header fields
    const abbrHeader = row.querySelector('.savedMetabolite-header .abbr');
    if (abbrHeader) {
        abbrHeader.innerHTML = updatedMetabolite.vmh_abbr || '<span style="color: #999;">no abbr</span>';
    }
    const nameHeader = row.querySelector('.savedMetabolite-header .name');
    if (nameHeader) {
        nameHeader.textContent = updatedMetabolite.name;
    }
    
    // Update input fields in the details section, if they exist
    const abbrInput = row.querySelector('.abbr-input');
    if (abbrInput) {
        abbrInput.value = updatedMetabolite.vmh_abbr;
    }
    const nameInput = row.querySelector('.name-input');
    if (nameInput) {
        nameInput.value = updatedMetabolite.name;
    }
    
}

async function loadSavedMetabolites() {
    const listContainer = document.getElementById('metabolitesList');
    listContainer.innerHTML = `
        <div class="ui active inverted dimmer">
            <div class="ui text loader">Loading saved metabolites...</div>
        </div>`;
    
    try {
        const userId = sessionStorage.getItem('userID');
        const response = await fetch(`/get_saved_metabolites/?user_id=${userId}`);
        const { metabolites } = await response.json();
        listContainer.innerHTML = '';  // Clear the loader
        metabolites.forEach(met => {
            console.log('Metabolite:', met);
            const metaboliteDiv = document.createElement('div');
            metaboliteDiv.className = 'item';
            const externalLinksHtml = Object.entries(met.external_links || {}).map(([key, value]) => `
            <div class="field">
                <label>${key}</label>
                <div class="ui action input">
                    <input type="text" value="${value || ''}" class="readonly-field" readonly data-id="${key}" data-metabolite-id="${met.id}">
                    <button class="ui button edit-button" onclick="toggleEditField(this, '${key}')">Edit</button>
                </div>
            </div>
            `).join('');
            const reactionsHtml = `
            <div class="reactions-toggle-container">
                <button class="ui button toggle-reactions" onclick="toggleReactions(this)">Show Reactions</button>
                <div class="reactions-container" style="display: none;">
                    ${met.reactions.length > 0 
                        ? `<strong>In reactions:</strong>
                        <ul>
                            ${met.reactions.map(r => `<li><a href="/?reaction_id=${r.id}">${r.name}</a></li>`).join('')}
                        </ul>`
                        : '<strong>Not used in any reactions</strong>'}
                </div>
            </div>
            `;
            metaboliteDiv.innerHTML = `
            <div class="content" data-metabolite-id="${met.id}">
                <div class="header savedMetabolite-header">
                    <div class="left-section">
                        <div class="share-selection">
                            <input type="checkbox" class="share-checkbox" data-metabolite-id="${met.id}">
                        </div>
                        <span class="toggle-arrow">▶</span>
                        <span class="abbr">${met.vmh_abbr || '<span style="color: #999;">no abbr</span>'}</span>
                        <span class="divider"> | </span>
                        <span class="name">${met.name}</span>
                    </div>
                    <div class="buttons">
                        <button class="ui primary button download-btn" 
                                data-mol-file="${btoa(met.mol_file)}" 
                                data-filename="${met.name.replace(/ /g, '_')}.mol" 
                                onclick="downloadMolFile(this)">
                            Download
                        </button>
                        <button class="ui negative button delete-btn" data-metabolite-id="${met.id}" onclick="deleteMetabolite(${met.id})">
                            Delete
                        </button>
                    </div>
                </div>
        
                <div class="description metabolite-details" style="display: none;">
                    <div class="ui grid compact-info-grid">
                        <!-- Left Column (Editable Info) -->
                        <div class="eight wide column">
                            <div class="ui form small-form">
                                <div class="twofields">
                                    <div class="field">
                                        <label>Abbreviation</label>
                                        <div class="ui action input">
                                            <input type="text" value="${met.vmh_abbr || ''}" 
                                                class="abbr-input readonly-field" readonly data-metabolite-id="${met.id}">
                                            <button class="ui button generate-abbr-btn" 
                                                    onclick="generateAbbreviation(this)" style="display: none;">
                                                Generate
                                            </button>
                                            <button class="ui button edit-button" 
                                                    onclick="toggleEditField(this, 'abbr')">
                                                Edit
                                            </button>
                                        </div>
                                    </div>
                                    <div class="field">
                                        <label>Name</label>
                                        <div class="ui action input">
                                            <input type="text" value="${met.name}" 
                                                class="name-input readonly-field" readonly data-metabolite-id="${met.id}">
                                            <button class="ui button edit-button" 
                                                    onclick="toggleEditField(this, 'name')">
                                                Edit
                                            </button>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            ${reactionsHtml}
                            <div class="external-links-container">
                                <button class="ui button toggle-links" onclick="toggleExternalLinks(this)">Show External Links</button>
                                <div class="links-container" style="display: none;">
                                    ${externalLinksHtml}
                                </div>
                            </div>
                        </div>
        
                        <!-- Right Column (Molecular Info) -->
                        <div class="eight wide column">
                            <div class="ui list compact-info-list">
                                <div class="item">
                                    <strong>InChI Key:</strong> <span>${met.inchi_key || 'N/A'}</span>
                                </div>
                                <div class="item">
                                    <strong>InChI:</strong> <span>${met.inchi || 'N/A'}</span>
                                </div>
                                <div class="item">
                                    <strong>SMILES:</strong> <span>${met.smiles || 'N/A'}</span>
                                </div>
                                <div class="item">
                                    <strong>Molecular Weight:</strong> <span>${met.mol_w || 'N/A'}</span>
                                </div>
                                <div class="item">
                                    <strong>Formula:</strong> <span>${met.mol_formula || 'N/A'}</span>
                                </div>
                            </div>
        
                            <div class="structure-controls">
                                <button class="ui button toggle-3d" 
                                        onclick="toggle3DView(this, '${btoa(met.mol_file)}')">
                                    Show 3D Structure
                                </button>
                            </div>
                            <div class="structure-container" style="display: none; height: 400px; width: 100%; position: relative;">
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;    
            listContainer.appendChild(metaboliteDiv);
            document.querySelectorAll('.share-checkbox').forEach(checkbox => {
                checkbox.addEventListener('click', function (e) {
                    e.stopPropagation(); // This stops the click from bubbling up
                });
            });
        });

        // Toggle animation and details
        document.querySelectorAll('.savedMetabolite-header').forEach(header => {
            header.addEventListener('click', function(e) {
                if(!e.target.closest('.delete-btn, .edit-button')) {
                    this.classList.toggle('active');
                    const details = this.parentNode.querySelector('.metabolite-details');
                    details.style.display = details.style.display === 'none' ? 'block' : 'none';
                }
            });
        });
    }
    catch (error) {
        listContainer.innerHTML = '<p>Error loading metabolites.</p>';
        console.error(error);
    }
}

function toggleExternalLinks(button) {
    const container = button.nextElementSibling;
    container.style.display = container.style.display === 'none' ? 'flex' : 'none';
    button.textContent = container.style.display === 'none' ? 'Show External Links' : 'Hide External Links';
}
function toggleReactions(button) {
    const container = button.nextElementSibling;
    container.style.display = container.style.display === 'none' ? 'block' : 'none';
    button.textContent = container.style.display === 'none' ? 'Show Reactions' : 'Hide Reactions';
}
function downloadMolFile(buttonElement) {
    const molFileBase64 = buttonElement.getAttribute('data-mol-file');
    const filename = buttonElement.getAttribute('data-filename');
    const molFileContent = atob(molFileBase64);
    
    const blob = new Blob([molFileContent], { type: 'chemical/x-mdl-molfile' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
}

// 3D Structure Viewer
function toggle3DView(button, molFileBase64) {
    const container = button.parentNode.nextElementSibling; // Target the structure-container
    if (container.style.display === 'none') {
        container.style.display = 'block';
        button.textContent = 'Hide 3D Structure';
        
        if (!container.hasAttribute('data-initialized')) {
            const molFile = atob(molFileBase64); // Decode base64 string
            const viewer = $3Dmol.createViewer(container, {
                backgroundColor: 'white'
            });
            viewer.addModel(molFile, "mol");
            viewer.setStyle({}, { 
                stick: { radius: 0.15 }, 
                sphere: { scale: 0.25 } 
            });
            viewer.zoomTo();
            viewer.render();
            container.setAttribute('data-initialized', 'true');
        }
    } else {
        container.style.display = 'none';
        button.textContent = 'Show 3D Structure';
    }
}
// Handle metabolite deletion
async function deleteMetabolite(metaboliteId) {
    try {
        const response = await fetch(`/delete_metabolite/${metaboliteId}/`, {
            method: 'DELETE',
            headers: {
                'X-CSRFToken': csrfToken,
            },
            body: JSON.stringify({ user_id: sessionStorage.getItem('userID') })
        });
        const data = await response.json();

        if (data.status === 'needs_confirmation') {
            // Build confirmation message with reactions list
            let message = data.message + '\n\nReactions to be deleted:\n';
            data.reactions.forEach(r => {
                message += `• ${r.short_name || 'Reaction ' + r.id}\n`;
            });
            const isConfirmed = await Notify.confirm({ title: 'Delete metabolite', message: message, confirmText: 'Delete', cancelText: 'Cancel', danger: true });
            if (isConfirmed) {
                const confirmResponse = await fetch(`/delete_metabolite/${metaboliteId}/?confirm=true`, {
                    method: 'DELETE',
                    headers: {
                        'X-CSRFToken': csrfToken
                    },
                    body: JSON.stringify({ user_id: sessionStorage.getItem('userID') })
                });
                const confirmData = await confirmResponse.json();
                if (confirmData.status === 'success') {
                    loadSavedMetabolites();
                    showToast('Metabolite and related reactions deleted successfully.');
                } else {
                    showToast('Deletion failed: ' + (confirmData.message || 'Unknown error'), '#cc0000');
                }
            }
        } else if (data.status === 'success') {
            loadSavedMetabolites();
            showToast(data.message || 'Metabolite deleted successfully.');
        } else {
            showToast('Error: ' + (data.message || 'Failed to delete metabolite.'), '#cc0000');
        }
    } catch (error) {
        showToast('Error: ' + error.message, '#cc0000');
    }
}

// New metabolite creation
const createNewMetaboliteBtn = document.getElementById('createNewMetaboliteBtn');
if (createNewMetaboliteBtn) {
    createNewMetaboliteBtn.addEventListener('click', function() {
        const form = document.getElementById('newMetaboliteForm');
        if (form) form.style.display = 'block';
        this.style.display = 'none';
    });
}

async function createNewMetabolite() {
    const nameInput = document.getElementById('newMetaboliteName');
    const abbrInput = document.getElementById('newMetaboliteAbbr');
    const molFileInput = document.getElementById('newMetaboliteMolFile');
    if (!nameInput || !abbrInput || !molFileInput) {
        Notify.error('The new metabolite form is not available on this page.');
        return;
    }
    const formData = new FormData();
    formData.append('name', nameInput.value);
    formData.append('vmh_abbr', abbrInput.value);
    formData.append('mol_file', molFileInput.files[0]);
    formData.append('user_id', sessionStorage.getItem('userID'));

    const response = await fetch('/create_metabolite/', {
        method: 'POST',
        headers: {
            'X-CSRFToken': csrfToken
        },
        body: formData
    });

    if(response.ok) {
        loadSavedMetabolites();
        cancelNewMetabolite();
        showToast('Metabolite created successfully');
    } else {
        showToast('Error creating metabolite', 'error');
    }
}

function cancelNewMetabolite() {
    const form = document.getElementById('newMetaboliteForm');
    const createBtn = document.getElementById('createNewMetaboliteBtn');
    const nameInput = document.getElementById('newMetaboliteName');
    const abbrInput = document.getElementById('newMetaboliteAbbr');
    const molFileInput = document.getElementById('newMetaboliteMolFile');
    if (form) form.style.display = 'none';
    if (createBtn) createBtn.style.display = 'block';
    if (nameInput) nameInput.value = '';
    if (abbrInput) abbrInput.value = '';
    if (molFileInput) molFileInput.value = '';
}

async function generateAbbreviation(button) {
    try {
        // Get the closest metabolite header
        const metaboliteHeader = button.closest('.content').querySelector('.savedMetabolite-header');
        if (!metaboliteHeader) {
            showToast('Could not find metabolite header', '#cc0000');
            return;
        }

        // Get the metabolite name from the header
        const nameElement = metaboliteHeader.querySelector('.name');
        if (!nameElement) {
            showToast('Metabolite name element not found', '#cc0000');
            return;
        }
        const metaboliteName = nameElement.textContent.trim();

        // Get the abbreviation input field
        const input = button.parentNode.querySelector('input');
        if (!input) {
            showToast('Abbreviation input field not found', '#cc0000');
            return;
        }

        // Disable button to prevent multiple clicks
        button.disabled = true;
        button.textContent = 'Generating...';
        saveBtn = button.parentNode.querySelector('.edit-button');
        saveBtn.disabled = true;
        const response = await fetch('/create-formula-abbr/', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({
                metabolite: null,
                mtype: 'saved',
                metabolite_name: metaboliteName
            })
        });

        const data = await response.json();
        if (data.abbr) {
            input.value = data.abbr;
            input.dispatchEvent(new Event('change'));
            showToast('Abbreviation generated successfully');
        } else {
            showToast('Failed to generate abbreviation', '#cc0000');
        }
    } catch (error) {
        console.error('Error:', error);
        showToast('Error generating abbreviation', '#cc0000');
    } finally {
        // Re-enable button and clean up
        button.disabled = false;
        saveBtn.disabled = false;
        button.textContent = 'Generate';
    }
}
