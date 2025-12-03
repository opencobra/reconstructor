let eventListenersAttached = false; // Flag to track if event listeners have been attached

function createSectionHTML(sectionTitle, className, items, reactionId, isExtLink = false, isRef = false, isGene = false) {
	if (isGene) {
		return createGeneSectionHTML(sectionTitle, className, items, reactionId);
	}

	let sectionHTML = `
        <div class="${className}-section">
            <p>${sectionTitle}</p>
            ${items
							.map(
								(item, index) => `
                <div class="${className}-item" data-reaction-id="${reactionId}" data-index="${index}"> <!-- Ensure this div represents an ext-link-item for external links -->
                    ${isExtLink ? createExtLinkSelect(item, reactionId, index) : ''}
                    ${isRef ? createRefSelect(item, reactionId, index) : ''}
                    <input type="text" class="${className}-input" placeholder="Enter ${sectionTitle.toLowerCase()}" value="${
									item.info
								}" data-reaction-id="${reactionId}" data-index="${index}">
                    <button class="remove-${className}" data-reaction-id="${reactionId}" data-index="${index}">Remove</button>
                </div>
            `
							)
							.join('')}
            <button class="add-${className}" data-reaction-id="${reactionId}">Add ${sectionTitle.slice(0, -1)}</button>
        </div>
    `;
	return sectionHTML;
}

function createExtLinkSelect(item, reactionId, index) {
	const linkTypes = ['CHO Models', 'COG', 'EC Number', 'KEGG orthology', 'KEGG reaction', 'MetanetX', 'Rhea', 'SEED', 'Wikipedia'];
	let selectHTML = `
        <select class="ext-link-type-select" data-reaction-id="${reactionId}" data-index="${index}">
            ${linkTypes.map((type) => `<option value="${type}" ${item.ext_link_type === type ? 'selected' : ''}>${type}</option>`).join('')}
        </select>
    `;
	return selectHTML;
}
function createRefSelect(item, reactionId, index) {
	const refTypes = ['DOI', 'PMID'];
	let selectHTML = `
        <select class="ref-type-select" data-reaction-id="${reactionId}" data-index="${index}">
            ${refTypes.map((type) => `<option value="${type}" ${item.ref_type === type ? 'selected' : ''}>${type}</option>`).join('')}
        </select>
    `;
	return selectHTML;
}

function createGeneSectionHTML(sectionTitle, className, items, reactionId) {
	console.log(items);
	console.log(reactionId);
	let sectionHTML = `
        <div class="${className}-section">
            <p>${sectionTitle}</p>
            ${items
							.map(
								(item, index) => `
                <div class="${className}-item" data-reaction-id="${reactionId}" data-index="${index}">
                <input type="text" class="${className}-input" placeholder="Enter ${sectionTitle.toLowerCase()}" value="${
									item.info
								}" data-reaction-id="${reactionId}" data-index="${index}" readonly>
                <button class="remove-${className}" data-reaction-id="${reactionId}" data-index="${index}">Remove</button>
                </div>
            `
							)
							.join('')}
        </div>
    `;
	return sectionHTML;
}

async function callPrepareAddToVMH(reactionIds) {
	try {
		const response = await fetch(prepAddToVMHUrl, {
			method: 'POST',
			headers: {
				'X-Requested-With': 'XMLHttpRequest',
				'Content-Type': 'application/json',
				'X-CSRFToken': csrfToken,
			},
			body: JSON.stringify({
				reactionIds: reactionIds,
			}),
		});

		if (!response.ok) {
			// Handle HTTP errors
			console.error(`HTTP error! status: ${response.status}`);
			return null;
		}

		const data = await response.json();
		return data;
	} catch (error) {
		// Handle network or other errors
		console.error('Error fetching data: ', error);
		return null; // Or handle error appropriately
	}
}

function toggleDetails(event, headerElement) {
	// Check if the click target is one of the input fields or the select dropdown
	if (
		event.target.classList.contains('reaction-name-input') ||
		event.target.classList.contains('reaction-abbreviation-input') ||
		event.target.tagName.toLowerCase() === 'select' ||
		event.target.classList.contains('cs-info-button')
	) {
		return; // Do nothing if clicked inside the input fields or select dropdown
	}

	// Find the toggle icon within the header element
	const toggleIcon = headerElement.querySelector('.toggle-icon');

	// Assume the detailed sections are initially visible; we need to check and set their initial state if undefined
	let initialStateSet = false; // Flag to track if we've set the initial state

	// The detailed sections are immediate siblings of the header, so we toggle their display.
	let currentElement = headerElement.nextElementSibling;
	while (currentElement) {
		// If initial state isn't set, explicitly set it to "none" to ensure consistent behavior
		if (!initialStateSet && (currentElement.style.display === '' || currentElement.style.display === 'block')) {
			currentElement.style.display = 'none';
			initialStateSet = true; // Mark that we've set the initial state
		} else {
			// Toggle the display style
			if (currentElement.style.display === 'none') {
				currentElement.style.display = 'block'; // Show the section
			} else {
				currentElement.style.display = 'none'; // Hide the section
			}
		}
		// Move to the next sibling
		currentElement = currentElement.nextElementSibling;
	}

	// After setting the initial state or toggling, update the icon accordingly
	if (initialStateSet) {
		// If we just set the initial state, ensure the icon is set for the next action
		toggleIcon.textContent = '+'; // Set to plus since we've just hidden the sections
	} else {
		// Toggle the icon between "+" and "−"
		toggleIcon.textContent = toggleIcon.textContent === '+' ? '−' : '+';
	}
}

function toggleInfo() {
	// Create modal div
	var modal = document.createElement('div');
	modal.id = 'myModal';
	modal.className = 'modal';

	// Create modal content div
	var modalContent = document.createElement('div');
	modalContent.className = 'modal-content';

	// Create close button
	// var closeButton = document.createElement('span');
	// closeButton.innerHTML = '&times;';
	// closeButton.className = 'close';
	// closeButton.onclick = function() {
	//     modal.style.display = 'none';
	// };

	// // Append close button to modal content
	// modalContent.appendChild(closeButton);

	// Create modal header
	var modalHeader = document.createElement('h2');
	modalHeader.innerText = 'Confidence score information';
	modalContent.appendChild(modalHeader);

	// Create modal body text
	var modalBody = document.createElement('p');
	modalBody.innerText =
		'At the core of the human metabolism resource lies the human metabolic reconstruction, which is amenable to metabolic modeling. Consequently, metabolic functions without genetic evidence but with physiological evidence have also been included in this resource. Many of the reactions are associated with a confidence score indicating the evidence supporting their inclusion into the human metabolic reconstruction.';
	modalContent.appendChild(modalBody);

	// Create table
	var table = document.createElement('table');
	var tableHtml = `
    <table>
    <tr>
        <th>Evidence type</th>
        <th style="text-align: center;">Confidence score</th>
        <th>Examples</th>
    </tr>
    <tr>
        <td>Biochemical data</td>
        <td style="text-align: center;">4</td>
        <td>Direct evidence for gene product function and biochemical reaction: Protein purification, biochemical assays, experimentally solved protein structures, and comparative gene-expression studies.</td>
    </tr>
    <tr>
        <td>Genetic data</td>
        <td style="text-align: center;">3</td>
        <td>Direct and indirect evidence for gene function: Knock-out characterization, knock-in characterization, and over-expression.</td>
    </tr>
    <tr>
        <td>Physiological data</td>
        <td style="text-align: center;">2</td>
        <td>Indirect evidence for biochemical reactions based on physiological data: secretion products or defined medium components serve as evidence for transport and metabolic reactions.</td>
    </tr>
    <tr>
        <td>Sequence data</td>
        <td style="text-align: center;">2</td>
        <td>Evidence for gene function: Genome annotation.</td>
    </tr>
    <tr>
        <td>Modeling data</td>
        <td style="text-align: center;">1</td>
        <td>No evidence is available but reaction is required for modeling. The included function is a hypothesis and needs experimental verification. The reaction mechanism may be different from the included reaction(s).</td>
    </tr>
    <tr>
        <td>Not evaluated</td>
        <td style="text-align: center;">0</td>
        <td></td>
    </tr>
</table>

    `;
	table.innerHTML = tableHtml;
	modalContent.appendChild(table);

	// Create citation
	var citation = document.createElement('p');
	citation.innerText =
		'Taken from Thiele, I., Palsson, B. O., "A protocol for generating a high-quality genome-scale metabolic reconstruction.", Nat Protocols, 5(1): 93 - 121 (2010).';
	modalContent.appendChild(citation);

	// Create OK button
	var okButton = document.createElement('button');
	okButton.innerText = 'Ok';
	okButton.onclick = function () {
		modal.style.display = 'none';
	};
	modalContent.appendChild(okButton);

	// Append modal content to modal
	modal.appendChild(modalContent);

	// Append modal to body
	document.body.appendChild(modal);

	// Display the modal
	modal.style.display = 'block';
}

// VMH Submission Modal Functions
function showVMHModal() {
	const modal = document.getElementById('vmhSubmissionModal');
	modal.classList.add('active');
}

function hideVMHModal() {
	const modal = document.getElementById('vmhSubmissionModal');
	modal.classList.remove('active');
}

function updateVMHModalStatus(state, data = {}) {
	const statusContainer = document.getElementById('vmhSubmissionStatus');
	const footer = document.getElementById('vmhSubmissionFooter');
	const title = document.getElementById('vmhSubmissionTitle');
	
	switch(state) {
		case 'loading':
			title.textContent = 'Adding to VMH';
			footer.style.display = 'none';
			statusContainer.innerHTML = `
				<div class="vmh-loading-spinner"></div>
				<div class="vmh-status-text">${data.text || 'Processing...'}</div>
				<div class="vmh-status-subtext">${data.subtext || 'Please wait while we process your request'}</div>
				<div class="vmh-progress-steps">
					<div class="vmh-step ${data.step >= 1 ? (data.step > 1 ? 'vmh-step-complete' : 'vmh-step-active') : 'vmh-step-pending'}">
						<div class="vmh-step-icon"><i class="fas ${data.step > 1 ? 'fa-check' : 'fa-clipboard-check'}"></i></div>
						<div class="vmh-step-label">Validating requirements</div>
					</div>
					<div class="vmh-step ${data.step >= 2 ? (data.step > 2 ? 'vmh-step-complete' : 'vmh-step-active') : 'vmh-step-pending'}">
						<div class="vmh-step-icon"><i class="fas ${data.step > 2 ? 'fa-check' : 'fa-cogs'}"></i></div>
						<div class="vmh-step-label">Initializing Cobra ToolBox</div>
					</div>
					<div class="vmh-step ${data.step >= 3 ? 'vmh-step-active' : 'vmh-step-pending'}">
						<div class="vmh-step-icon"><i class="fas fa-database"></i></div>
						<div class="vmh-step-label">Adding to VMH database</div>
					</div>
				</div>
			`;
			break;
			
		case 'success':
			title.textContent = 'Success!';
			footer.style.display = 'flex';
			
			let reactionsHTML = '';
			let metabolitesHTML = '';
			
			if (data.reactions && data.reactions.length > 0) {
				reactionsHTML = `
					<div class="vmh-result-section">
						<h4><i class="fas fa-flask"></i> Reactions Added (${data.reactions.length})</h4>
						${data.reactions.map(([abbr, [id, formula]]) => `
							<div class="vmh-result-item" onclick="this.classList.toggle('expanded')">
								<div style="display: flex; align-items: center; gap: 8px; flex: 1;">
									<span class="vmh-result-abbr">${abbr}</span>
									<span class="vmh-result-id">(ID: ${id})</span>
								</div>
								<i class="fas fa-chevron-down" style="font-size: 10px; color: #999;"></i>
								<div class="vmh-result-formula">${formula}</div>
							</div>
						`).join('')}
					</div>
				`;
			}
			
			if (data.metabolites && data.metabolites.length > 0) {
				metabolitesHTML = `
					<div class="vmh-result-section">
						<h4><i class="fas fa-atom"></i> Metabolites Added (${data.metabolites.length})</h4>
						${data.metabolites.map(([abbr, [id, formula, inchiKey]]) => `
							<div class="vmh-result-item" onclick="this.classList.toggle('expanded')">
								<div style="display: flex; align-items: center; gap: 8px; flex: 1;">
									<span class="vmh-result-abbr">${abbr}</span>
									<span class="vmh-result-id">(ID: ${id})</span>
								</div>
								<i class="fas fa-chevron-down" style="font-size: 10px; color: #999;"></i>
								<div class="vmh-result-formula">Formula: ${formula}<br>InChIKey: ${inchiKey}</div>
							</div>
						`).join('')}
					</div>
				`;
			}
			
			statusContainer.innerHTML = `
				<div class="vmh-success-icon"><i class="fas fa-check"></i></div>
				<div class="vmh-success-title">Successfully Added to VMH!</div>
				<div class="vmh-success-subtitle">${data.reactions?.length || 0} reaction(s) and ${data.metabolites?.length || 0} metabolite(s) were added</div>
				<div class="vmh-result-details">
					${reactionsHTML}
					${metabolitesHTML}
				</div>
			`;
			break;
			
		case 'error':
			title.textContent = 'Unable to Add';
			footer.style.display = 'flex';
			statusContainer.innerHTML = `
				<div class="vmh-error-icon"><i class="fas fa-exclamation-circle"></i></div>
				<div class="vmh-error-message-main">${data.message || 'An unexpected error occurred. Please try again.'}</div>
				<div class="vmh-error-hint">Please fix the issue above and try again.</div>
			`;
			break;
	}
}

function addToVMH() {
	setButtonState(true);

	if (!validateInputs()) {
		displayValidationMessage(true, 'Fill in all non-VMH metabolite names');
		setButtonState(false);
		
		// Show inline error message in workspace details
		const detailPanel = document.getElementById('wsAvailableReactionDetails');
		if (detailPanel) {
			let errorMsg = detailPanel.querySelector('.ws-validation-message');
			if (!errorMsg) {
				errorMsg = document.createElement('div');
				errorMsg.className = 'ws-validation-message';
				errorMsg.innerHTML = '<i class="fas fa-exclamation-circle"></i> Please fill in all required fields highlighted in red';
				detailPanel.insertBefore(errorMsg, detailPanel.firstChild);
			}
			// Scroll to top of detail panel
			detailPanel.scrollTop = 0;
		}
		
		return;
	} else {
		displayValidationMessage(false);
		// Remove inline error message if exists
		const errorMsg = document.querySelector('.ws-validation-message');
		if (errorMsg) errorMsg.remove();
	}

	// Show the professional modal
	showVMHModal();
	updateVMHModalStatus('loading', { 
		text: 'Validating requirements...', 
		subtext: 'Checking if all fields are properly filled',
		step: 1 
	});

	// Progress through steps
	let currentStep = 1;
	const stepTexts = [
		{ text: 'Validating requirements...', subtext: 'Checking if all fields are properly filled', step: 1 },
		{ text: 'Initializing Cobra ToolBox...', subtext: 'Preparing the metabolic modeling tools', step: 2 },
		{ text: 'Adding to VMH database...', subtext: 'This may take a moment', step: 3 }
	];
	
	const stepInterval = setInterval(() => {
		currentStep++;
		if (currentStep <= 3) {
			updateVMHModalStatus('loading', stepTexts[currentStep - 1]);
		}
	}, 7500);
	
	const updatedReactions = checkedReactions.map((reactionId) => {
		const reactionIdNum = Number(reactionId);
		const reaction = reactions_active.find((r) => r.pk === reactionIdNum);
		
		if (!reaction) {
			console.error(`Reaction with pk ${reactionId} not found`);
			return null;
		}

		// Check if form inputs exist for this reaction (i.e., it's currently displayed)
		const nameInputField = document.querySelector(`.reaction-name-input[data-reaction-id="${reactionId}"]`);
		const hasFormInputs = nameInputField !== null;

		let reactionAbbr = '';
		let subsDetails = [];
		let prodsDetails = [];
		let referencesData = [];
		let extLinksData = [];
		let commentsData = [];
		let selectedValue = reaction.fields.confidence_score || ' ';

		if (hasFormInputs) {
			// Gather from form inputs (reaction is currently displayed)
			const abbrInputField = document.querySelector(`.reaction-abbreviation-input[data-reaction-id="${reactionId}"]`);
			if (nameInputField && nameInputField.value.trim() !== '') {
				reaction.fields.description = nameInputField.value.trim();
			}
			if (abbrInputField && abbrInputField.value.trim() !== '') {
				reactionAbbr = abbrInputField.value.trim();
			}
			
			// Gather Substrates and Products Details (support both old and new class names)
			const detailRows = document.querySelectorAll(`tr.ws-met-row[data-reaction-id="${reactionId}"], tr.detail-item[data-reaction-id="${reactionId}"]`);
			detailRows.forEach((row) => {
				const nameInput = row.querySelector('input[type="text"][name*="NameInput"]');
				const abbrInput = row.querySelector('input[type="text"][name*="AbbrInput"]');
				const detail = {
					name: nameInput.value.trim(),
					abbreviation: abbrInput.value.trim(),
				};

				if (nameInput.name === 'subsNameInput') {
					subsDetails.push(detail);
				} else {
					prodsDetails.push(detail);
				}
			});

			referencesData = Array.from(document.querySelectorAll(`.reference-input[data-reaction-id="${reactionId}"]`)).map((input) => {
				const select = input.previousElementSibling;
				return {
					ref_type: select ? select.value : '',
					info: input.value.trim(),
				};
			});
			
			extLinksData = Array.from(document.querySelectorAll(`.ext-link-item[data-reaction-id="${reactionId}"]`)).map((container) => {
				const select = container.querySelector('.ext-link-type-select');
				const input = container.querySelector('.ext-link-input');
				return {
					ext_link_type: select ? select.value : '',
					info: input ? input.value.trim() : '',
				};
			});
			
			commentsData = Array.from(document.querySelectorAll(`.comment-input[data-reaction-id="${reactionId}"]`)).map((input) => ({
				info: input.value.trim(),
			}));

			// Support both old and new class names for confidence dropdown
			const confidenceDropdown = document.querySelector(`.ws-confidence-select[data-reaction-id="${reactionId}"]`) 
				|| document.querySelector(`.confidencedropdown[data-reaction-id="${reactionId}"]`);
			if (confidenceDropdown) {
				selectedValue = confidenceDropdown.value;
			}
		} else {
			// Use original reaction data (reaction not currently displayed)
			reactionAbbr = reaction.fields.short_name || '';
			
			// Get substrate info from reaction fields
			const subsNames = JSON.parse(reaction.fields.substrates_names || '[]');
			const subsAbbrs = JSON.parse(reaction.fields.subs_found || '[]'); // These might need adjustment based on actual data
			subsDetails = subsNames.map((name, idx) => ({
				name: name,
				abbreviation: '' // Will be filled by backend if needed
			}));
			
			// Get product info from reaction fields
			const prodsNames = JSON.parse(reaction.fields.products_names || '[]');
			prodsDetails = prodsNames.map((name, idx) => ({
				name: name,
				abbreviation: '' // Will be filled by backend if needed
			}));
			
			referencesData = reaction.fields.references || [];
			extLinksData = reaction.fields.ext_links || [];
			commentsData = reaction.fields.comments || [];
		}

		// Returning the updated reaction object
		return {
			pk: reaction.pk,
			description: reaction.fields.description,
			abbreviation: reactionAbbr,
			substrates_info: JSON.stringify(subsDetails),
			products_info: JSON.stringify(prodsDetails),
			references: JSON.stringify(referencesData),
			ext_links: JSON.stringify(extLinksData),
			comments: JSON.stringify(commentsData),
			confidence_score: JSON.stringify(selectedValue),
		};
	}).filter(r => r !== null); // Filter out any null entries

	if (updatedReactions.length === 0) {
		clearInterval(stepInterval);
		setButtonState(false);
		updateVMHModalStatus('error', {
			message: 'No valid reactions to add. Please try again.'
		});
		return;
	}

	fetch(addToVMHUrl, {
		method: 'POST',
		headers: {
			'Content-Type': 'application/json',
			'X-CSRFToken': csrfToken,
		},
		body: JSON.stringify({
			userID: userID,
			reactions: updatedReactions,
		}),
	})
		.then((response) => response.json())
		.then((data) => {
			clearInterval(stepInterval); // Stop updating the steps
			setButtonState(false);
			
			if (data.status === 'success') {
				const reactions = Object.entries(data.rxn_added_info); // Array of [abbr, [rxn_ID, reaction_formula]]
				const metabolites = Object.entries(data.met_added_info); // Array of [abbr, [met_ID, formula, inchiKey]]

				// Update modal with success state
				updateVMHModalStatus('success', {
					reactions: reactions,
					metabolites: metabolites
				});

				// Remove added reactions from DOM and data
				reactions.forEach(([abbr, [id]]) => {
					// Remove from DOM
					const reactionItem = document.querySelector(`.item[data-pk="${id}"]`);
					if (reactionItem) {
						reactionItem.remove();
					}
					// Remove from reactions_active array
					const index = reactions_active.findIndex((r) => r.pk === id);
					if (index !== -1) {
						reactions_active.splice(index, 1);
					}
					// Remove from selected set
					if (window.wsSelectedReactions) {
						window.wsSelectedReactions.delete(String(id));
					}
				});
				
				// Clear workspace selections and update UI
				if (typeof window.clearWorkspaceSelections === 'function') {
					window.clearWorkspaceSelections();
				}
				
				// Clear the details panel
				const detailPanel = document.getElementById('wsAvailableReactionDetails');
				if (detailPanel) {
					detailPanel.innerHTML = `
						<div class="ws-status-card ws-status-success">
							<div class="ws-status-icon"><i class="fas fa-check-circle"></i></div>
							<h3>Successfully Added!</h3>
							<p>${reactions.length} reaction(s) have been added to VMH.</p>
						</div>
					`;
				}
				
				// Show empty state if no more reactions
				if (reactions_active.length === 0) {
					const selectAllRow = document.querySelector('.ws-select-all-row');
					if (selectAllRow) selectAllRow.style.display = 'none';
					
					const listEl = document.getElementById('wsAvailableReactionList');
					if (listEl) {
						listEl.innerHTML = `
							<div class="ws-empty-list">
								<i class="fas fa-check-circle" style="color: #43a047;"></i>
								<p>All reactions have been added!</p>
								<p style="margin-top: 8px; font-size: 12px;">Go to Saved Reactions to add more</p>
							</div>
						`;
					}
				}
			} else {
				// Handle error case with the new modal
				updateVMHModalStatus('error', {
					message: data.message || 'An error occurred while adding to VMH.'
				});
			}
		})
		.catch((error) => {
			clearInterval(stepInterval); // Stop updating the steps
			console.error('Error:', error);
			setButtonState(false);
			
			// Show error in the modal
			updateVMHModalStatus('error', {
				message: 'A network error occurred. Please check your connection and try again.'
			});
		});
}

const confirmAddToVMHBtn = document.getElementById('confirmAddToVMH');
confirmAddToVMHBtn && confirmAddToVMHBtn.addEventListener('click', addToVMH);

// VMH Submission Modal close button
const vmhModalCloseBtn = document.getElementById('vmhModalClose');
vmhModalCloseBtn && vmhModalCloseBtn.addEventListener('click', function() {
	hideVMHModal();
});

// Close VMH modal when clicking overlay
const vmhModalOverlay = document.querySelector('.vmh-submission-overlay');
vmhModalOverlay && vmhModalOverlay.addEventListener('click', function() {
	// Only allow closing if not in loading state
	const footer = document.getElementById('vmhSubmissionFooter');
	if (footer && footer.style.display !== 'none') {
		hideVMHModal();
	}
});

// Close the alert modal when the user clicks on <span> (x)
document.querySelector('.close-alert-btn').addEventListener('click', function () {
	document.getElementById('alertModal').style.display = 'none';

	if (window.pkPendingRemoval !== null) {
		const li = document.querySelector(`#wsAvailableReactionList li[data-pk="${window.pkPendingRemoval}"]`);
		if (li) li.remove();
		window.pkPendingRemoval = null;
	}
});

// Also close the modal if the user clicks anywhere outside of the modal
window.onclick = function (event) {
	if (event.target == document.getElementById('alertModal')) {
		document.getElementById('alertModal').style.display = 'none';

		if (window.pkPendingRemoval !== null) {
			const li = document.querySelector(`#wsAvailableReactionList li[data-pk="${window.pkPendingRemoval}"]`);
			if (li) li.remove();
			window.pkPendingRemoval = null;
		}
	}
};

// Close the modal when the user clicks on <span> (x)
const closeModalBtn = document.getElementsByClassName('close')[0];
closeModalBtn &&
	closeModalBtn.addEventListener('click', function (event) {
		var modal = document.getElementById('reactionModal');
		modal.style.display = 'none';
		window.location.reload();
	});

const closeResponseModalBtn = document.getElementById('responseModalClose');
closeResponseModalBtn &&
	closeResponseModalBtn.addEventListener('click', function (event) {
		if (event.target.id === 'responseModalClose' || event.target.classList.contains('close')) {
			var responseModal = document.getElementById('responseModal');
			responseModal.style.display = 'none';
		}
		window.location.reload();
	});

// Navigation buttons
const homeBtn = document.getElementById('backToHome');
homeBtn && homeBtn.addEventListener('click', () => (window.location.href = '/'));

const wsBtn = document.getElementById('toWorkspace');
wsBtn && wsBtn.addEventListener('click', () => (window.location.href = '/VMH_Workspace/'));

const srBtn = document.getElementById('backToSavedReactions');
srBtn && srBtn.addEventListener('click', () => (window.location.href = '/saved_reactions/'));

// send_to_workspace
const sendBtn = document.getElementById('sendToWorkspace');
sendBtn &&
	sendBtn.addEventListener('click', () => {
		if (checkedReactions.length === 0) {
			showToast('Please select at least one reaction', '#f44336');
			return;
		}
		document.getElementById('confirmSendText').textContent = `Move ${checkedReactions.length} selected reaction(s) to the VMH Workspace?`;
		$('.ui.small.modal#confirmSendModal').modal('show');
	});

const confirmSendBtn = document.getElementById('confirmSendBtn');
confirmSendBtn &&
	confirmSendBtn.addEventListener('click', () => {
		// POST list of ids then go to workspace
		fetch('/send_to_workspace/', {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json',
				'X-CSRFToken': csrfToken,
			},
			body: JSON.stringify({
				userID: userID,
				reactionIds: checkedReactions,
			}),
		})
			.then(() => (window.location.href = '/VMH_Workspace/'))
			.catch((err) => showToast('Could not open workspace', '#f44336'));
	});

const cancelSendBtn = document.getElementById('cancelSendBtn');
cancelSendBtn &&
	cancelSendBtn.addEventListener('click', () => {
		$('.ui.small.modal#confirmSendModal').modal('hide');
	});
