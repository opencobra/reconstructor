document.querySelectorAll('.tab-link').forEach((btn) => {
	btn.addEventListener('click', () => {
		const t = btn.dataset.tab;
		// header buttons
		document.querySelectorAll('.tab-link').forEach((b) => b.classList.remove('active'));
		btn.classList.add('active');
		// containers
		document.querySelectorAll('.tab-content').forEach((c) => c.classList.remove('active'));
		document.getElementById(t).classList.add('active');
	});
});

// sub-tabs (delegated per parent .subtabs block)
document.querySelectorAll('.subtabs').forEach((bar) => {
	bar.addEventListener('click', (e) => {
		if (!e.target.classList.contains('subtab-link')) return;
		const btn = e.target;
		const sub = btn.dataset.subtab;
		const pane = bar.closest('.tab-content');

		// activate clicked sub-tab
		pane.querySelectorAll('.subtab-link').forEach((b) => b.classList.remove('active'));
		btn.classList.add('active');

		// show matching panel
		pane.querySelectorAll('.subtab-content').forEach((p) => p.classList.remove('active'));
		pane.querySelector('#' + sub).classList.add('active');
	});
});

/* global reactions, createSectionHTML, attachDynamicEventListeners,
          displayValidationMessage, validateInputs, setButtonState,
          createExtLinkSelect, createRefSelect, createGeneSectionHTML */

(function () {
	// 1 – draw the list of abbreviations
	const availableListEl = document.getElementById('wsAvailableReactionList');
	const availableDetailEl = document.getElementById('wsAvailableReactionDetails');
	window.pkPendingRemoval = null;
	
	// Track selected reactions for batch operations
	window.wsSelectedReactions = new Set();

	// Show empty state if no reactions
	if (reactions_active.length === 0) {
		const selectAllRow = document.querySelector('.ws-select-all-row');
		if (selectAllRow) selectAllRow.style.display = 'none';
		
		availableListEl.innerHTML = `
			<div class="ws-empty-list">
				<i class="fas fa-inbox"></i>
				<p>No reactions in workspace</p>
				<p style="margin-top: 8px; font-size: 12px;">Go to Saved Reactions and send some reactions here</p>
			</div>
		`;
	}

	// Load available reactions with checkboxes
	reactions_active.forEach((r) => {
		const li = document.createElement('li');
		li.className = 'item';
		li.dataset.pk = r.pk;
		li.innerHTML = `
			<label class="ws-checkbox-container" onclick="event.stopPropagation();">
				<input type="checkbox" class="ws-reaction-checkbox" data-pk="${r.pk}">
				<span class="ws-checkmark"></span>
			</label>
			<span class="ws-reaction-name">${r.fields.short_name}</span>
		`;
		availableListEl.appendChild(li);
	});

	// Handle checkbox changes
	availableListEl.addEventListener('change', (e) => {
		if (e.target.classList.contains('ws-reaction-checkbox')) {
			const pk = e.target.dataset.pk;
			if (e.target.checked) {
				window.wsSelectedReactions.add(pk);
			} else {
				window.wsSelectedReactions.delete(pk);
			}
			updateBatchActionsUI();
		}
	});

	// Select All functionality
	const selectAllCheckbox = document.getElementById('wsSelectAll');
	if (selectAllCheckbox) {
		selectAllCheckbox.addEventListener('change', (e) => {
			const isChecked = e.target.checked;
			document.querySelectorAll('.ws-reaction-checkbox').forEach((cb) => {
				cb.checked = isChecked;
				const pk = cb.dataset.pk;
				if (isChecked) {
					window.wsSelectedReactions.add(pk);
				} else {
					window.wsSelectedReactions.delete(pk);
				}
			});
			updateBatchActionsUI();
		});
	}

	// Update batch actions UI
	function updateBatchActionsUI() {
		const batchActions = document.getElementById('wsBatchActions');
		const selectedCount = document.getElementById('wsSelectedCount');
		const count = window.wsSelectedReactions.size;
		
		if (count > 0) {
			batchActions.style.display = 'flex';
			selectedCount.textContent = `${count} selected`;
		} else {
			batchActions.style.display = 'none';
		}

		// Update select all checkbox state
		const allCheckboxes = document.querySelectorAll('.ws-reaction-checkbox');
		const checkedCheckboxes = document.querySelectorAll('.ws-reaction-checkbox:checked');
		const selectAllCb = document.getElementById('wsSelectAll');
		if (selectAllCb) {
			if (allCheckboxes.length > 0 && checkedCheckboxes.length === allCheckboxes.length) {
				selectAllCb.checked = true;
				selectAllCb.indeterminate = false;
			} else if (checkedCheckboxes.length > 0) {
				selectAllCb.checked = false;
				selectAllCb.indeterminate = true;
			} else {
				selectAllCb.checked = false;
				selectAllCb.indeterminate = false;
			}
		}
	}

	// Make updateBatchActionsUI available globally
	window.updateBatchActionsUI = updateBatchActionsUI;

	// Remove a reaction from the list
	function removeReactionFromList(pk) {
		const listItem = document.querySelector(`#wsAvailableReactionList .item[data-pk="${pk}"]`);
		if (listItem) {
			listItem.classList.add('ws-item-removing');
			setTimeout(() => {
				listItem.remove();
				// Clear the details panel
				availableDetailEl.innerHTML = '<p>Select a reaction to view details</p>';
				// Remove from selected set
				window.wsSelectedReactions.delete(String(pk));
				updateBatchActionsUI();
				// Remove from reactions_active array
				const index = reactions_active.findIndex(r => r.pk === parseInt(pk));
				if (index !== -1) {
					reactions_active.splice(index, 1);
				}
			}, 300);
		}
	}

	// Make removeReactionFromList available globally
	window.removeReactionFromList = removeReactionFromList;

	// Batch add button handler
	const batchAddBtn = document.getElementById('wsBatchAddBtn');
	if (batchAddBtn) {
		batchAddBtn.addEventListener('click', () => {
			if (window.wsSelectedReactions.size === 0) {
				showWorkspaceToast('Please select at least one reaction', 'error');
				return;
			}
			addSelectedToVMH();
		});
	}

	// 2 – on click → render editable form (click on the reaction name, not checkbox)
	availableListEl.addEventListener('click', async (e) => {
		// Only respond to clicks on the item itself or the reaction name, not the checkbox
		const clickedItem = e.target.closest('.item');
		if (!clickedItem) return;
		if (e.target.classList.contains('ws-reaction-checkbox') || e.target.classList.contains('ws-checkmark') || e.target.classList.contains('ws-checkbox-container')) return;
		
		const pk = +clickedItem.dataset.pk;
		const rxn = reactions_active.find((r) => r.pk === pk);
		if (!rxn) return;

		availableDetailEl.innerHTML = '';

		const card = await buildEditableCard(rxn);
		if (card) availableDetailEl.appendChild(card);

		document.querySelectorAll('#wsAvailableReactionList .item').forEach((el) => el.classList.remove('active'));
		clickedItem.classList.add('active');
	});

	// 3 – helper that re-uses modal-builders
	async function buildEditableCard(reaction) {
		let vmhResponse;

		const reactionIndex = reaction.pk;
		const reactionId = reaction.pk;

		document.getElementById('loadingIndicator').style.display = 'flex';
		document.getElementById('loadingText').textContent = 'Gathering Data and (if needed) Generating Abbreviations for Selected Reactions';

		try {
			vmhResponse = await callPrepareAddToVMH([reaction.pk]);
		} catch (error) {
			console.error('Error fetching VMH preparation data: ', error);
			document.getElementById('loadingIndicator').style.display = 'none';
			// Show inline error instead of modal
			const errorCard = document.createElement('div');
			errorCard.className = 'ws-error-card';
			errorCard.innerHTML = `
				<div class="ws-error-icon"><i class="fas fa-exclamation-triangle"></i></div>
				<h3>Error Loading Reaction</h3>
				<p>An error occurred while fetching data from VMH. Please try again later.</p>
				<button class="ws-error-retry-btn" onclick="location.reload()">
					<i class="fas fa-redo"></i> Retry
				</button>
			`;
			return errorCard;
		}

		if (vmhResponse.status === 'error') {
			document.getElementById('loadingIndicator').style.display = 'none';
			
			// Check if it's an "already in VMH" error
			const isAlreadyInVMH = vmhResponse.message && vmhResponse.message.includes('already in VMH');
			
			if (isAlreadyInVMH) {
				// Show a friendly "already in VMH" card instead of removing
				const inVMHCard = document.createElement('div');
				inVMHCard.className = 'ws-status-card ws-status-success';
				inVMHCard.innerHTML = `
					<div class="ws-status-icon"><i class="fas fa-check-circle"></i></div>
					<h3>Already in VMH</h3>
					<p>The reaction <strong>"${reaction.fields.short_name}"</strong> is already in the VMH database.</p>
					<div class="ws-status-actions">
						<button class="ws-status-btn ws-status-btn-secondary" data-remove-pk="${reaction.pk}">
							<i class="fas fa-trash"></i> Remove from List
						</button>
					</div>
				`;
				
				// Add event listener for remove button
				inVMHCard.querySelector('[data-remove-pk]').addEventListener('click', (e) => {
					const pk = e.target.closest('[data-remove-pk]').dataset.removePk;
					removeReactionFromList(pk);
				});
				
				// Mark the reaction item in the list
				const listItem = document.querySelector(`#wsAvailableReactionList .item[data-pk="${reaction.pk}"]`);
				if (listItem) {
					listItem.classList.add('ws-item-in-vmh');
				}
				
				return inVMHCard;
			} else {
				// Show generic error card
				const errorCard = document.createElement('div');
				errorCard.className = 'ws-error-card';
				errorCard.innerHTML = `
					<div class="ws-error-icon"><i class="fas fa-exclamation-circle"></i></div>
					<h3>Unable to Process</h3>
					<p>${vmhResponse.message}</p>
					<button class="ws-error-retry-btn" onclick="location.reload()">
						<i class="fas fa-redo"></i> Retry
					</button>
				`;
				return errorCard;
			}
		}

		document.getElementById('loadingIndicator').style.display = 'none';
		document.getElementById('loadingText').textContent = '';

		subsInVMH = vmhResponse.subs_in_vmh;
		prodsInVMH = vmhResponse.prods_in_vmh;
		subsAbbr = vmhResponse.subs_abbr;
		prodsAbbr = vmhResponse.prods_abbr;
		subsNeedNewNames = vmhResponse.subs_need_new_names;
		prodsNeedNewNames = vmhResponse.prods_need_new_names;
		reactionAbbrs = vmhResponse.reaction_abbrs;

		let substrates_names = JSON.parse(reaction.fields.substrates_names);
		let products_names = JSON.parse(reaction.fields.products_names);
		let subs_comps = JSON.parse(reaction.fields.subs_comps);
		let prods_comps = JSON.parse(reaction.fields.prods_comps);
		let subs_stoich = JSON.parse(reaction.fields.subs_sch);
		let prods_stoich = JSON.parse(reaction.fields.prods_sch);
		let substrates = JSON.parse(reaction.fields.substrates);
		let products = JSON.parse(reaction.fields.products);
		let substrates_types = JSON.parse(reaction.fields.substrates_types);
		let products_types = JSON.parse(reaction.fields.products_types);
		let subsInVMHForReaction = subsInVMH[0];
		let prodsInVMHForReaction = prodsInVMH[0];
		let subsAbbrForReaction = subsAbbr[0];
		let prodsAbbrForReaction = prodsAbbr[0];
		let subsNeedNewNamesForReaction = subsNeedNewNames[0];
		let prodsNeedNewNamesForReaction = prodsNeedNewNames[0];
		let reactionAbbrForReaction = reactionAbbrs[0];
		let confidenceScore = reaction.fields.confidence_score || ' ';

		const listItem = document.createElement('div');
		listItem.className = 'modal-reaction-entry';

		listItem.innerHTML = `
                <div class="reaction-header">

                    <div class="reaction-field">
                        <label>Description:</label>
                        <input type="text" class="reaction-name-input" placeholder="Reaction Name in VMH"
                            value="${reaction.fields.description}" data-reaction-id="${reaction.pk}">
                    </div>

                    <div class="reaction-field">
                        <label>Reaction Abbreviation:</label>
                        <input type="text" class="reaction-abbreviation-input" placeholder="Enter reaction abbreviation"
                            value="${reaction.fields.short_name}" data-reaction-id="${reaction.pk}">
                    </div>

                    <div class="reaction-field">
                        <label for="confidencedropdown-${reactionId}">Confidence Score:</label>
                        <select class="confidencedropdown" id="confidencedropdown-${reaction.pk}" data-reaction-id="${reaction.pk}">
                            <option value=" " ${confidenceScore === ' ' ? 'selected' : ''}>-</option>
                            <option value="1" ${confidenceScore === '1' ? 'selected' : ''}>1</option>
                            <option value="2" ${confidenceScore === '2' ? 'selected' : ''}>2</option>
                            <option value="3" ${confidenceScore === '3' ? 'selected' : ''}>3</option>
                            <option value="4" ${confidenceScore === '4' ? 'selected' : ''}>4</option>
                        </select>
                    </div>

                    <div class="cs-info">
                        <button class="cs-info-button" onclick="toggleInfo()">i</button>
                    </div>
                </div>


                <div class="reaction-details">
                    <p>Substrates:</p>
                    <table>
                        <thead>
                            <tr>
                                <th>Stoichiometry</th>
                                <th>Comp</th>
                                <th>Name</th>
                                <th>Abbreviation</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${substrates_names
															.map(
																(name, index) => `
                                <tr class="detail-item" data-reaction-id="${reaction.pk}" data-tooltip-content="${formatTooltipContent(
																	substrates[index],
																	substrates_types[index],
																	subs_comps[index]
																)}">
                                    <td>
                                        <text>${subs_stoich[index]}</text>
                                    </td>
                                    <td>
                                        <text>${subs_comps[index]}</text>
                                    </td>
                                    <td>
                                        <input type="text" name=subsNameInput placeholder="Name" value="${name}" ${
																	subsInVMHForReaction[index] ? 'readonly' : ''
																}>
                                        ${
																					subsNeedNewNamesForReaction[index]
																						? '<span class="info-icon">&#63;</span><div class="tooltip-content">The name ' +
																						  name +
																						  ' is already assigned in VMH, assign another for this metabolite. Please also check that the metabolite you are adding is not already in VMH.</div>'
																						: ''
																				}
                                    </td>
                                    <td>
                                        <input type="text" name=subsAbbrInput placeholder="Abbreviation" value="${subsAbbrForReaction[index]}" ${
																	subsInVMHForReaction[index] ? 'readonly' : ''
																}>
                                    </td>
                                </tr>
                            `
															)
															.join('')}
                        </tbody>
                    </table>
                </div>
                <div class="reaction-details">
                    <p>Products:</p>
                    <table>
                        <thead>
                            <tr>
                                <th>Stoichiometry</th>
                                <th>Comp</th>
                                <th>Name</th>
                                <th>Abbreviation</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${products_names
															.map(
																(name, index) => `
                                <tr class="detail-item" data-reaction-id="${reaction.pk}" data-tooltip-content="${formatTooltipContent(
																	products[index],
																	products_types[index],
																	prods_comps[index]
																)}">
                                    <td>
                                        <text>${prods_stoich[index]}</text>
                                    </td>
                                    <td>
                                        <text>${prods_comps[index]}</text>
                                    </td>
                                    <td>
                                        <input type="text" name=prodsNameInput placeholder="Name" value="${name}" ${
																	prodsInVMHForReaction[index] ? 'readonly' : ''
																}>
                                        ${
																					prodsNeedNewNamesForReaction[index]
																						? '<span class="info-icon">&#63;</span><div class="tooltip-content">The name ' +
																						  name +
																						  ' is already assigned in VMH, assign another for this metabolite. Please also check that the metabolite you are adding is not already in VMH.</div>'
																						: ''
																				}
                                    </td>
                                    <td>
                                        <input type="text" name=prodsAbbrInput placeholder="Abbreviation" value="${prodsAbbrForReaction[index]}" ${
																	prodsInVMHForReaction[index] ? 'readonly' : ''
																}>
                                    </td>
                                </tr>
                            `
															)
															.join('')}
                        </tbody>
                    </table>
                </div>
                `;
		// Add extra sections
		let references = reaction.fields.references || [];
		let ext_links = reaction.fields.ext_links || [];
		let comments = reaction.fields.comments || [];
		let gene_info = reaction.fields.gene_info || [];

		gene_info = gene_info.map((item) => {
			if (item.info) {
				let first = item.info.split(';')[0].trim();

				if (first.startsWith("GPR: ")) {
					first = first.slice("GPR: ".length);
				}

				item.info = first;
			}
			return item;
		});

		// Build GPR summary string
		const gprItems = gene_info.filter(item => item.info && item.info.trim() !== '').map(item => item.info);
		const gprSummary = gprItems.length > 0 
			? (gprItems.length > 1 
				? gprItems.map(item => `(${item})`).join(' OR ') 
				: gprItems[0]) 
			: '';

		listItem.innerHTML += createSectionHTML('References', 'reference', references, reaction.pk, false, true);
		listItem.innerHTML += createSectionHTML('External Links', 'ext-link', ext_links, reaction.pk, true, false);
		listItem.innerHTML += createSectionHTML('Comments', 'comment', comments, reaction.pk);
		listItem.innerHTML += createSectionHTML('Gene Info', 'gene-info', gene_info, reaction.pk, false, false, true);
		
		// Add GPR summary after Gene Info section
		if (gprSummary) {
			listItem.innerHTML += `<div class="gpr-summary"><span class="gpr-label">GPR:</span> <code class="gpr-formula">${gprSummary}</code></div>`;
		}

		listItem.innerHTML += `<div class="ws-rxn-actions">
        <button class="ui primary button" id="submitVMHBtn" data-pk="${reaction.pk}">Add to VMH</button>
    </div>`;

		return listItem;
	}

	// 4 – Save & Submit stubs
	availableDetailEl.addEventListener('click', (e) => {
		if (e.target.id === 'submitVMHBtn') {
			const pk = +e.target.dataset.pk;
			// reuse addToVMH() – but restrict to this one pk
			addToVMHforSingle(pk);
		}
	});

	function addToVMHforSingle(pk) {
		window.checkedReactions = [String(pk)];
		addToVMH(); // comes from handleAdd2VMH.js
	}

	// Batch add selected reactions to VMH
	async function addSelectedToVMH() {
		const selectedPks = Array.from(window.wsSelectedReactions);
		if (selectedPks.length === 0) {
			showWorkspaceToast('Please select at least one reaction', 'error');
			return;
		}

		// Show confirmation modal with selected reactions
		const selectedReactions = selectedPks.map(pk => {
			const rxn = reactions_active.find(r => r.pk === parseInt(pk));
			return rxn ? rxn.fields.short_name : pk;
		});

		const confirmMessage = `Add ${selectedPks.length} reaction(s) to VMH?\n\n${selectedReactions.join(', ')}`;
		
		if (!confirm(confirmMessage)) {
			return;
		}

		// Use the existing addToVMH flow with multiple reactions
		window.checkedReactions = selectedPks;
		addToVMH(); // This will process all selected reactions
	}

	// Make addSelectedToVMH available 
	window.addSelectedToVMH = addSelectedToVMH;
})();

const modalList = document.getElementById('wsAvailableReactionDetails');

modalList.addEventListener('click', function (e) {
	// Handle the "Add" button clicks using event delegation
	if (e.target && e.target.matches('.add-reference, .add-ext-link, .add-comment')) {
		const type = e.target.classList.contains('add-reference') ? 'reference' : e.target.classList.contains('add-ext-link') ? 'ext-link' : 'comment';
		const parentSection = e.target.parentNode;
		const reactionId = e.target.getAttribute('data-reaction-id');
		const newItem = document.createElement('div');
		newItem.className = `${type}-item`;
		const items = parentSection.querySelectorAll(`.${type}-item`);
		const newIndex = items.length; // Calculate new index based on existing items
		newItem.setAttribute('data-reaction-id', reactionId);
		newItem.setAttribute('data-index', newIndex);
		newItem.innerHTML = type === 'ext-link' ? createExtLinkSelect({}, reactionId, newIndex) : '';
		newItem.innerHTML += type === 'reference' ? createRefSelect({}, reactionId, newIndex) : '';
		newItem.innerHTML += `
            <input type="text" class="${type}-input" placeholder="Enter ${type}" value="" data-reaction-id="${reactionId}" data-index="${newIndex}">
            <button class="remove-${type}" data-reaction-id="${reactionId}" data-index="${newIndex}">Remove</button>
        `;
		parentSection.insertBefore(newItem, e.target);
	}

	// Handle the "Remove" button clicks using event delegation
	if (e.target && e.target.matches('.remove-reference, .remove-ext-link, .remove-comment,.remove-gene-info')) {
		const isGeneInfo = e.target.matches('.remove-gene-info');
		e.target.parentElement.remove();
		
		// Update GPR summary if a gene info was removed
		if (isGeneInfo) {
			updateGPRSummary();
		}
	}
});

// Function to update GPR summary based on current gene-info inputs
function updateGPRSummary() {
	const gprSummaryEl = document.querySelector('.gpr-summary');
	const geneInfoInputs = document.querySelectorAll('.gene-info-input');
	
	const gprItems = [];
	geneInfoInputs.forEach(input => {
		const value = input.value.trim();
		if (value) {
			gprItems.push(value);
		}
	});
	
	if (gprItems.length > 0) {
		// Add brackets around each item when joining with OR (only if multiple items)
		const gprText = gprItems.length > 1 
			? gprItems.map(item => `(${item})`).join(' OR ')
			: gprItems[0];
		if (gprSummaryEl) {
			gprSummaryEl.innerHTML = `<span class="gpr-label">GPR:</span> <code class="gpr-formula">${gprText}</code>`;
			gprSummaryEl.style.display = '';
		}
	} else {
		if (gprSummaryEl) {
			gprSummaryEl.style.display = 'none';
		}
	}
}

const parsedAddedReactions = reactions_added;

const tableBody = document.getElementById('addedReactionTableBody');
if (parsedAddedReactions.length === 0) {
	const row = document.createElement('tr');
	row.innerHTML = `<td colspan="3" style="text-align: center; color: #999;">No reactions have been added yet.</td>`;
	tableBody.appendChild(row);
} else {
	parsedAddedReactions.forEach((entry) => {
		const f = entry.fields;
		const row = document.createElement('tr');
		row.innerHTML = `
			<td>${f.reaction_abbr}</td>
			<td style="white-space: pre-wrap">${f.reaction_formula}</td>
			<td>${new Date(f.created_at).toLocaleString()}</td>
		`;
		tableBody.appendChild(row);
	});
}

// Toast notification function for workspace
function showWorkspaceToast(message, type = 'info') {
	// Remove existing toasts
	const existingToasts = document.querySelectorAll('.ws-toast');
	existingToasts.forEach(t => t.remove());

	const toast = document.createElement('div');
	toast.className = `ws-toast ws-toast-${type}`;
	toast.innerHTML = `
		<i class="fas ${type === 'success' ? 'fa-check-circle' : type === 'error' ? 'fa-exclamation-circle' : 'fa-info-circle'}"></i>
		<span>${message}</span>
	`;
	document.body.appendChild(toast);

	// Animate in
	setTimeout(() => toast.classList.add('show'), 10);

	// Auto remove after 4 seconds
	setTimeout(() => {
		toast.classList.remove('show');
		setTimeout(() => toast.remove(), 300);
	}, 4000);
}

// Clear selection after successful add (hook into existing addToVMH success)
window.clearWorkspaceSelections = function() {
	window.wsSelectedReactions.clear();
	document.querySelectorAll('.ws-reaction-checkbox').forEach(cb => {
		cb.checked = false;
	});
	if (typeof window.updateBatchActionsUI === 'function') {
		window.updateBatchActionsUI();
	}
};

/* ========== Workspace State Persistence ========== */

// Save workspace state to localStorage
function saveWorkspaceState() {
	const state = {
		activeTab: document.querySelector('.tab-link.active')?.dataset.tab,
		selectedReactionPk: document.querySelector('#wsAvailableReactionList .item.active')?.dataset.pk,
		checkedReactions: Array.from(window.wsSelectedReactions || []),
		timestamp: Date.now()
	};
	localStorage.setItem('vmhWorkspaceState', JSON.stringify(state));
}

// Restore workspace state from localStorage
function restoreWorkspaceState() {
	try {
		const stateJson = localStorage.getItem('vmhWorkspaceState');
		if (!stateJson) return;
		
		const state = JSON.parse(stateJson);
		
		// Only restore if state is less than 1 hour old
		if (Date.now() - state.timestamp > 3600000) {
			localStorage.removeItem('vmhWorkspaceState');
			return;
		}

		// Restore active tab
		if (state.activeTab) {
			const tabBtn = document.querySelector(`.tab-link[data-tab="${state.activeTab}"]`);
			if (tabBtn) tabBtn.click();
		}

		// Restore checked reactions
		if (state.checkedReactions && state.checkedReactions.length > 0) {
			state.checkedReactions.forEach(pk => {
				const checkbox = document.querySelector(`.ws-reaction-checkbox[data-pk="${pk}"]`);
				if (checkbox) {
					checkbox.checked = true;
					window.wsSelectedReactions.add(pk);
				}
			});
			if (typeof window.updateBatchActionsUI === 'function') {
				window.updateBatchActionsUI();
			}
		}

		// Restore selected reaction (click to load details)
		if (state.selectedReactionPk) {
			const reactionItem = document.querySelector(`#wsAvailableReactionList .item[data-pk="${state.selectedReactionPk}"]`);
			if (reactionItem) {
				// Small delay to ensure DOM is ready
				setTimeout(() => {
					const nameSpan = reactionItem.querySelector('.ws-reaction-name');
					if (nameSpan) nameSpan.click();
				}, 100);
			}
		}

		showWorkspaceToast('Workspace state restored', 'info');
	} catch (e) {
		console.error('Error restoring workspace state:', e);
	}
}

// Auto-save state on important actions
function initStateTracking() {
	// Save when tab changes
	document.querySelectorAll('.tab-link').forEach(btn => {
		btn.addEventListener('click', () => {
			setTimeout(saveWorkspaceState, 100);
		});
	});

	// Save when a reaction is selected
	const reactionList = document.getElementById('wsAvailableReactionList');
	if (reactionList) {
		reactionList.addEventListener('click', (e) => {
			if (e.target.closest('.item') && !e.target.classList.contains('ws-reaction-checkbox')) {
				setTimeout(saveWorkspaceState, 100);
			}
		});
	}

	// Save when checkboxes change
	document.addEventListener('change', (e) => {
		if (e.target.classList.contains('ws-reaction-checkbox') || e.target.id === 'wsSelectAll') {
			setTimeout(saveWorkspaceState, 100);
		}
	});

	// Save before page unload
	window.addEventListener('beforeunload', saveWorkspaceState);
}

// Initialize state tracking and restore state
document.addEventListener('DOMContentLoaded', () => {
	initStateTracking();
	// Small delay to let the list render first
	setTimeout(restoreWorkspaceState, 200);
});
