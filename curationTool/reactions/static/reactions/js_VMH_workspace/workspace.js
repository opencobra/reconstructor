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

/**
 * VMH Preparation Cache Manager
 * 
 * Handles caching of reaction preparation data to avoid redundant API calls.
 * Provides methods for storing, retrieving, and invalidating cached data.
 */
class VMHPrepCache {
	constructor(options = {}) {
		this.cache = new Map();
		this.maxAge = options.maxAge || 180 * 60 * 1000; // Default: 180 minutes
		this.maxSize = options.maxSize || 500; // Maximum number of cached items
		this.storageKey = 'vmh_prep_cache';
		this.useSessionStorage = options.useSessionStorage !== false;
		
		// Try to restore cache from sessionStorage on initialization
		this._restoreFromStorage();
	}

	/**
	 * Generate a cache key for a reaction
	 * @param {number} reactionPk - The reaction primary key
	 * @returns {string} Cache key
	 */
	_getCacheKey(reactionPk) {
		return `reaction_${reactionPk}`;
	}

	/**
	 * Check if a cached entry is still valid (not expired)
	 * @param {Object} entry - Cache entry with timestamp
	 * @returns {boolean} True if valid, false if expired
	 */
	_isValid(entry) {
		if (!entry || !entry.timestamp) return false;
		return (Date.now() - entry.timestamp) < this.maxAge;
	}

	/**
	 * Store data in the cache
	 * @param {number} reactionPk - The reaction primary key
	 * @param {Object} data - The data to cache
	 */
	set(reactionPk, data) {
		const key = this._getCacheKey(reactionPk);
		
		// Enforce max size by removing oldest entries
		if (this.cache.size >= this.maxSize) {
			const oldestKey = this.cache.keys().next().value;
			this.cache.delete(oldestKey);
		}

		const entry = {
			data: data,
			timestamp: Date.now()
		};
		
		this.cache.set(key, entry);
		this._persistToStorage();
		
		console.debug(`[VMHPrepCache] Cached data for reaction ${reactionPk}`);
	}

	/**
	 * Retrieve data from the cache
	 * @param {number} reactionPk - The reaction primary key
	 * @returns {Object|null} Cached data or null if not found/expired
	 */
	get(reactionPk) {
		const key = this._getCacheKey(reactionPk);
		const entry = this.cache.get(key);

		if (!entry) {
			console.debug(`[VMHPrepCache] Cache miss for reaction ${reactionPk}`);
			return null;
		}

		if (!this._isValid(entry)) {
			console.debug(`[VMHPrepCache] Cache expired for reaction ${reactionPk}`);
			this.cache.delete(key);
			this._persistToStorage();
			return null;
		}

		console.debug(`[VMHPrepCache] Cache hit for reaction ${reactionPk}`);
		return entry.data;
	}

	/**
	 * Check if valid cached data exists for a reaction
	 * @param {number} reactionPk - The reaction primary key
	 * @returns {boolean} True if valid cache exists
	 */
	has(reactionPk) {
		const key = this._getCacheKey(reactionPk);
		const entry = this.cache.get(key);
		return entry && this._isValid(entry);
	}

	/**
	 * Invalidate (remove) cached data for a specific reaction
	 * @param {number} reactionPk - The reaction primary key
	 */
	invalidate(reactionPk) {
		const key = this._getCacheKey(reactionPk);
		this.cache.delete(key);
		this._persistToStorage();
		console.debug(`[VMHPrepCache] Invalidated cache for reaction ${reactionPk}`);
	}

	/**
	 * Invalidate all cached data
	 */
	invalidateAll() {
		this.cache.clear();
		this._persistToStorage();
		console.debug('[VMHPrepCache] Invalidated all cache entries');
	}

	/**
	 * Persist cache to sessionStorage for page refreshes
	 * @private
	 */
	_persistToStorage() {
		if (!this.useSessionStorage) return;
		
		try {
			const serializable = {};
			this.cache.forEach((value, key) => {
				serializable[key] = value;
			});
			sessionStorage.setItem(this.storageKey, JSON.stringify(serializable));
		} catch (e) {
			console.warn('[VMHPrepCache] Failed to persist cache to sessionStorage:', e);
		}
	}

	/**
	 * Restore cache from sessionStorage
	 * @private
	 */
	_restoreFromStorage() {
		if (!this.useSessionStorage) return;
		
		try {
			const stored = sessionStorage.getItem(this.storageKey);
			if (stored) {
				const parsed = JSON.parse(stored);
				Object.entries(parsed).forEach(([key, entry]) => {
					// Only restore valid (non-expired) entries
					if (this._isValid(entry)) {
						this.cache.set(key, entry);
					}
				});
				console.debug(`[VMHPrepCache] Restored ${this.cache.size} entries from sessionStorage`);
			}
		} catch (e) {
			console.warn('[VMHPrepCache] Failed to restore cache from sessionStorage:', e);
		}
	}

	/**
	 * Get cache statistics
	 * @returns {Object} Cache stats
	 */
	getStats() {
		let validCount = 0;
		let expiredCount = 0;
		
		this.cache.forEach((entry) => {
			if (this._isValid(entry)) {
				validCount++;
			} else {
				expiredCount++;
			}
		});

		return {
			total: this.cache.size,
			valid: validCount,
			expired: expiredCount,
			maxSize: this.maxSize,
			maxAgeMs: this.maxAge
		};
	}
}

// Create global cache instance
window.vmhPrepCache = new VMHPrepCache({
	maxAge: 10 * 60 * 1000, // 10 minutes
	maxSize: 50,
	useSessionStorage: true
});

(function () {
	// 1 – draw the list of abbreviations
	const availableListEl = document.getElementById('wsAvailableReactionList');
	const availableDetailEl = document.getElementById('wsAvailableReactionDetails');
	window.pkPendingRemoval = null;
	
	// Track selected reactions for batch operations
	window.wsSelectedReactions = new Set();

	function readBalanceValue(rawValue) {
		if (rawValue === null || rawValue === undefined || rawValue === '') return null;

		if (Array.isArray(rawValue)) {
			return readBalanceValue(rawValue[0]);
		}

		if (typeof rawValue === 'boolean') return rawValue;

		if (typeof rawValue === 'number') {
			if (rawValue === 1) return true;
			if (rawValue === 0) return false;
			return null;
		}

		if (typeof rawValue === 'string') {
			const normalized = rawValue.trim().toLowerCase();
			if (normalized === 'true') return true;
			if (normalized === 'false') return false;

			try {
				return readBalanceValue(JSON.parse(rawValue));
			} catch (error) {
				return null;
			}
		}

		return null;
	}

	function formatBalanceIssues(issues) {
		if (issues.length === 1) return issues[0];
		return `${issues.slice(0, -1).join(', ')} and ${issues[issues.length - 1]}`;
	}

	function getBalanceWarning(reaction) {
		const balancedCount = readBalanceValue(reaction.fields.balanced_count);
		const balancedCharge = readBalanceValue(reaction.fields.balanced_charge);
		const issues = [];

		if (balancedCount === false) issues.push('atom count');
		if (balancedCharge === false) issues.push('charge');

		if (issues.length === 0) return null;

		return {
			issues,
			label: formatBalanceIssues(issues)
		};
	}

	function createBalanceWarningHTML(reaction) {
		const warning = getBalanceWarning(reaction);
		if (!warning) return '';

		return `
			<div class="ws-balance-warning" role="status">
				<div class="ws-balance-warning-icon">
					<i class="fas fa-exclamation-triangle"></i>
				</div>
				<div class="ws-balance-warning-copy">
					<div class="ws-balance-warning-title">Balance warning</div>
					<div class="ws-balance-warning-text">
						This reaction is unbalanced by ${warning.label}. You can still add it to VMH, but please review the details.
					</div>
				</div>
			</div>
		`;
	}

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
		const balanceWarning = getBalanceWarning(r);
		const li = document.createElement('li');
		li.className = `item${balanceWarning ? ' ws-item-balance-warning' : ''}`;
		li.dataset.pk = r.pk;
		li.innerHTML = `
			<label class="ws-checkbox-container" onclick="event.stopPropagation();">
				<input type="checkbox" class="ws-reaction-checkbox" data-pk="${r.pk}">
				<span class="ws-checkmark"></span>
			</label>
			<span class="ws-reaction-name">${r.fields.short_name}</span>
			${balanceWarning ? `
				<span class="ws-list-balance-warning" title="Unbalanced by ${balanceWarning.label}">
					<i class="fas fa-exclamation-triangle"></i>
					<span>Unbalanced</span>
				</span>
			` : ''}
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
				// Invalidate cache for this reaction
				if (window.vmhPrepCache) {
					window.vmhPrepCache.invalidate(pk);
				}
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
		let fromCache = false;

		const reactionIndex = reaction.pk;
		const reactionId = reaction.pk;

		// Check cache first
		const cachedData = window.vmhPrepCache.get(reaction.pk);
		
		if (cachedData) {
			// Use cached data - no need to show loading indicator
			vmhResponse = cachedData;
			fromCache = true;
			console.debug(`[Workspace] Using cached data for reaction ${reaction.pk}`);
		} else {
			// No cache - fetch from server
			document.getElementById('loadingIndicator').style.display = 'flex';
			document.getElementById('loadingText').textContent = 'Gathering Data and (if needed) Generating Abbreviations for Selected Reactions';

			try {
				vmhResponse = await callPrepareAddToVMH([reaction.pk]);
				
				// Cache successful responses (not errors)
				if (vmhResponse && vmhResponse.status === 'success') {
					window.vmhPrepCache.set(reaction.pk, vmhResponse);
				}
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
		}

		if (vmhResponse.status === 'error') {
			if (!fromCache) {
				document.getElementById('loadingIndicator').style.display = 'none';
			}
			
			// Invalidate cache for this reaction since it has an error status
			window.vmhPrepCache.invalidate(reaction.pk);
			
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

		// Hide loading indicator only if we showed it (not from cache)
		if (!fromCache) {
			document.getElementById('loadingIndicator').style.display = 'none';
			document.getElementById('loadingText').textContent = '';
		}

		subsInVMH = vmhResponse.subs_in_vmh;
		prodsInVMH = vmhResponse.prods_in_vmh;
		subsAbbr = vmhResponse.subs_abbr;
		prodsAbbr = vmhResponse.prods_abbr;
		subsNeedNewNames = vmhResponse.subs_need_new_names;
		prodsNeedNewNames = vmhResponse.prods_need_new_names;
		reactionAbbrs = vmhResponse.reaction_abbrs;

		// Check if there's a cached form state to restore
		const cachedFormState = vmhResponse._formState || null;

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
		
		// Use cached form state values if available, otherwise use reaction fields
		let displayDescription = cachedFormState ? cachedFormState.description : reaction.fields.description;
		let displayAbbreviation = cachedFormState ? cachedFormState.abbreviation : reaction.fields.short_name;
		let confidenceScore = cachedFormState ? cachedFormState.confidence_score : (reaction.fields.confidence_score || ' ');

		const listItem = document.createElement('div');
		listItem.className = 'ws-reaction-card';

		listItem.innerHTML = `
			<div class="ws-card-header">
				<div class="ws-card-title-row">
					<h3 class="ws-card-title">${reaction.fields.short_name}</h3>
					<div class="ws-confidence-badge" data-score="${confidenceScore}">
						<select class="ws-confidence-select" id="confidencedropdown-${reaction.pk}" data-reaction-id="${reaction.pk}">
							<option value=" " ${confidenceScore === ' ' ? 'selected' : ''}>–</option>
							<option value="1" ${confidenceScore === '1' ? 'selected' : ''}>1</option>
							<option value="2" ${confidenceScore === '2' ? 'selected' : ''}>2</option>
							<option value="3" ${confidenceScore === '3' ? 'selected' : ''}>3</option>
							<option value="4" ${confidenceScore === '4' ? 'selected' : ''}>4</option>
						</select>
						<button class="ws-cs-info-btn" onclick="toggleInfo()" title="What is confidence score?">
							<i class="fas fa-info-circle"></i>
						</button>
					</div>
				</div>
			</div>
			
			<div class="ws-card-body">
				${createBalanceWarningHTML(reaction)}
				<div class="ws-form-section">
					<div class="ws-form-row">
						<div class="ws-form-group ws-form-group-lg">
							<label class="ws-form-label">Description</label>
							<input type="text" class="ws-form-input reaction-name-input" 
								placeholder="Enter reaction description" 
								value="${displayDescription}" 
								data-reaction-id="${reaction.pk}">
						</div>
						<div class="ws-form-group">
							<label class="ws-form-label">Abbreviation</label>
							<input type="text" class="ws-form-input reaction-abbreviation-input" 
								placeholder="e.g., RXN001" 
								value="${displayAbbreviation}" 
								data-reaction-id="${reaction.pk}">
						</div>
					</div>
				</div>

				<div class="ws-metabolites-section">
					<div class="ws-metabolites-panel">
						<div class="ws-panel-header">
							<span class="ws-panel-title">Substrates</span>
						</div>
						<div class="ws-metabolites-table-wrapper">
							<table class="ws-metabolites-table">
								<thead>
									<tr>
										<th class="ws-col-stoich">Stoich</th>
										<th class="ws-col-comp">Comp</th>
										<th class="ws-col-name">Name</th>
										<th class="ws-col-abbr">Abbreviation</th>
									</tr>
								</thead>
								<tbody>
									${substrates_names.map((name, index) => `
										<tr class="ws-met-row ${subsInVMHForReaction[index] ? 'ws-met-vmh' : 'ws-met-new'}" 
											data-reaction-id="${reaction.pk}" 
											data-tooltip-content="${formatTooltipContent(substrates[index], substrates_types[index], subs_comps[index])}">
											<td class="ws-col-stoich"><span class="ws-stoich-value">${subs_stoich[index]}</span></td>
											<td class="ws-col-comp"><span class="ws-comp-badge">${subs_comps[index]}</span></td>
											<td class="ws-col-name">
												<div class="ws-name-cell">
													<input type="text" name="subsNameInput" class="ws-met-input ${subsInVMHForReaction[index] ? 'ws-input-readonly' : ''}" 
														placeholder="Name" value="${name}" ${subsInVMHForReaction[index] ? 'readonly' : ''}>
													${subsNeedNewNamesForReaction[index] ? `
														<span class="ws-warning-icon" title="Name already in VMH">
															<i class="fas fa-exclamation-triangle"></i>
														</span>
													` : ''}
													${subsInVMHForReaction[index] ? '<span class="ws-vmh-tag">VMH</span>' : ''}
												</div>
											</td>
											<td class="ws-col-abbr">
												<input type="text" name="subsAbbrInput" class="ws-met-input ${subsInVMHForReaction[index] ? 'ws-input-readonly' : ''}" 
													placeholder="Abbr" value="${subsAbbrForReaction[index]}" ${subsInVMHForReaction[index] ? 'readonly' : ''}>
											</td>
										</tr>
									`).join('')}
								</tbody>
							</table>
						</div>
					</div>

					<div class="ws-reaction-arrow">
						<i class="fas fa-long-arrow-alt-right"></i>
					</div>

					<div class="ws-metabolites-panel">
						<div class="ws-panel-header">
							<span class="ws-panel-title">Products</span>
						</div>
						<div class="ws-metabolites-table-wrapper">
							<table class="ws-metabolites-table">
								<thead>
									<tr>
										<th class="ws-col-stoich">Stoich</th>
										<th class="ws-col-comp">Comp</th>
										<th class="ws-col-name">Name</th>
										<th class="ws-col-abbr">Abbreviation</th>
									</tr>
								</thead>
								<tbody>
									${products_names.map((name, index) => `
										<tr class="ws-met-row ${prodsInVMHForReaction[index] ? 'ws-met-vmh' : 'ws-met-new'}" 
											data-reaction-id="${reaction.pk}" 
											data-tooltip-content="${formatTooltipContent(products[index], products_types[index], prods_comps[index])}">
											<td class="ws-col-stoich"><span class="ws-stoich-value">${prods_stoich[index]}</span></td>
											<td class="ws-col-comp"><span class="ws-comp-badge">${prods_comps[index]}</span></td>
											<td class="ws-col-name">
												<div class="ws-name-cell">
													<input type="text" name="prodsNameInput" class="ws-met-input ${prodsInVMHForReaction[index] ? 'ws-input-readonly' : ''}" 
														placeholder="Name" value="${name}" ${prodsInVMHForReaction[index] ? 'readonly' : ''}>
													${prodsNeedNewNamesForReaction[index] ? `
														<span class="ws-warning-icon" title="Name already in VMH">
															<i class="fas fa-exclamation-triangle"></i>
														</span>
													` : ''}
													${prodsInVMHForReaction[index] ? '<span class="ws-vmh-tag">VMH</span>' : ''}
												</div>
											</td>
											<td class="ws-col-abbr">
												<input type="text" name="prodsAbbrInput" class="ws-met-input ${prodsInVMHForReaction[index] ? 'ws-input-readonly' : ''}" 
													placeholder="Abbr" value="${prodsAbbrForReaction[index]}" ${prodsInVMHForReaction[index] ? 'readonly' : ''}>
											</td>
										</tr>
									`).join('')}
								</tbody>
							</table>
						</div>
					</div>
				</div>
		`;
		// Add extra sections - use cached form state if available
		let references = cachedFormState ? cachedFormState.references : (reaction.fields.references || []);
		let ext_links = cachedFormState ? cachedFormState.ext_links : (reaction.fields.ext_links || []);
		let comments = cachedFormState ? cachedFormState.comments : (reaction.fields.comments || []);
		let gene_info = cachedFormState ? cachedFormState.gene_info : (reaction.fields.gene_info || []);

		// Only process gene_info if it came from the original reaction fields (not cached)
		if (!cachedFormState) {
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
		}

		// Build GPR summary string
		const gprItems = gene_info.filter(item => item.info && item.info.trim() !== '').map(item => item.info);
		const gprSummary = gprItems.length > 0 
			? (gprItems.length > 1 
				? gprItems.map(item => `(${item})`).join(' OR ') 
				: gprItems[0]) 
			: '';

		// Build annotation cards as a single string
		const annotationCards = 
			createSectionHTMLNew('References', 'reference', references, reaction.pk, false, true) +
			createSectionHTMLNew('External Links', 'ext-link', ext_links, reaction.pk, true, false) +
			createSectionHTMLNew('Comments', 'comment', comments, reaction.pk);

		// Build the entire annotations section as ONE complete string to avoid browser auto-closing unclosed tags
		listItem.innerHTML += `
				<div class="ws-annotations-section">
					<div class="ws-section-header">
						<h4 class="ws-section-title"><i class="fas fa-tags"></i> Annotations & Metadata</h4>
					</div>
					<div class="ws-annotations-grid">
						${annotationCards}
					</div>
				</div>
		`;
		
		// Gene Info section with GPR summary
		listItem.innerHTML += `
				<div class="ws-gene-section">
					<div class="ws-section-header">
						<h4 class="ws-section-title"><i class="fas fa-dna"></i> Gene Association</h4>
					</div>
					${gprSummary ? `
						<div class="ws-gpr-display">
							<span class="ws-gpr-label">GPR Rule:</span>
							<code class="ws-gpr-code">${gprSummary}</code>
						</div>
					` : `
						<div class="ws-gpr-empty">
							<i class="fas fa-info-circle"></i>
							<span>No gene associations defined</span>
						</div>
					`}
					<div class="ws-gene-items">
						${gene_info.map((item, index) => `
							<div class="gene-info-item ws-gene-chip" data-reaction-id="${reaction.pk}" data-index="${index}">
								<input type="text" class="gene-info-input" value="${item.info}" 
									data-reaction-id="${reaction.pk}" data-index="${index}" readonly>
								<button class="ws-chip-remove remove-gene-info" data-reaction-id="${reaction.pk}" data-index="${index}">
									<i class="fas fa-times"></i>
								</button>
							</div>
						`).join('')}
					</div>
				</div>
			</div>

			<div class="ws-card-footer">
				<div class="ws-footer-left">
					<button class="ws-footer-btn ws-reset-btn" id="resetChangesBtn" data-pk="${reaction.pk}" title="Discard all changes and reload from server">
						<i class="fas fa-undo"></i>
						<span>Reset</span>
					</button>
					<button class="ws-footer-btn ws-save-local-btn" id="saveLocalBtn" data-pk="${reaction.pk}" title="Save changes to local database">
						<i class="fas fa-save"></i>
						<span>Save Draft</span>
					</button>
				</div>
				<button class="ws-submit-btn" id="submitVMHBtn" data-pk="${reaction.pk}">
					<i class="fas fa-cloud-upload-alt"></i>
					<span>Add to VMH</span>
				</button>
			</div>
		`;

		return listItem;
	}

	/**
	 * Collect current form state for a reaction
	 * @param {number} reactionPk - The reaction primary key
	 * @returns {Object|null} Current form state or null if form not found
	 */
	function collectFormState(reactionPk) {
		const nameInput = document.querySelector(`.reaction-name-input[data-reaction-id="${reactionPk}"]`);
		if (!nameInput) return null;

		const abbrInput = document.querySelector(`.reaction-abbreviation-input[data-reaction-id="${reactionPk}"]`);
		const confidenceSelect = document.querySelector(`.ws-confidence-select[data-reaction-id="${reactionPk}"]`);

		// Collect substrate details
		const subsRows = document.querySelectorAll(`tr.ws-met-row[data-reaction-id="${reactionPk}"] input[name="subsNameInput"]`);
		const subsDetails = [];
		subsRows.forEach((nameInput) => {
			const row = nameInput.closest('tr');
			const abbrInput = row.querySelector('input[name="subsAbbrInput"]');
			subsDetails.push({
				name: nameInput.value.trim(),
				abbreviation: abbrInput ? abbrInput.value.trim() : ''
			});
		});

		// Collect product details
		const prodsRows = document.querySelectorAll(`tr.ws-met-row[data-reaction-id="${reactionPk}"] input[name="prodsNameInput"]`);
		const prodsDetails = [];
		prodsRows.forEach((nameInput) => {
			const row = nameInput.closest('tr');
			const abbrInput = row.querySelector('input[name="prodsAbbrInput"]');
			prodsDetails.push({
				name: nameInput.value.trim(),
				abbreviation: abbrInput ? abbrInput.value.trim() : ''
			});
		});

		// Collect references
		const references = [];
		document.querySelectorAll(`.reference-item[data-reaction-id="${reactionPk}"]`).forEach((item) => {
			const select = item.querySelector('.ref-type-select');
			const input = item.querySelector('.reference-input');
			if (input && input.value.trim()) {
				references.push({
					ref_type: select ? select.value : 'DOI',
					info: input.value.trim()
				});
			}
		});

		// Collect external links
		const ext_links = [];
		document.querySelectorAll(`.ext-link-item[data-reaction-id="${reactionPk}"]`).forEach((item) => {
			const select = item.querySelector('.ext-link-type-select');
			const input = item.querySelector('.ext-link-input');
			if (input && input.value.trim()) {
				ext_links.push({
					ext_link_type: select ? select.value : 'KEGG reaction',
					info: input.value.trim()
				});
			}
		});

		// Collect comments
		const comments = [];
		document.querySelectorAll(`.comment-item[data-reaction-id="${reactionPk}"]`).forEach((item) => {
			const input = item.querySelector('.comment-input');
			if (input && input.value.trim()) {
				comments.push({
					info: input.value.trim()
				});
			}
		});

		// Collect gene info
		const gene_info = [];
		document.querySelectorAll(`.gene-info-item[data-reaction-id="${reactionPk}"]`).forEach((item) => {
			const input = item.querySelector('.gene-info-input');
			if (input && input.value.trim()) {
				gene_info.push({
					info: input.value.trim()
				});
			}
		});

		return {
			description: nameInput.value.trim(),
			abbreviation: abbrInput ? abbrInput.value.trim() : '',
			confidence_score: confidenceSelect ? confidenceSelect.value : '',
			substrates_info: subsDetails,
			products_info: prodsDetails,
			references: references,
			ext_links: ext_links,
			comments: comments,
			gene_info: gene_info
		};
	}

	/**
	 * Update the cache with current form state
	 * @param {number} reactionPk - The reaction primary key
	 */
	function updateCacheFromForm(reactionPk) {
		const cachedData = window.vmhPrepCache.get(reactionPk);
		if (!cachedData) return;

		const formState = collectFormState(reactionPk);
		if (!formState) return;

		// Update the reaction object in reactions_active with form state
		const reaction = reactions_active.find(r => r.pk === reactionPk);
		if (reaction) {
			reaction.fields.description = formState.description;
			reaction.fields.short_name = formState.abbreviation;
			reaction.fields.confidence_score = formState.confidence_score;
			reaction.fields.references = formState.references;
			reaction.fields.ext_links = formState.ext_links;
			reaction.fields.comments = formState.comments;
			reaction.fields.gene_info = formState.gene_info;
			
			// Store form state in cache for persistence
			cachedData._formState = formState;
			window.vmhPrepCache.set(reactionPk, cachedData);
		}

		console.debug(`[Workspace] Updated cache for reaction ${reactionPk}`);
	}

	/**
	 * Reset reaction to original state from server
	 * @param {number} reactionPk - The reaction primary key
	 */
	async function resetReactionChanges(reactionPk) {
		// Show loading state
		const availableDetailEl = document.getElementById('wsAvailableReactionDetails');
		availableDetailEl.innerHTML = '<div class="ws-loading-inline"><i class="fas fa-spinner fa-spin"></i> Resetting...</div>';

		try {
			// Invalidate cache to force fresh fetch
			window.vmhPrepCache.invalidate(reactionPk);
			
			// Fetch fresh reaction data from server (basic fields)
			const [reactionResponse, detailsResponse] = await Promise.all([
				fetch(`/get_reaction/${reactionPk}/`),
				fetch('/get_reaction_details/', {
					method: 'POST',
					headers: {
						'Content-Type': 'application/json',
						'X-CSRFToken': csrfToken
					},
					body: JSON.stringify(reactionPk)
				})
			]);
			
			const freshData = await reactionResponse.json();
			const detailsData = await detailsResponse.json();
			
			// Check for error response
			if (freshData.error) {
				showWorkspaceToast(freshData.error, 'error');
				return;
			}
			
			// Update the reaction in reactions_active with fresh data
			const reactionIndex = reactions_active.findIndex(r => r.pk === reactionPk);
			if (reactionIndex !== -1) {
				// Update the fields in the existing reaction object
				const reaction = reactions_active[reactionIndex];
				
				// Update editable fields from fresh server data
				reaction.fields.description = freshData.description || '';
				reaction.fields.short_name = freshData.short_name || '';
				reaction.fields.confidence_score = freshData.confidence_score || '';
				reaction.fields.substrates_names = JSON.stringify(freshData.substrates_names || []);
				reaction.fields.products_names = JSON.stringify(freshData.products_names || []);
				
				// Update metadata from details endpoint
				reaction.fields.references = detailsData.references || null;
				reaction.fields.ext_links = detailsData.external_links || null;
				reaction.fields.comments = detailsData.comments || null;
				reaction.fields.gene_info = detailsData.gene_info || null;
				
				// Rebuild the card with fresh data
				const card = await buildEditableCard(reaction);
				availableDetailEl.innerHTML = '';
				if (card) availableDetailEl.appendChild(card);
				
			} else {
				showWorkspaceToast('Reaction not found in list', 'error');
			}
		} catch (error) {
			console.error('Error resetting reaction:', error);
			showWorkspaceToast('Failed to reset changes', 'error');
			
			// Try to at least show the current state
			const reaction = reactions_active.find(r => r.pk === reactionPk);
			if (reaction) {
				const card = await buildEditableCard(reaction);
				availableDetailEl.innerHTML = '';
				if (card) availableDetailEl.appendChild(card);
			}
		}
	}

	/**
	 * Save current form changes to the local database
	 * @param {number} reactionPk - The reaction primary key
	 */
	async function saveChangesToLocal(reactionPk) {
		const formState = collectFormState(reactionPk);
		if (!formState) {
			showWorkspaceToast('No form data to save', 'error');
			return;
		}

		// Find the reaction
		const reaction = reactions_active.find(r => r.pk === reactionPk);
		if (!reaction) {
			showWorkspaceToast('Reaction not found', 'error');
			return;
		}

		// Show loading state on button
		const saveBtn = document.getElementById('saveLocalBtn');
		if (saveBtn) {
			saveBtn.disabled = true;
			saveBtn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> <span>Saving...</span>';
		}

		try {
			const response = await fetch('/save_reaction_draft/', {
				method: 'POST',
				headers: {
					'Content-Type': 'application/json',
					'X-CSRFToken': csrfToken,
					'X-Requested-With': 'XMLHttpRequest'
				},
				body: JSON.stringify({
					reactionId: reactionPk,
					description: formState.description,
					abbreviation: formState.abbreviation,
					confidence_score: formState.confidence_score,
					substrates_info: formState.substrates_info,
					products_info: formState.products_info,
					references: formState.references,
					ext_links: formState.ext_links,
					comments: formState.comments,
					gene_info: formState.gene_info
				})
			});

			const data = await response.json();

			if (data.status === 'success') {
				// Update the reaction object in reactions_active
				reaction.fields.description = formState.description;
				reaction.fields.short_name = formState.abbreviation;
				reaction.fields.confidence_score = formState.confidence_score;
				reaction.fields.references = formState.references;
				reaction.fields.ext_links = formState.ext_links;
				reaction.fields.comments = formState.comments;
				reaction.fields.gene_info = formState.gene_info;

				// Invalidate cache so next load gets fresh data
				window.vmhPrepCache.invalidate(reactionPk);

				showWorkspaceToast('Draft saved successfully', 'success');
			} else {
				showWorkspaceToast(data.message || 'Failed to save draft', 'error');
			}
		} catch (error) {
			console.error('Error saving draft:', error);
			showWorkspaceToast('Failed to save draft', 'error');
		} finally {
			// Restore button state
			if (saveBtn) {
				saveBtn.disabled = false;
				saveBtn.innerHTML = '<i class="fas fa-save"></i> <span>Save Draft</span>';
			}
		}
	}

	// Make functions available globally
	window.collectFormState = collectFormState;
	window.updateCacheFromForm = updateCacheFromForm;
	window.resetReactionChanges = resetReactionChanges;
	window.saveChangesToLocal = saveChangesToLocal;

	// New section HTML generator for the redesigned layout
	function createSectionHTMLNew(sectionTitle, className, items, reactionId, isExtLink = false, isRef = false) {
		const iconMap = {
			'References': 'fa-book',
			'External Links': 'fa-external-link-alt',
			'Comments': 'fa-comment-alt'
		};
		const icon = iconMap[sectionTitle] || 'fa-tag';
		
		return `
			<div class="ws-annotation-card ${className}-section">
				<div class="ws-annotation-header">
					<span class="ws-annotation-icon"><i class="fas ${icon}"></i></span>
					<span class="ws-annotation-title">${sectionTitle}</span>
				</div>
				<div class="ws-annotation-body">
					${items.map((item, index) => `
						<div class="${className}-item ws-annotation-item" data-reaction-id="${reactionId}" data-index="${index}">
							${isExtLink ? createExtLinkSelect(item, reactionId, index) : ''}
							${isRef ? createRefSelect(item, reactionId, index) : ''}
							<input type="text" class="${className}-input ws-annotation-input" 
								placeholder="Enter ${sectionTitle.toLowerCase().slice(0, -1)}" 
								value="${item.info}" 
								data-reaction-id="${reactionId}" data-index="${index}">
							<button class="ws-annotation-remove remove-${className}" data-reaction-id="${reactionId}" data-index="${index}">
								<i class="fas fa-trash-alt"></i>
							</button>
						</div>
					`).join('')}
					<button class="ws-annotation-add add-${className}" data-reaction-id="${reactionId}">
						<i class="fas fa-plus"></i>
						<span>Add ${sectionTitle.slice(0, -1)}</span>
					</button>
				</div>
			</div>
		`;
	}

	// 4 – Save & Submit stubs
	availableDetailEl.addEventListener('click', (e) => {
		// Handle Submit to VMH button
		if (e.target.id === 'submitVMHBtn' || e.target.closest('#submitVMHBtn')) {
			const btn = e.target.id === 'submitVMHBtn' ? e.target : e.target.closest('#submitVMHBtn');
			const pk = +btn.dataset.pk;
			// Update cache before submitting
			updateCacheFromForm(pk);
			// reuse addToVMH() – but restrict to this one pk
			addToVMHforSingle(pk);
		}
		
		// Handle Reset Changes button
		if (e.target.id === 'resetChangesBtn' || e.target.closest('#resetChangesBtn')) {
			const btn = e.target.id === 'resetChangesBtn' ? e.target : e.target.closest('#resetChangesBtn');
			const pk = +btn.dataset.pk;
			resetReactionChanges(pk);
		}
		
		// Handle Save Draft button
		if (e.target.id === 'saveLocalBtn' || e.target.closest('#saveLocalBtn')) {
			const btn = e.target.id === 'saveLocalBtn' ? e.target : e.target.closest('#saveLocalBtn');
			const pk = +btn.dataset.pk;
			saveChangesToLocal(pk);
		}
	});

	// Auto-save to cache on form changes (debounced)
	let saveTimeout = null;
	availableDetailEl.addEventListener('input', (e) => {
		// Check if the input is part of the reaction form
		const reactionInput = e.target.closest('[data-reaction-id]');
		if (!reactionInput) return;
		
		const reactionPk = +reactionInput.dataset.reactionId;
		if (!reactionPk) return;

		// Debounce the cache update
		clearTimeout(saveTimeout);
		saveTimeout = setTimeout(() => {
			updateCacheFromForm(reactionPk);
		}, 500);
	});

	// Also save on select changes
	availableDetailEl.addEventListener('change', (e) => {
		const reactionInput = e.target.closest('[data-reaction-id]');
		if (!reactionInput) return;
		
		const reactionPk = +reactionInput.dataset.reactionId;
		if (!reactionPk) return;

		updateCacheFromForm(reactionPk);
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
	// Check both the target and its parent (for clicks on icon/span inside button)
	const addBtn = e.target.closest('.add-reference, .add-ext-link, .add-comment, .ws-annotation-add');
	if (addBtn) {
		const type = addBtn.classList.contains('add-reference') ? 'reference' : 
					 addBtn.classList.contains('add-ext-link') ? 'ext-link' : 'comment';
		const parentSection = addBtn.parentNode;
		const reactionId = addBtn.getAttribute('data-reaction-id');
		const newItem = document.createElement('div');
		newItem.className = `${type}-item ws-annotation-item`;
		const items = parentSection.querySelectorAll(`.${type}-item`);
		const newIndex = items.length; // Calculate new index based on existing items
		newItem.setAttribute('data-reaction-id', reactionId);
		newItem.setAttribute('data-index', newIndex);
		
		// Build the inner HTML with new styling
		let innerHTML = '';
		if (type === 'ext-link') {
			innerHTML += createExtLinkSelect({}, reactionId, newIndex);
		}
		if (type === 'reference') {
			innerHTML += createRefSelect({}, reactionId, newIndex);
		}
		innerHTML += `
			<input type="text" class="${type}-input ws-annotation-input" 
				placeholder="Enter ${type === 'ext-link' ? 'external link' : type}" 
				value="" data-reaction-id="${reactionId}" data-index="${newIndex}">
			<button class="ws-annotation-remove remove-${type}" data-reaction-id="${reactionId}" data-index="${newIndex}">
				<i class="fas fa-trash-alt"></i>
			</button>
		`;
		newItem.innerHTML = innerHTML;
		parentSection.insertBefore(newItem, addBtn);
	}

	// Handle the "Remove" button clicks using event delegation
	// Check both the target and its parent (for clicks on icon inside button)
	const removeBtn = e.target.closest('.remove-reference, .remove-ext-link, .remove-comment, .remove-gene-info, .ws-annotation-remove');
	if (removeBtn) {
		const isGeneInfo = removeBtn.classList.contains('remove-gene-info');
		// Find the parent item - could be .ws-annotation-item, .gene-info-item, .reference-item, etc.
		const parentItem = removeBtn.closest('.ws-annotation-item, .gene-info-item, .reference-item, .ext-link-item, .comment-item');
		if (parentItem) {
			parentItem.remove();
		}
		
		// Update GPR summary if a gene info was removed
		if (isGeneInfo) {
			updateGPRSummary();
		}
	}
});

// Function to update GPR summary based on current gene-info inputs
function updateGPRSummary() {
	const gprDisplayEl = document.querySelector('.ws-gpr-display');
	const gprEmptyEl = document.querySelector('.ws-gpr-empty');
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
		
		if (gprEmptyEl) {
			gprEmptyEl.style.display = 'none';
		}
		
		if (gprDisplayEl) {
			gprDisplayEl.innerHTML = `<span class="ws-gpr-label">GPR Rule:</span><code class="ws-gpr-code">${gprText}</code>`;
			gprDisplayEl.style.display = '';
		} else {
			// Create new display element if it doesn't exist
			const geneSection = document.querySelector('.ws-gene-section');
			if (geneSection && gprEmptyEl) {
				const newDisplay = document.createElement('div');
				newDisplay.className = 'ws-gpr-display';
				newDisplay.innerHTML = `<span class="ws-gpr-label">GPR Rule:</span><code class="ws-gpr-code">${gprText}</code>`;
				geneSection.insertBefore(newDisplay, gprEmptyEl);
			}
		}
	} else {
		if (gprDisplayEl) {
			gprDisplayEl.style.display = 'none';
		}
		if (gprEmptyEl) {
			gprEmptyEl.style.display = '';
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
