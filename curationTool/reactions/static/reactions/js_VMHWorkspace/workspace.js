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

	// Load available reactions
	reactions_active.forEach((r) => {
		const li = document.createElement('li');
		li.className = 'item';
		li.dataset.pk = r.pk;
		li.textContent = r.fields.short_name;
		availableListEl.appendChild(li);
	});

	// 2 – on click → render editable form
	availableListEl.addEventListener('click', async (e) => {
		// ← mark as async
		if (!e.target.matches('.item')) return;
		const pk = +e.target.dataset.pk;
		const rxn = reactions_active.find((r) => r.pk === pk);
		if (!rxn) return;

		availableDetailEl.innerHTML = ''; // clear

		const card = await buildEditableCard(rxn); // ← await it
		if (card) availableDetailEl.appendChild(card);

		document.querySelectorAll('#wsAvailableReactionList .item').forEach((el) => el.classList.remove('active'));
		e.target.classList.add('active');
	});

	// 3 – helper that re-uses your modal-builders
	async function buildEditableCard(reaction) {
		let vmhResponse;

		const reactionIndex = reaction.pk; // ← Add this!
		const reactionId = reaction.pk; // ← Add this!

		document.getElementById('loadingIndicator').style.display = 'flex';
		document.getElementById('loadingText').textContent = 'Gathering Data and (if needed) Generating Abbreviations for Selected Reactions';

		try {
			vmhResponse = await callPrepareAddToVMH([reaction.pk]);
		} catch (error) {
			console.error('Error fetching VMH preparation data: ', error);
			document.getElementById('alertMessage').textContent =
				'An error occurred while fetching data from VMH. Check console for error. Please try again later.';
			document.getElementById('alertModal').style.display = 'block';
			document.getElementById('loadingIndicator').style.display = 'none';
			return '';
		}

		if (vmhResponse.status === 'error') {
			document.getElementById('alertMessage').textContent = vmhResponse.message;
			document.getElementById('alertModal').style.display = 'block';
			document.getElementById('loadingIndicator').style.display = 'none';
			return '';
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
                        <input type="text" class="reaction-name-input" placeholder="Reaction Description"
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
			if (item.info) item.info = item.info.split(';')[0];
			return item;
		});

		listItem.innerHTML += createSectionHTML('References', 'reference', references, reaction.pk, false, true);
		listItem.innerHTML += createSectionHTML('External Links', 'ext-link', ext_links, reaction.pk, true, false);
		listItem.innerHTML += createSectionHTML('Comments', 'comment', comments, reaction.pk);
		listItem.innerHTML += createSectionHTML('Gene Info', 'gene-info', gene_info, reaction.pk, false, false, true);

		listItem.innerHTML += `<div class="ws-rxn-actions">
        <button class="ui primary button" id="submitVMHBtn" data-pk="${reaction.pk}">Add to VMH</button>
    </div>`;

		return listItem;
	}

	// 4 – Save & Submit stubs  (wire into your existing endpoints)
	availableDetailEl.addEventListener('click', (e) => {
		if (e.target.id === 'submitVMHBtn') {
			const pk = +e.target.dataset.pk;
			// reuse addToVMH() – but restrict to this one pk
			addToVMHforSingle(pk);
		}
	});

	function addToVMHforSingle(pk) {
		// thin wrapper around your big addToVMH(); adapt as you wish
		window.checkedReactions = [String(pk)];
		addToVMH(); // comes from handleAdd2VMH.js
	}
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
		e.target.parentElement.remove();
	}
});

const parsedAddedReactions = reactions_added;

const tableBody = document.getElementById('addedReactionTableBody');
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
