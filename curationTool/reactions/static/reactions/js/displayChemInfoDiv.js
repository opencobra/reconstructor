function loadChemInfoDiv(reactionData) {
	const contentDiv = document.querySelector('div.content-div[name="cheminfo-div"]');
	if (!contentDiv) return;

	// Don’t nuke the whole div; leave the placeholder structure in place.
	// Just fill it from reactionData:

	// Header lines + charges
	renderChemInfoHeader({
		massBalanced: !!(reactionData.balanced_count && reactionData.balanced_count[0]),
		chargeBalanced: !!(reactionData.balanced_charge && reactionData.balanced_charge[0]),
		subsAtoms: reactionData.subs_atoms ? reactionData.subs_atoms[0] : {},
		prodsAtoms: reactionData.prods_atoms ? reactionData.prods_atoms[0] : {},
		subsCharge: reactionData.subs_charge ? reactionData.subs_charge[0] : 0,
		prodsCharge: reactionData.prods_charge ? reactionData.prods_charge[0] : 0,
	});

	// Original table (in-place)
	renderAtomComparisonTableInPlace(
		reactionData.subs_atoms ? reactionData.subs_atoms[0] : {},
		reactionData.prods_atoms ? reactionData.prods_atoms[0] : {},
		reactionData.symb_to_name ? reactionData.symb_to_name[0] : {}
	);

	// VMH Formula (make it idempotent and render into a fixed section)
	appendVMHFormulaSection(contentDiv, reactionData);
}

// Helper function to create the balanced text
function createBalancedText(reactionData) {
	const balancedText = document.createElement('p');
	const isBalanced = reactionData.balanced_count[0] && reactionData.balanced_charge[0];
	balancedText.textContent = isBalanced ? 'Balanced' : 'Not Balanced';
	balancedText.className = `balanced-text ${isBalanced ? 'balanced' : 'not-balanced'}`;
	return balancedText;
}

// Helper function to create molecular formula
function createMolecularFormula(reactionData) {
	const formulaParagraph = document.createElement('p');
	formulaParagraph.textContent = 'Molecular Formula: ';
	const atomColorMap = generateAtomColorMap(reactionData);
	const formulaSpan = document.createElement('span');
	formulaSpan.innerHTML = reactionData.molc_formula[0].replace(/([A-Z][a-z]?)/g, (match) => {
		return atomColorMap[match] ? `<span style="color: ${atomColorMap[match]}">${match}</span>` : match;
	});
	formulaParagraph.appendChild(formulaSpan);
	return formulaParagraph;
}

// Helper function to create generic information divs
function createInfoDiv(text, className) {
	const infoDiv = document.createElement('div');
	infoDiv.className = className;
	infoDiv.textContent = text;
	return infoDiv;
}

// Helper function to generate atom color map
function generateAtomColorMap(reactionData) {
	const atomTypes = new Set([...Object.keys(reactionData.subs_atoms[0]), ...Object.keys(reactionData.prods_atoms[0])]);
	const colors = getDistinctColors(atomTypes.size);
	const atomColorMap = {};
	let i = 0;
	atomTypes.forEach((atom) => {
		atomColorMap[atom] = colors[i++];
	});
	return atomColorMap;
}

function getDistinctColors(count) {
	// Function to generate distinct colors for atom types
	const colors = [];
	for (let i = 0; i < count; i++) {
		const hue = (i * 137.508) % 360; // Use golden angle approximation
		colors.push(`hsl(${hue}, 100%, 50%)`);
	}
	return colors;
}

function appendVMHFormulaSection(contentDiv, reactionData) {
	// Reuse (or create) a single container so we don't duplicate the section
	let formulaSection = contentDiv.querySelector('#vmh-formula-section');
	if (!formulaSection) {
		formulaSection = document.createElement('div');
		formulaSection.id = 'vmh-formula-section';
		contentDiv.appendChild(formulaSection);
	}
	// Clear and rebuild
	formulaSection.innerHTML = '';

	// Card wrapper (just for styling)
	const card = document.createElement('div');
	card.className = 'vmh-formula-section';
	card.style.marginTop = '30px';
	card.style.padding = '20px';
	card.style.border = '1px solid #ccc';
	card.style.borderRadius = '8px';
	card.style.backgroundColor = '#f4f4f4';
	card.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.1)';
	card.style.textAlign = 'center';

	const formulaTitle = document.createElement('h3');
	formulaTitle.textContent = 'VMH Formula';
	formulaTitle.style.marginBottom = '20px';
	card.appendChild(formulaTitle);

	// Formula display
	const formulaDisplay = document.createElement('div');
	formulaDisplay.className = 'vmh-formula-display';
	formulaDisplay.style.fontSize = '18px';
	formulaDisplay.style.fontFamily = 'monospace';
	formulaDisplay.style.padding = '10px';
	formulaDisplay.style.border = '1px dashed #ccc';
	formulaDisplay.style.borderRadius = '4px';
	formulaDisplay.style.backgroundColor = '#fff';
	formulaDisplay.style.cursor = 'pointer';

	// Use saved formula if present; otherwise build from arrays
	if (reactionData.rxn_formula) {
		formulaDisplay.textContent = reactionData.rxn_formula;
	} else {
		const formulaText = [];
		// Substrates
		(reactionData.subs_sch || []).forEach((stoich, index) => {
			if (stoich > 1) formulaText.push(parseFloat(stoich).toFixed(1) + ' ');
			const found = reactionData.subs_found ? reactionData.subs_found[index] : false;
			const comp = (reactionData.subs_comps && reactionData.subs_comps[index]) || '';
			const token = found ? `${reactionData.substrates[index]}[${comp}]` : `${reactionData.substrates_names[index]}[${comp}]`;
			formulaText.push(token);
			if (index < reactionData.subs_sch.length - 1) formulaText.push(' + ');
		});

		// Direction
		const dir = reactionData.direction;
		formulaText.push(dir === 'forward' ? ' -> ' : dir === 'bidirectional' ? ' <=> ' : ' ');

		// Products
		(reactionData.prod_sch || []).forEach((stoich, index) => {
			if (stoich > 1) formulaText.push(parseFloat(stoich).toFixed(1) + ' ');
			const found = reactionData.prod_found ? reactionData.prod_found[index] : false;
			const comp = (reactionData.prods_comps && reactionData.prods_comps[index]) || '';
			const token = found ? `${reactionData.products[index]}[${comp}]` : `${reactionData.products_names[index]}[${comp}]`;
			formulaText.push(token);
			if (index < reactionData.prod_sch.length - 1) formulaText.push(' + ');
		});

		formulaDisplay.textContent = formulaText.join('');
	}

	// Click-to-copy
	formulaDisplay.addEventListener('click', () => {
		const text = formulaDisplay.textContent;
		if (navigator.clipboard && navigator.clipboard.writeText) {
			navigator.clipboard
				.writeText(text)
				.then(() => {
					if (typeof showToast === 'function') showToast('Formula copied to clipboard!');
				})
				.catch((err) => console.error('Could not copy text: ', err));
		} else {
			// Fallback
			const textarea = document.createElement('textarea');
			textarea.value = text;
			document.body.appendChild(textarea);
			textarea.select();
			try {
				document.execCommand('copy');
				if (typeof showToast === 'function') showToast('Formula copied to clipboard!');
			} catch (err) {
				console.error('Could not copy text: ', err);
			}
			document.body.removeChild(textarea);
		}
	});

	card.appendChild(formulaDisplay);

	// Generate Abbreviations button
	const generateAbbrButton = document.createElement('button');
	generateAbbrButton.textContent = 'Generate Abbreviations';
	generateAbbrButton.style.marginTop = '20px';
	generateAbbrButton.style.padding = '10px 20px';
	generateAbbrButton.style.fontSize = '16px';
	generateAbbrButton.style.borderRadius = '4px';
	generateAbbrButton.style.border = 'none';
	generateAbbrButton.style.backgroundColor = '#007BFF';
	generateAbbrButton.style.color = '#fff';
	generateAbbrButton.style.cursor = 'pointer';

	generateAbbrButton.addEventListener('click', () => {
		// Opens your existing modal that handles non-VMH abbrev generation & saving
		openAbbreviationModal(reactionData);
	});

	card.appendChild(generateAbbrButton);

	// Attach the card into the stable section container
	formulaSection.appendChild(card);
}

function openAbbreviationModal(reactionData) {
	// Create the modal container
	const modal = document.createElement('div');
	modal.style.position = 'fixed';
	modal.style.top = '50%';
	modal.style.left = '50%';
	modal.style.transform = 'translate(-50%, -50%)';
	modal.style.zIndex = '1000';
	modal.style.backgroundColor = '#fff';
	modal.style.padding = '20px';
	modal.style.borderRadius = '8px';
	modal.style.boxShadow = '0 4px 8px rgba(0, 0, 0, 0.2)';
	modal.style.width = '400px';
	modal.style.textAlign = 'center';

	// Modal title
	const modalTitle = document.createElement('h4');
	modalTitle.textContent = 'Generate Abbreviation for Non-VMH Metabolites';
	modalTitle.style.marginBottom = '20px';
	modal.appendChild(modalTitle);

	// Create the loader element
	const loader = document.createElement('div');
	loader.style.display = 'none';
	loader.style.position = 'absolute';
	loader.style.top = '50%';
	loader.style.left = '50%';
	loader.style.transform = 'translate(-50%, -50%)';
	loader.style.padding = '20px';
	loader.style.backgroundColor = 'white';
	loader.style.color = 'black';
	loader.style.borderRadius = '5px';
	loader.style.zIndex = '1000';
	loader.style.textAlign = 'center';
	loader.style.boxShadow = '0px 10px 15px rgba(0, 0, 0, 0.3)';
	loader.textContent = 'Loading...';
	modal.appendChild(loader);
	// Collect unique metabolites not found in VMH
	const uniqueNotFoundMetabolites = new Map();

	reactionData.substrates_names.forEach((name, index) => {
		if (!reactionData.subs_found[index]) {
			uniqueNotFoundMetabolites.set(name, {
				metabolite: reactionData.substrates[index],
				mtype: reactionData.subs_types[index],
			});
		}
	});

	reactionData.products_names.forEach((name, index) => {
		if (!reactionData.prod_found[index]) {
			// Avoid overwriting if the metabolite is already in the map
			if (!uniqueNotFoundMetabolites.has(name)) {
				uniqueNotFoundMetabolites.set(name, {
					metabolite: reactionData.products[index],
					mtype: reactionData.prods_types[index],
				});
			}
		}
	});
	// Add close button
	const closeButton = document.createElement('button');
	closeButton.textContent = 'Close';
	closeButton.style.marginTop = '20px';
	closeButton.style.padding = '10px 20px';
	closeButton.style.fontSize = '16px';
	closeButton.style.borderRadius = '4px';
	closeButton.style.border = 'none';
	closeButton.style.backgroundColor = '#dc3545';
	closeButton.style.color = '#fff';
	closeButton.style.cursor = 'pointer';

	closeButton.addEventListener('click', () => {
		modal.remove();
	});

	// Generate abbreviations for unique metabolites
	uniqueNotFoundMetabolites.forEach((details, name) => {
		const metaboliteDiv = document.createElement('div');
		metaboliteDiv.style.marginBottom = '10px';

		const nameSpan = document.createElement('span');
		nameSpan.textContent = name;
		metaboliteDiv.appendChild(nameSpan);

		const generateButton = document.createElement('button');
		generateButton.textContent = 'Generate';
		generateButton.style.marginLeft = '10px';
		generateButton.style.padding = '5px 10px';
		generateButton.style.fontSize = '14px';
		generateButton.style.borderRadius = '4px';
		generateButton.style.border = 'none';
		generateButton.style.backgroundColor = '#28a745';
		generateButton.style.color = '#fff';
		generateButton.style.cursor = 'pointer';
		function escapeRegExp(string) {
			return string.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
		}

		generateButton.addEventListener('click', () => {
			showLoaderdiv(generateButton, closeButton); // Pass the closeButton
			fetch('/create-formula-abbr/', {
				method: 'POST',
				headers: {
					'Content-Type': 'application/json',
				},
				body: JSON.stringify({
					metabolite: details.metabolite,
					mtype: details.mtype,
					metabolite_name: name,
				}),
			})
				.then((response) => response.json())
				.then((data) => {
					console.log(data);
					if (data.abbr) {
						hideLoaderdiv(generateButton, closeButton); // Pass the closeButton

						// Update the name span to show the abbreviation
						nameSpan.textContent = data.abbr;

						const formulaDisplay = document.querySelector('.vmh-formula-display');

						if (formulaDisplay) {
							const escapedName = escapeRegExp(name);
							const regex = new RegExp(`\\b${escapedName}\\[(c|g|e|l|m|n|x|r)\\]`, 'gi');

							// Replace the matched metabolites with the abbreviation and the same compartment
							const updatedFormula = formulaDisplay.textContent.replace(regex, `${data.abbr}[$1]`);

							formulaDisplay.textContent = updatedFormula;
						}
						// Save the updated formula to the backend
						const updatedFormulaText = document.querySelector('.vmh-formula-display').textContent;
						const reactionId = reactionData.reaction_id; // Assume the reaction ID is available
						fetch('/save-formula/', {
							method: 'POST',
							headers: {
								'Content-Type': 'application/json',
								'X-CSRFToken': csrfToken, // Replace csrfToken with the actual CSRF token variable
							},
							body: JSON.stringify({
								formula: updatedFormulaText,
								reaction_id: reactionId,
							}),
						})
							.then((saveResponse) => {
								if (saveResponse.ok) {
									console.log('Formula successfully saved.');
								} else {
									console.error('Failed to save the formula.');
								}
							})
							.catch((saveError) => console.error('Error saving formula:', saveError));

						// Remove the generate button
						generateButton.remove();
					}
				})
				.catch((error) => {
					console.error('Error:', error);
					hideLoaderdiv(generateButton, closeButton);
				});
		});

		metaboliteDiv.appendChild(generateButton);
		modal.appendChild(metaboliteDiv);
	});

	modal.appendChild(closeButton);

	// Append modal to document body
	document.body.appendChild(modal);
}

function populateAtomComparisonTable(substratesData, productsData, listDiv, charge, symb_to_name) {
	listDiv.style.display = 'flex';
	listDiv.style.justifyContent = 'space-between';
	listDiv.style.alignItems = 'flex-start';
	listDiv.style.padding = '20px';
	listDiv.style.border = '1px solid #ddd';
	listDiv.style.borderRadius = '8px';
	listDiv.style.marginTop = '5%';
	listDiv.style.backgroundColor = '#f9f9f9';
	listDiv.style.fontSize = '16px'; // Replace '14px' with your desired font size

	const allAtoms = { ...substratesData, ...productsData };

	const table = document.createElement('table');
	table.className = 'atom-comparison-table';
	table.style.width = '70%';
	table.style.borderCollapse = 'collapse';

	const thead = document.createElement('thead');
	const headerRow = document.createElement('tr');
	['Atom', 'Substrates', 'Products'].forEach((text) => {
		const th = document.createElement('th');
		th.textContent = text;
		th.style.borderBottom = '2px solid #000';
		th.style.padding = '10px';
		th.style.backgroundColor = '#eaeaea';
		th.style.textAlign = 'center';
		th.style.fontWeight = 'bold';
		th.style.fontSize = '19px'; // Replace '16px' with your desired font size
		headerRow.appendChild(th);
	});
	thead.appendChild(headerRow);
	table.appendChild(thead);

	const tbody = document.createElement('tbody');
	Object.keys(allAtoms).forEach((atom) => {
		const row = document.createElement('tr');
		const substratesCount = substratesData[atom] || 0;
		const productsCount = productsData[atom] || 0;

		if (substratesCount !== productsCount) {
			row.style.backgroundColor = '#ffcccc'; // Highlight in red
		}

		[symb_to_name[atom] || atom, substratesCount, productsCount].forEach((text, index) => {
			const cell = document.createElement('td');
			cell.textContent = text;
			cell.style.padding = '8px';
			cell.style.borderBottom = '1px solid #ddd';
			cell.style.textAlign = 'center';
			if (index === 0) {
				cell.style.fontWeight = 'bold';
			}
			row.appendChild(cell);
		});
		tbody.appendChild(row);
	});
	table.appendChild(tbody);

	listDiv.appendChild(table);

	const chargeDiv = document.createElement('div');
	chargeDiv.innerHTML = `<strong>Charge Info</strong><br>Substrates: ${charge.substrates}<br>Products: ${charge.products}`;
	chargeDiv.style.marginLeft = '20px';
	chargeDiv.style.padding = '10px';
	chargeDiv.style.border = '1px solid #ccc';
	chargeDiv.style.borderRadius = '8px';
	chargeDiv.style.backgroundColor = '#fff';
	chargeDiv.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
	listDiv.appendChild(chargeDiv);
}
function showLoaderdiv(button, closeButton) {
	button.disabled = true; // Disable the button
	button.style.backgroundColor = '#ccc';
	button.style.cursor = 'not-allowed';
	button.textContent = 'Generating...';
	// Add a spinner and loading message over the button
	const spinner = document.createElement('span');
	spinner.classList.add('spinner');
	spinner.style.marginLeft = '10px';
	spinner.style.display = 'inline-block'; // Ensure it is inline with the button text
	spinner.style.border = '3px solid rgba(0, 0, 0, 0.1)'; // Light gray border
	spinner.style.borderTop = '3px solid #3498db'; // Contrasting blue for spinning effect
	spinner.style.borderRadius = '50%';
	spinner.style.width = '16px';
	spinner.style.height = '16px';
	spinner.style.animation = 'spin 1s linear infinite';
	button.appendChild(spinner);

	// Disable close button
	closeButton.disabled = true;
	closeButton.style.backgroundColor = '#aaa';
	closeButton.style.cursor = 'not-allowed';
}

function hideLoaderdiv(button, closeButton) {
	button.disabled = false; // Re-enable the button
	button.style.backgroundColor = '#28a745'; // Restore original color
	button.style.cursor = 'pointer'; // Restore cursor

	// Remove the spinner
	const spinner = button.querySelector('.spinner');
	if (spinner) {
		spinner.remove();
	}

	// Re-enable close button
	closeButton.disabled = false;
	closeButton.style.backgroundColor = '#dc3545';
	closeButton.style.cursor = 'pointer';
}

// Aggregate totals from the current form (same logic you already use)
function getTotalsFromForm() {
	const totalAtomsSubs = {};
	const totalAtomsProds = {};
	let totalChargeSubs = 0;
	let totalChargeProds = 0;

	const substratesDiv = document.getElementById('substratesDiv');
	const productsDiv = document.getElementById('productsDiv');
	const substrateGroups = substratesDiv ? substratesDiv.querySelectorAll('.inputs-group') : [];
	const productGroups = productsDiv ? productsDiv.querySelectorAll('.inputs-group') : [];

	substrateGroups.forEach((group) => {
		if (group.dataset.atomCounts) {
			const counts = JSON.parse(group.dataset.atomCounts);
			for (const elem in counts) {
				totalAtomsSubs[elem] = (totalAtomsSubs[elem] || 0) + counts[elem];
			}
			totalChargeSubs += parseFloat(group.dataset.charge) || 0;
		}
	});

	productGroups.forEach((group) => {
		if (group.dataset.atomCounts) {
			const counts = JSON.parse(group.dataset.atomCounts);
			for (const elem in counts) {
				totalAtomsProds[elem] = (totalAtomsProds[elem] || 0) + counts[elem];
			}
			totalChargeProds += parseFloat(group.dataset.charge) || 0;
		}
	});

	return {
		subsAtoms: totalAtomsSubs,
		prodsAtoms: totalAtomsProds,
		subsCharge: totalChargeSubs,
		prodsCharge: totalChargeProds,
	};
}

// Render the live comparison into cheminfo-div using your existing table builder
function renderLiveChemInfoFromForm() {
	const holder = ensureChemLivePanel();
	if (!holder) return;

	// Clear previous
	holder.innerHTML = '';

	const { subsAtoms, prodsAtoms, subsCharge, prodsCharge } = getTotalsFromForm();

	// Balanced badges
	const allAtoms = new Set([...Object.keys(subsAtoms), ...Object.keys(prodsAtoms)]);
	let massBalanced = true;
	for (const a of allAtoms) {
		if ((subsAtoms[a] || 0) !== (prodsAtoms[a] || 0)) {
			massBalanced = false;
			break;
		}
	}
	const chargeBalanced = subsCharge === prodsCharge;

	const badge = document.createElement('p');
	badge.className = `balanced-text ${massBalanced && chargeBalanced ? 'balanced' : 'not-balanced'}`;
	badge.textContent = massBalanced && chargeBalanced ? 'Balanced' : 'Not Balanced';
	badge.style.fontWeight = 'bold';
	badge.style.marginBottom = '12px';
	holder.appendChild(badge);

	// Build the same nice table+charge using your existing function
	// Note: we pass an empty symb_to_name map so it falls back to the symbol (your function already does that).
	populateAtomComparisonTable(
		subsAtoms,
		prodsAtoms,
		holder,
		{ substrates: subsCharge, products: prodsCharge },
		{} // symb_to_name
	);
}

// Convert atom counts to a formula string like C6H12O6 (with <sub>n</sub>)
function countsToFormula(counts) {
	if (!counts) return '—';
	const keys = Object.keys(counts);
	if (!keys.length) return '—';

	// Hill order: C, H, then alphabetical
	const order = (a, b) => {
		if (a === 'C') return -1;
		if (b === 'C') return 1;
		if (a === 'H' && b !== 'C') return -1;
		if (b === 'H' && a !== 'C') return 1;
		return a.localeCompare(b);
	};

	return keys
		.sort(order)
		.map((k) => `${k}<sub>${counts[k]}</sub>`)
		.join('');
}

// Fill the existing .atom-comparison-table (tbody) with current counts
function renderAtomComparisonTableInPlace(substratesData, productsData, symb_to_name) {
	const root = document.querySelector('div.content-div[name="cheminfo-div"]');
	if (!root) return;
	const tbody = root.querySelector('.atom-comparison-table tbody');
	if (!tbody) return;

	const allAtoms = new Set([...Object.keys(substratesData || {}), ...Object.keys(productsData || {})]);
	const mapName = symb_to_name || {};

	const frag = document.createDocumentFragment();
	Array.from(allAtoms)
		.sort()
		.forEach((atom) => {
			const row = document.createElement('tr');
			const subs = (substratesData && substratesData[atom]) || 0;
			const prods = (productsData && productsData[atom]) || 0;

			if (subs !== prods) row.style.backgroundColor = '#ffcccc';

			const tdAtom = document.createElement('td');
			tdAtom.textContent = mapName[atom] || atom;
			tdAtom.style.fontWeight = 'bold';

			const tdSubs = document.createElement('td');
			tdSubs.textContent = subs;
			tdSubs.style.textAlign = 'center';

			const tdProds = document.createElement('td');
			tdProds.textContent = prods;
			tdProds.style.textAlign = 'center';

			row.appendChild(tdAtom);
			row.appendChild(tdSubs);
			row.appendChild(tdProds);
			frag.appendChild(row);
		});

	// Replace whole tbody
	tbody.innerHTML = '';
	tbody.appendChild(frag);
}

// Updates Balanced / Mass? / Charge? + Molecular Formula (per-metabolite, not aggregated)
// Reuses a stable color map so element colors don't jump between updates.
function renderChemInfoHeader({ massBalanced, chargeBalanced, subsAtoms, prodsAtoms, subsCharge, prodsCharge }) {
	const root = document.querySelector('div.content-div[name="cheminfo-div"]');
	if (!root) return;

	// Balanced
	const statusSpan = root.querySelector('.balanced-status .status-text');
	if (statusSpan) statusSpan.textContent = massBalanced && chargeBalanced ? 'Yes' : 'No';

	// Mass Balanced?
	const massSpan = root.querySelector('.balance-info.mass span');
	if (massSpan) massSpan.textContent = massBalanced ? 'Yes' : 'No';

	// Charge Balanced?
	const chargeSpan = root.querySelector('.balance-info.charge span');
	if (chargeSpan) chargeSpan.textContent = chargeBalanced ? 'Yes' : 'No';

	// Charge info box
	const chargeSub = root.querySelector('.charge-info .charge-sub');
	const chargeProd = root.querySelector('.charge-info .charge-prod');
	if (chargeSub) chargeSub.textContent = subsCharge ?? 0;
	if (chargeProd) chargeProd.textContent = prodsCharge ?? 0;

	// Molecular Formula (PER metabolite: e.g., H2O + O2 | ...), not aggregated
	const formulaSpan = root.querySelector('.formula-placeholder');
	if (!formulaSpan) return;

	// Gather all element symbols currently present so we can build a stable color map
	const atomsSet = collectAtomsFromGroups('substratesDiv');
	for (const a of collectAtomsFromGroups('productsDiv')) atomsSet.add(a);

	const colorMap = ensureStableAtomColorMap(formulaSpan, atomsSet);
	const subsHTML = buildSideFormula('substratesDiv', colorMap);
	const prodsHTML = buildSideFormula('productsDiv', colorMap);

	formulaSpan.innerHTML = `Subs: ${subsHTML} &nbsp; | &nbsp; Prods: ${prodsHTML}`;
}

// Collect all element symbols present in a side by reading each group's data-atomCounts
function collectAtomsFromGroups(sideDivId) {
	const out = new Set();
	const container = document.getElementById(sideDivId);
	if (!container) return out;
	container.querySelectorAll('.inputs-group').forEach((group) => {
		if (!group.dataset.atomCounts) return;
		try {
			const counts = JSON.parse(group.dataset.atomCounts);
			Object.keys(counts || {}).forEach((k) => {
				if (counts[k] > 0) out.add(k);
			});
		} catch (_) {}
	});
	return out;
}

// Build the per-metabolite formula string (e.g., "2 H2O + O2") from groups in a side.
function buildSideFormula(sideDivId, colorMap) {
	const container = document.getElementById(sideDivId);
	if (!container) return '—';

	const parts = [];
	container.querySelectorAll('.inputs-group').forEach((group) => {
		if (!group.dataset.atomCounts) return;
		let counts;
		try {
			counts = JSON.parse(group.dataset.atomCounts);
		} catch {
			return;
		}

		// Read the stoichiometry input inside this group (your markup uses duplicate ids, so scope to group)
		const stoichInput = group.querySelector('input[type="number"]');
		const coeff = stoichInput ? parseFloat(stoichInput.value) : 1;
		const coeffTxt = coeff && coeff !== 1 ? `${stripTrailingZeros(coeff)} ` : '';

		const formulaHTML = countsToFormulaHTML(counts, colorMap);
		parts.push(`${coeffTxt}${formulaHTML}`);
	});

	return parts.length ? parts.join(' + ') : '—';
}

// Turn an element-count object into Hill-ordered, colored HTML with subscripts (e.g., C6H12O6)
function countsToFormulaHTML(counts, colorMap) {
	if (!counts) return '—';
	const keys = Object.keys(counts).filter((k) => counts[k] > 0);
	if (!keys.length) return '—';

	keys.sort(hillOrder);
	return keys
		.map((k) => {
			const sym = colorMap && colorMap[k] ? `<span style="color:${colorMap[k]}">${k}</span>` : k;
			const n = counts[k];
			return n === 1 ? sym : `${sym}<sub>${n}</sub>`;
		})
		.join('');
}

// Hill-system ordering: C, H, then alphabetical
function hillOrder(a, b) {
	if (a === 'C') return -1;
	if (b === 'C') return 1;
	if (a === 'H' && b !== 'C') return -1;
	if (b === 'H' && a !== 'C') return 1;
	return a.localeCompare(b);
}

// Make a stable color map for current atoms; reuse previous one if stored on the span
function ensureStableAtomColorMap(formulaSpan, atomsSet) {
	// Reuse prior map if present
	try {
		if (formulaSpan.dataset.atomColorMap) {
			const prev = JSON.parse(formulaSpan.dataset.atomColorMap);
			// Ensure it covers any new atoms (assign colors deterministically if needed)
			const missing = [...atomsSet].filter((a) => !prev[a]);
			if (missing.length === 0) return prev;

			const addColors =
				typeof getDistinctColors === 'function' ? getDistinctColors(missing.length) : missing.map((_, i) => `hsl(${(i * 137.508) % 360} 100% 50%)`);

			missing.forEach((el, i) => {
				prev[el] = addColors[i];
			});
			formulaSpan.dataset.atomColorMap = JSON.stringify(prev);
			return prev;
		}
	} catch (_) {}

	// Build fresh
	const atoms = [...atomsSet].sort(hillOrder);
	const colors =
		typeof getDistinctColors === 'function' ? getDistinctColors(atoms.length) : atoms.map((_, i) => `hsl(${(i * 137.508) % 360} 100% 50%)`);

	const map = {};
	atoms.forEach((el, i) => {
		map[el] = colors[i];
	});
	try {
		formulaSpan.dataset.atomColorMap = JSON.stringify(map);
	} catch (_) {}
	return map;
}

function stripTrailingZeros(n) {
	const s = String(n);
	return s.includes('.') ? s.replace(/\.?0+$/, '') : s;
}
