document.getElementById('subsystemField').addEventListener('keyup', function (event) {
	const inputVal = this.value.toLowerCase();
	const dropdown = document.getElementById('subsystemDropdown');
	dropdown.innerHTML = ''; // Clear current dropdown content
	const matches = subsystemList.filter((subsystem) => subsystem.toLowerCase().includes(inputVal));

	matches.forEach((match) => {
		const element = document.createElement('div');
		element.textContent = match;
		element.addEventListener('click', function () {
			document.getElementById('subsystemField').value = this.textContent; // Fill input field on click
			dropdown.style.display = 'none'; // Hide the dropdown after selection
		});
		dropdown.appendChild(element);
	});

	dropdown.style.display = matches.length > 0 ? 'block' : 'none';
});
document.addEventListener('click', function (event) {
	const dropdown = document.getElementById('subsystemDropdown');
	const inputField = document.getElementById('subsystemField');

	// Check if the click is outside the dropdown and input field
	if (!dropdown.contains(event.target) && !inputField.contains(event.target)) {
		dropdown.style.display = 'none';
	}
});
function updatePlaceholder(selectElement, inputElement) {
	selectElement.addEventListener('change', function () {
		var placeholderText = '';
		switch (this.value) {
			case 'VMH':
				placeholderText = 'VMH Abbreviation';
				break;
			case 'ChEBI ID':
				placeholderText = 'ChEBI ID';
				break;
			case 'ChEBI Name':
				placeholderText = 'ChEBI Name';
				break;
			case 'SwissLipids':
				placeholderText = 'SwissLipids ID';
				break;
			case 'PubChem ID':
				placeholderText = 'PubChem ID';
				break;
			case 'MDL Mol file':
				break;
			case 'Saved':
				placeholderText = 'Search...';
				break;
			case 'Draw':
				placeholderText = ''; // No placeholder for these types
				break;
		}
		console.log('Setting placeholder to:', placeholderText);
		inputElement.placeholder = placeholderText;
	});

	// Immediately update placeholder for the current selection
	selectElement.dispatchEvent(new Event('change'));
}
document.addEventListener('DOMContentLoaded', function () {
	// const reactionInput = document.getElementById('reactionAbbreviation');
	// reactionInput.placeholder = 'Select a source first';
	// reactionInput.disabled = true;

	// Attach the placeholder update logic to the initial fields
	var initialSubstrateTypeSelect = document.querySelector('select[name="substrates_type"]');
	var initialSubstrateInput = document.querySelector('input[name="substrates"]');
	updatePlaceholder(initialSubstrateTypeSelect, initialSubstrateInput);

	var initialProductTypeSelect = document.querySelector('select[name="products_type"]');
	var initialProductInput = document.querySelector('input[name="products"]');
	updatePlaceholder(initialProductTypeSelect, initialProductInput);

	attachEventListenersToDoneButtons();
	attachEventListenersToDoneAllButtons();
});
// Adds a new substrate field to the reaction form when the 'Add Substrate' button is clicked.
document.getElementById('addSubstrate').addEventListener('click', function () {
	addField('substratesDiv', 'substrates', 'subs_sch');
});
// Adds a new product field to the reaction form when the 'Add Product' button is clicked.
document.getElementById('addProduct').addEventListener('click', function () {
	addField('productsDiv', 'products', 'prod_sch');
});
// Sets up event listeners for all 'Remove' buttons to delete their parent input group.
var removeButtons = document.querySelectorAll('.remove-field-btn');
removeButtons.forEach(function (button) {
	button.onclick = function () {
		removeField(this);
	};
});

document.getElementById('applyAllSubsComps').addEventListener('click', function () {
	var selectedValue = document.getElementById('subsCompsApplyAllSelect').value;
	// Apply the selected compartment to all substrate compartment selects
	var substrateCompsSelects = document.querySelectorAll('#substratesDiv .inputs-group select[name="subs_comps"]');
	substrateCompsSelects.forEach(function (select) {
		if (!select.disabled) {
			select.value = selectedValue;
		}
	});
});

document.getElementById('applyAllSubsType').addEventListener('click', function () {
	var selectedValue = document.getElementById('subsTypeApplyAllSelect').value;
	var substrateTypeSelects = document.querySelectorAll('#substratesDiv .inputs-group select[name="substrates_type"]');

	substrateTypeSelects.forEach(function (select, index) {
		try {
			if (!select.disabled) {
				select.value = selectedValue;
				toggleFileInput(select.closest('.inputs-group'), selectedValue);
				handleSelectChange({ target: select });
				const identifierInput = select.closest('.inputs-group')?.querySelector('.cell-identifier input[type="text"]');
				if (identifierInput) {
					updatePlaceholder(select, identifierInput);
				}
			}
		} catch (error) {
			console.error('Error processing select element at index:', index, error);
		}
	});
});

document.getElementById('applyAllProdsComps').addEventListener('click', function () {
	var selectedValue = document.getElementById('prodsCompsApplyAllSelect').value;
	var productCompsSelects = document.querySelectorAll('#productsDiv .inputs-group select[name="prod_comps"]');
	productCompsSelects.forEach(function (select) {
		if (!select.disabled) {
			select.value = selectedValue;
		}
	});
});
document.getElementById('applyAllProdsType').addEventListener('click', function () {
	var selectedValue = document.getElementById('prodsTypeApplyAllSelect').value;
	var productTypeSelects = document.querySelectorAll('#productsDiv .inputs-group select[name="products_type"]');

	// Logging the NodeList to verify its content

	productTypeSelects.forEach(function (select, index) {
		try {
			if (!select.disabled) {
				select.value = selectedValue;
				toggleFileInput(select.closest('.inputs-group'), selectedValue);
				handleSelectChange({ target: select });
				const identifierInput = select.closest('.inputs-group')?.querySelector('.cell-identifier input[type="text"]');
				if (identifierInput) {
					updatePlaceholder(select, identifierInput);
				}
			}
		} catch (error) {
			console.error('Error processing select element at index:', index, error);
		}
		mi;
	});
});

function removeField(button) {
	let inputsGroup = button.closest('.inputs-group');
	inputsGroup.remove();
	updateAtomChargeCounters();
}

function toggleRotation(event) {
	const iconElement = event.currentTarget.querySelector('img, i');
	if (iconElement) {
		iconElement.classList.toggle('rotated');
	}
}
document.getElementById('toggleApplyAllsubs').addEventListener('click', function () {
	var applyOptions = document.getElementById('applyAllOptionsSubs');
	applyOptions.style.display = applyOptions.style.display === 'none' ? 'grid' : 'none';
});
document.getElementById('toggleApplyAllprods').addEventListener('click', function () {
	var applyOptions = document.getElementById('applyAllOptionsProds');
	applyOptions.style.display = applyOptions.style.display === 'none' ? 'grid' : 'none';
});

document.getElementById('toggleApplyAllsubs').addEventListener('click', toggleRotation);

document.getElementById('toggleApplyAllprods').addEventListener('click', toggleRotation);

// Initializes event listeners on existing select elements for file input toggling upon DOM content load.
document.addEventListener('DOMContentLoaded', function () {
	// Select both 'products_type' and 'substrates_type' select elements
	document.querySelectorAll('select[name="products_type"], select[name="substrates_type"]').forEach((selectElement) => {
		// Initial call to handle the current state
		toggleFileInput(selectElement.closest('.inputs-group'), selectElement.value);
		handleMetaboliteTypeChange(selectElement);

		// Attach change event listener to these select elements
		selectElement.addEventListener('change', function () {
			toggleFileInput(selectElement.closest('.inputs-group'), selectElement.value);
			handleMetaboliteTypeChange(selectElement);
		});
	});
});

function buildReactantRow(containerId, inputName, numberName) {
	const row = document.createElement('div');
	row.className = 'inputs-group reactant-grid';

	const removeCell = document.createElement('div');
	removeCell.className = 'cell cell-remove';
	const removeBtn = document.createElement('button');
	removeBtn.type = 'button';
	removeBtn.className = 'ui inverted red button remove-field-btn';
	removeBtn.setAttribute('aria-label', containerId === 'substratesDiv' ? 'Remove substrate' : 'Remove product');
	removeBtn.innerHTML = '<span aria-hidden="true">×</span>';
	removeBtn.onclick = function () {
		removeField(this);
	};
	removeCell.appendChild(removeBtn);
	row.appendChild(removeCell);

	const stoichCell = document.createElement('div');
	stoichCell.className = 'cell cell-stoich';
	const numberInput = document.createElement('input');
	numberInput.type = 'number';
	numberInput.name = numberName;
	numberInput.id = numberName;
	numberInput.min = '1';
	numberInput.value = '1';
	stoichCell.appendChild(numberInput);
	row.appendChild(stoichCell);

	const identifierCell = document.createElement('div');
	identifierCell.className = 'cell cell-identifier';
	const mainInput = document.createElement('input');
	mainInput.type = 'text';
	mainInput.name = inputName;
	mainInput.id = inputName;
	identifierCell.appendChild(mainInput);
	row.appendChild(identifierCell);

	const compartmentCell = document.createElement('div');
	compartmentCell.className = 'cell cell-compartment';
	const compartmentSelect = createCompartmentSelect(containerId === 'substratesDiv' ? 'subs_comps' : 'prod_comps');
	compartmentCell.appendChild(compartmentSelect);
	row.appendChild(compartmentCell);

	const typeCell = document.createElement('div');
	typeCell.className = 'cell cell-type';
	const selectInputType = document.createElement('select');
	selectInputType.name = `${inputName}_type`;
	selectInputType.id = `${inputName}_type`;
	['VMH', 'ChEBI ID', 'ChEBI Name', 'SwissLipids', 'MDL Mol file', 'Draw', 'PubChem ID', 'Saved'].forEach(function (type) {
		const option = document.createElement('option');
		option.value = type;
		option.text = type === 'Saved' ? 'My Metabolites' : type;
		selectInputType.appendChild(option);
	});
	typeCell.appendChild(selectInputType);
	row.appendChild(typeCell);

	const verifyCell = document.createElement('div');
	verifyCell.className = 'cell cell-verify';
	const verifyField = document.createElement('div');
	verifyField.className = 'verify-field';
	const statusDot = document.createElement('span');
	statusDot.className = 'status-dot';
	statusDot.style.display = 'none';
	const doneButton = document.createElement('button');
	doneButton.type = 'button';
	doneButton.className = 'done-field-btn';
	doneButton.textContent = 'Verify';
	const hiddenInput = document.createElement('input');
	hiddenInput.type = 'text';
	hiddenInput.className = `${inputName}-name`;
	hiddenInput.name = `${inputName}_name`;
	hiddenInput.style.display = 'none';
	verifyField.appendChild(statusDot);
	verifyField.appendChild(doneButton);
	verifyField.appendChild(hiddenInput);
	verifyCell.appendChild(verifyField);
	row.appendChild(verifyCell);

	selectInputType.addEventListener('change', function () {
		toggleFileInput(row, this.value);
		handleMetaboliteTypeChange(this);
	});

	return {
		row,
		mainInput,
		numberInput,
		compartmentSelect,
		selectInputType,
		statusDot,
		doneButton,
		nameField: hiddenInput
	};
}

// Dynamically adds a new field container with specified inputs and select options to the given container ID.
function addField(containerId, inputName, numberName) {
	const { row, mainInput, selectInputType } = buildReactantRow(containerId, inputName, numberName);
	document.getElementById(containerId).appendChild(row);
	updatePlaceholder(selectInputType, mainInput);
	toggleFileInput(row, selectInputType.value);
	attachEventListenersToSelects();
	attachEventListenersToDoneButtons();
}

// Adds a field to the specified container with pre-filled data for substrates or products.
function addFieldWithData(container, name, schName, value, schValue, compValue, type, metab_name = null, metab_id = null) {
	return new Promise(async (resolve, reject) => {
		const { row, mainInput, numberInput, compartmentSelect, selectInputType, nameField } = buildReactantRow(container.id, name, schName);
		numberInput.value = schValue;
		mainInput.value = value;
		compartmentSelect.value = compValue;
		selectInputType.value = type;
		if (metab_name) {
			nameField.value = metab_name;
		}

		container.appendChild(row);
		updatePlaceholder(selectInputType, mainInput);
		toggleFileInput(row, selectInputType.value);
		handleMetaboliteTypeChange(selectInputType);

		if (type === 'MDL Mol file') {
			const fileInput = row.querySelector('.cell-identifier input[type="file"]');
			try {
				await loadFileToInputField(value, fileInput);
				resolve();
			} catch (error) {
				reject(error);
			}
		} else {
			resolve();
		}

		if (type === 'Saved') {
			const autocompleteContainer = row.querySelector('.autocomplete-container');
			if (autocompleteContainer) {
				const searchInput = autocompleteContainer.querySelector('.autocomplete-input');
				const hiddenInputField = autocompleteContainer.querySelector('input[type="hidden"]');
				if (searchInput) {
					searchInput.value = metab_name;
				}
				if (hiddenInputField) {
					hiddenInputField.value = value;
				}
			}
		}

		attachEventListenersToSelects();
		attachEventListenersToDoneButtons();
	});
}

function loadFileToInputField(url, fileInput) {
	// Return a new promise that resolves when the file has been loaded
	return new Promise((resolve, reject) => {
		fetchFileAsBlob(
			url,
			(blob) => {
				const fileName = url.split('/').pop();
				const file = new File([blob], fileName, { type: 'chemical/x-mdl-molfile', lastModified: new Date().getTime() });
				const dataTransfer = new DataTransfer();
				dataTransfer.items.add(file);
				fileInput.files = dataTransfer.files;
				resolve(); // Resolve the promise once the file is loaded
			},
			reject
		); // Reject the promise on error
	});
}

// Fetches a file as a blob from a URL
function fetchFileAsBlob(url, callback, errorCallback) {
	var xhr = new XMLHttpRequest();
	xhr.onload = function () {
		if (xhr.status === 200) {
			callback(xhr.response);
		} else {
			errorCallback(new Error('Failed to fetch blob'));
		}
	};
	xhr.onerror = errorCallback;
	xhr.open('GET', url);
	xhr.responseType = 'blob';
	xhr.send();
}

function createCompartmentSelect(prod_or_subs) {
	var compartmentSelect = document.createElement('select');
	compartmentSelect.id = 'compartment';
	compartmentSelect.name = prod_or_subs;
	['c', 'e', 'g', 'm', 'l', 'n', 'r', 'x', '-'].forEach(function (compartment) {
		var option = document.createElement('option');
		option.value = compartment;
		option.text = compartment;
		compartmentSelect.appendChild(option);
	});
	compartmentSelect.value = 'c'; // Set default value to 'c'
	return compartmentSelect;
}
// Updates the form fields with given substrates and products data.
async function updateFormFields(data) {
	let substratesDiv = document.getElementById('substratesDiv');
	let productsDiv = document.getElementById('productsDiv');

	function clearInputGroups(container) {
		let inputGroups = container.querySelectorAll('.inputs-group');
		inputGroups.forEach((inputGroup) => inputGroup.remove());
	}

	clearInputGroups(substratesDiv);
	clearInputGroups(productsDiv);
	// let substrateNames = document.getElementById('substrateNames');
	// let productNames = document.getElementById('productNames');
	// substrateNames.innerHTML = '';
	// productNames.innerHTML = '';

	// Update reaction direction
	let direction = document.getElementById('reactionDirection');
	direction.value = data.direction;

	let subsystem = document.getElementById('subsystemField');
	subsystem.value = data.subsystem;

	let promises = [];

	if (data.substrates_names) {
		data.substrates.forEach((substrate, index) => {
			if (substrate === 'empty') {
				addField('substratesDiv', 'substrates', 'subs_sch');
			} else {
				promises.push(
					addFieldWithData(
						substratesDiv,
						'substrates',
						'subs_sch',
						substrate,
						data.subs_sch[index],
						data.subs_comps[index],
						data.subs_types[index],
						data.substrates_names[index]
					)
				);
			}
		});

		data.products.forEach((product, index) => {
			if (product === 'empty') {
				addField('productsDiv', 'products', 'prod_sch');
			} else {
				promises.push(
					addFieldWithData(
						productsDiv,
						'products',
						'prod_sch',
						product,
						data.prod_sch[index],
						data.prods_comps[index],
						data.prods_types[index],
						data.products_names[index]
					)
				);
			}
		});
	} else {
		data.substrates.forEach((substrate, index) => {
			if (substrate === 'empty') {
				addField('substratesDiv', 'substrates', 'subs_sch');
			} else {
				promises.push(
					addFieldWithData(substratesDiv, 'substrates', 'subs_sch', substrate, data.subs_sch[index], data.subs_comps[index], data.subs_types[index])
				);
			}
		});

		data.products.forEach((product, index) => {
			if (product === 'empty') {
				addField('productsDiv', 'products', 'prod_sch');
			} else {
				promises.push(
					addFieldWithData(productsDiv, 'products', 'prod_sch', product, data.prod_sch[index], data.prods_comps[index], data.prods_types[index])
				);
			}
		});
	}

	await Promise.all(promises);

	// Check URL before hiding buttons
	if (window.location.href.includes('localhost/?reaction_id=')) {
		hideDoneFieldButtons(); // Call to hide buttons
		hideDoneFieldButtonall();
	}
}

function initAutocomplete(searchInput, hiddenInput, dropdown, metabolites) {
	// Filter and display options as the user types
	searchInput.addEventListener('input', function () {
		const query = this.value.toLowerCase();
		dropdown.innerHTML = ''; // Clear existing options

		// Filter metabolites by name or abbreviation
		const filtered = metabolites.filter((m) => m.name.toLowerCase().includes(query) || (m.vmh_abbr && m.vmh_abbr.toLowerCase().includes(query)));

		filtered.forEach((m) => {
			// Create a container for this option
			const optionDiv = document.createElement('div');
			optionDiv.classList.add('autocomplete-option');

			// Create a span for the metabolite abbreviation
			const abbrSpan = document.createElement('span');
			abbrSpan.classList.add('metab-abbr');
			abbrSpan.textContent = m.vmh_abbr ? m.vmh_abbr : 'no abbr';
			if (!m.vmh_abbr) {
				abbrSpan.style.color = 'grey';
			}

			// Create a span for the metabolite name
			const nameSpan = document.createElement('span');
			nameSpan.classList.add('metab-name');
			nameSpan.textContent = m.name;

			// Append both spans to the option container (with a separator)
			optionDiv.appendChild(abbrSpan);
			optionDiv.appendChild(document.createTextNode(' - '));
			optionDiv.appendChild(nameSpan);

			// When an option is clicked, update the fields accordingly
			optionDiv.addEventListener('click', function () {
				// Show the name in the visible search input
				searchInput.value = m.name;
				// Store the metabolite id in the hidden input (this is what will be submitted)
				hiddenInput.value = m.id;
				dropdown.innerHTML = ''; // Clear the dropdown
			});

			dropdown.appendChild(optionDiv);
		});
	});

	// Show the dropdown (if applicable) when the input is focused
	searchInput.addEventListener('focus', function () {
		searchInput.dispatchEvent(new Event('input'));
	});

	// Hide the dropdown when the input loses focus (with a slight delay to allow clicks)
	searchInput.addEventListener('blur', function () {
		setTimeout(() => {
			dropdown.innerHTML = '';
		}, 200);
	});
}

function handleMetaboliteTypeChange(selectElement) {
	const group = selectElement.closest('.inputs-group');
	if (!group) return;
	const identifierCell = group.querySelector('.cell-identifier');
	if (!identifierCell) return;

	const existingAutocomplete = identifierCell.querySelector('.autocomplete-container');
	let originalInput = identifierCell.querySelector('input[type="text"]');

	if (selectElement.value === 'Saved') {
		if (existingAutocomplete) {
			return;
		}
		if (!originalInput) {
			originalInput = document.createElement('input');
			originalInput.type = 'text';
			originalInput.name = selectElement.name.replace('_type', '');
			identifierCell.appendChild(originalInput);
		}
		const autocompleteContainer = document.createElement('div');
		autocompleteContainer.classList.add('autocomplete-container');
		const searchInput = document.createElement('input');
		searchInput.type = 'text';
		searchInput.classList.add('autocomplete-input');
		searchInput.placeholder = 'Search saved metabolites...';
		const hiddenInput = document.createElement('input');
		hiddenInput.type = 'hidden';
		hiddenInput.name = originalInput.name;
		autocompleteContainer.appendChild(searchInput);
		autocompleteContainer.appendChild(hiddenInput);
		const dropdown = document.createElement('div');
		dropdown.classList.add('autocomplete-dropdown');
		autocompleteContainer.appendChild(dropdown);
		identifierCell.replaceChild(autocompleteContainer, originalInput);
		let userId = sessionStorage.getItem('userID');
		if (!window.savedMetabolitesCache) {
			fetch(`/get_saved_metabolites/?user_id=${userId}`)
				.then((response) => response.json())
				.then((data) => {
					window.savedMetabolitesCache = data.metabolites;
					initAutocomplete(searchInput, hiddenInput, dropdown, window.savedMetabolitesCache);
				});
		} else {
			initAutocomplete(searchInput, hiddenInput, dropdown, window.savedMetabolitesCache);
		}
	} else {
		if (existingAutocomplete) {
			const hiddenInput = existingAutocomplete.querySelector('input[type="hidden"]');
			const replacementInput = document.createElement('input');
			replacementInput.type = 'text';
			replacementInput.name = hiddenInput ? hiddenInput.name : selectElement.name.replace('_type', '');
			replacementInput.required = true;
			identifierCell.replaceChild(replacementInput, existingAutocomplete);
			updatePlaceholder(selectElement, replacementInput);
		}
	}
}

// Toggles between text and file input based on the selected option in the corresponding select element.
function toggleFileInput(group, selectValue) {
	const identifierCell = group ? group.querySelector('.cell-identifier') : null;
	if (!identifierCell) return;
	const textInput = identifierCell.querySelector('input[type="text"]:not(.autocomplete-input)');
	let fileInput = identifierCell.querySelector('input[type="file"]');

	function updateTextInputWithFileName() {
		if (fileInput && fileInput.files.length > 0 && textInput) {
			textInput.value = fileInput.files[0].name;
		}
	}

	if (selectValue === 'MDL Mol file') {
		if (!fileInput) {
			fileInput = document.createElement('input');
			fileInput.type = 'file';
			fileInput.name = textInput ? textInput.name : '';
			fileInput.id = textInput ? textInput.id : '';
			fileInput.style.display = 'none';
			identifierCell.appendChild(fileInput);
			fileInput.addEventListener('change', updateTextInputWithFileName);
		}
		if (textInput) {
			textInput.style.display = 'none';
		}
		fileInput.style.display = 'block';
	} else {
		if (fileInput) {
			fileInput.style.display = 'none';
		}
		if (textInput) {
			textInput.style.display = 'block';
		}
	}
}

// Handles the 'Done' button in reaction form for each metabolite
function attachEventListenersToDoneButtons() {
	const doneButtons = document.querySelectorAll('.done-field-btn');
	doneButtons.forEach((button) => {
		button.addEventListener('click', handleDoneButtonClick);
	});
}
function attachEventListenersToDoneAllButtons() {
	const doneAllButtons = document.querySelectorAll('.done-field-btn-all');
	doneAllButtons.forEach((button) => {
		button.addEventListener('click', handleDoneAllButtonClick);
	});
}

function handleDoneButtonClick(event) {
	// Handles the 'Done' button click for each metabolite in the reaction form
	const button = event.target;
	const group = button.closest('.inputs-group');
	if (!group) return;
	const parentList = group.parentElement ? group.parentElement.id : '';
	let prefix, fullname;
	if (parentList === 'substratesDiv') {
		prefix = 'subs';
		fullname = 'substrates';
	} else {
		prefix = 'prod';
		fullname = 'products';
	}
	let stoichiometryField = group.querySelector(`input[name="${prefix}_sch"]`);
	let main_input = group.querySelector(`input[name="${fullname}"]`);
	let compartmentField = group.querySelector(`select[name="${prefix}_comps"]`);
	let typeField = group.querySelector(`select[name="${fullname}_type"]`);
	if (compartmentField.value === '-') {
		alert(`Please enter a compartment for ${fullname}`);
		return;
	}
	button.disabled = true;
	const originalButtonText = button.textContent;
	button.textContent = 'Loading...';
	let data = new FormData();
	data.append('metabolite', main_input.value);
	data.append('type', typeField.value);
	data.append('compartment', compartmentField.value);
	data.append('stoichiometry', stoichiometryField.value);
	data.append('userID', sessionStorage.getItem('userID'));
	if (typeField.value === 'Saved') {
		main_input = main_input.parentElement;
	}
	let fileInputField = group.querySelector('.cell-identifier input[type="file"]');
	if (fileInputField) {
		if (fileInputField.files.length > 0) {
			data.append('file', fileInputField.files[0]);
		}
	}

	// add file in main_input if file exists
	fetch(verifyMetabolite, {
		method: 'POST',
		headers: {
			'X-Requested-With': 'XMLHttpRequest',
			'X-CSRFToken': csrfToken,
		},
		body: data,
	})
		.then((response) => response.json())
		.then((data) => {
			if (data.error) {
				showErrorModal(data.message);
				window.scrollTo(0, 0);
			} else {
				// NEW: Warn the user if a saved metabolite exists
				if (data.saved_exists) {
					alert(
						`IMPORTANT: A saved metabolite named "${data.name_in_db}" already exists. \nIf you intend to use it, please select it as a 'My Metabolite' metabolite to ensure your reaction is set up correctly.`
					);
				}
				if (data.noStructure) {
					alert(
						`No structure found for ${data.name}, only the formula (${data.formula}) will be used. \n \nThis means that some operations which require a structure will not be skipped.`
					);
				}
				updateNameFields(data, main_input, typeField, button);
				// below code will add a hidden text div that stores the status of the metabolite
				group.dataset.atomCounts = JSON.stringify(data.atom_counts);
				group.dataset.charge = data.charge;

				let hidden_status = document.createElement('input');
				hidden_status.style.display = 'none';
				hidden_status.value = data.found;
				hidden_status.className = 'valid-status';
				group.appendChild(hidden_status);

				updateAtomChargeCounters();
			}
		})
		.finally(() => {
			button.disabled = false;
			button.textContent = originalButtonText;
		});
}

function handleDoneAllButtonClick(event) {
	const button = event.target;
	const section = button.closest('.reactant-section');
	if (!section) return;
	const container = section.querySelector('.reactant-rows');
	if (!container) return;
	const isSubstrates = container.id === 'substratesDiv';
	const prefix = isSubstrates ? 'subs' : 'prod';
	const fullname = isSubstrates ? 'substrates' : 'products';

	let compartmentFields = container.querySelectorAll(`select[name="${prefix}_comps"]`);

	// Check all compartment fields
	for (let i = 0; i < compartmentFields.length; i++) {
		if (compartmentFields[i].value === '-') {
			alert(`Please enter a compartment for ${fullname}`);
			return;
		}
	}

	const inputsGroupsInConrainer = Array.from(container.querySelectorAll('.inputs-group'));
	inputsGroupsInConrainer.forEach((group) => {
		const doneButton = group.querySelector('.done-field-btn');
		const validStatus = group.querySelector('.valid-status');
		if (doneButton && validStatus == null) {
			doneButton.click();
		}
	});
}

function updateNameFields(data = {}, main_input = null, typeField = null, button, showButton = false) {
	const group = button ? button.closest('.inputs-group') : null;
	if (!group) return;
	const nameField = group.querySelector('[class*="name"]');
	if (nameField) {
		nameField.style.display = 'block';
	}

	if (main_input && typeField) {
		if (main_input instanceof HTMLElement) {
			main_input.disabled = true;
		}
		typeField.disabled = true;
		if (typeField.value === 'Saved') {
			const textChild = main_input.querySelector && main_input.querySelector('input[type="text"]');
			if (textChild) {
				textChild.disabled = true;
			}
		}
	}
	button.style.display = showButton ? 'inline-flex' : 'none';

	let statusDot = group.querySelector('.status-dot');
	if (statusDot) {
		statusDot.style.display = 'block';
	}

	if (data && typeof updateStatusDot === 'function') {
		updateStatusDot(statusDot, data.found, data.miriam);
		if (data.name && nameField) {
			nameField.value = data.name;
		}
	}

	if (main_input && typeField) {
		let editDrawingButton = typeField.closest('.cell-type')?.querySelector('.edit-drawing-btn');
		let startDrawingButton = typeField.closest('.cell-type')?.querySelector('.start-drawing-btn');
		if (data.found) {
			if (typeField.value === 'Draw' && main_input.style) {
				main_input.style.visibility = 'visible';
			}
			if (main_input instanceof HTMLInputElement) {
				main_input.value = data.abbr;
			}
			typeField.value = 'VMH';
			toggleFileInput(group, 'VMH');
			if (nameField) {
				nameField.disabled = true;
			}
			if (editDrawingButton) {
				editDrawingButton.remove();
			}
			if (startDrawingButton) {
				startDrawingButton.remove();
			}
		} else {
			if (nameField) {
				if (typeField.value === 'Saved') {
					nameField.disabled = true;
				} else {
					nameField.disabled = false;
				}
			}
			if (editDrawingButton) {
				editDrawingButton.disabled = true;
			}
			if (startDrawingButton) {
				startDrawingButton.disabled = true;
			}
		}
	}
}
// Format an atom counts object into a string (e.g., "C: 6, H: 12, O: 6")
function formatAtomCounts(counts) {
	let parts = [];
	for (let elem in counts) {
		if (counts.hasOwnProperty(elem)) {
			parts.push(`${elem}: ${counts[elem]}`);
		}
	}
	return parts.join(', ');
}

// Aggregate and update the counters for substrates and products
function updateAtomChargeCounters() {
	// Aggregate from the form (same logic you had)
	const substratesDiv = document.getElementById('substratesDiv');
	const productsDiv = document.getElementById('productsDiv');

	const totalAtomsSubs = {};
	const totalAtomsProds = {};
	let totalChargeSubs = 0;
	let totalChargeProds = 0;

	(substratesDiv ? substratesDiv.querySelectorAll('.inputs-group') : []).forEach((group) => {
		if (group.dataset.atomCounts) {
			const counts = JSON.parse(group.dataset.atomCounts);
			for (const elem in counts) totalAtomsSubs[elem] = (totalAtomsSubs[elem] || 0) + counts[elem];
			totalChargeSubs += parseFloat(group.dataset.charge) || 0;
		}
	});

	(productsDiv ? productsDiv.querySelectorAll('.inputs-group') : []).forEach((group) => {
		if (group.dataset.atomCounts) {
			const counts = JSON.parse(group.dataset.atomCounts);
			for (const elem in counts) totalAtomsProds[elem] = (totalAtomsProds[elem] || 0) + counts[elem];
			totalChargeProds += parseFloat(group.dataset.charge) || 0;
		}
	});

	// Balanced?
	const atoms = new Set([...Object.keys(totalAtomsSubs), ...Object.keys(totalAtomsProds)]);
	let massBalanced = true;
	for (const a of atoms) {
		if ((totalAtomsSubs[a] || 0) !== (totalAtomsProds[a] || 0)) {
			massBalanced = false;
			break;
		}
	}
	const chargeBalanced = totalChargeSubs === totalChargeProds;

	// Update header lines + charges
	renderChemInfoHeader({
		massBalanced,
		chargeBalanced,
		subsAtoms: totalAtomsSubs,
		prodsAtoms: totalAtomsProds,
		subsCharge: totalChargeSubs,
		prodsCharge: totalChargeProds,
	});

	// Update the ORIGINAL table in-place
	renderAtomComparisonTableInPlace(totalAtomsSubs, totalAtomsProds, {});
}
