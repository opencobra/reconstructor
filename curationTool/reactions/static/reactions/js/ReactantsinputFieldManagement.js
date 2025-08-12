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
				toggleFileInput(select.parentElement, selectedValue);
				handleSelectChange({ target: select });
				updatePlaceholder(select, select.parentElement.querySelector('input[type="text"]'));
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
				toggleFileInput(select.parentElement, selectedValue);
				handleSelectChange({ target: select });
				updatePlaceholder(select, select.parentElement.querySelector('input[type="text"]'));
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
	const imgElement = event.currentTarget.querySelector('img');
	imgElement.classList.toggle('rotated');
}
document.getElementById('toggleApplyAllsubs').addEventListener('click', function () {
	var applyOptions = document.getElementById('applyAllOptionsSubs');
	applyOptions.style.display = applyOptions.style.display === 'none' ? 'flex' : 'none';
});
document.getElementById('toggleApplyAllprods').addEventListener('click', function () {
	var applyOptions = document.getElementById('applyAllOptionsProds');
	applyOptions.style.display = applyOptions.style.display === 'none' ? 'flex' : 'none';
});

document.getElementById('toggleApplyAllsubs').addEventListener('click', toggleRotation);

document.getElementById('toggleApplyAllprods').addEventListener('click', toggleRotation);

// Initializes event listeners on existing select elements for file input toggling upon DOM content load.
document.addEventListener('DOMContentLoaded', function () {
	// Select both 'products_type' and 'substrates_type' select elements
	document.querySelectorAll('select[name="products_type"], select[name="substrates_type"]').forEach((selectElement) => {
		// Initial call to handle the current state
		toggleFileInput(selectElement.parentElement, selectElement.value);
		handleMetaboliteTypeChange(selectElement);

		// Attach change event listener to these select elements
		selectElement.addEventListener('change', function () {
			toggleFileInput(selectElement.parentElement, selectElement.value);
			handleMetaboliteTypeChange(selectElement);
		});
	});
});

// Dynamically adds a new field container with specified inputs and select options to the given container ID.
function addField(containerId, inputName, numberName) {
	var inputsGroup = document.createElement('div');
	inputsGroup.className = 'inputs-group';

	var removeBtn = document.createElement('button');
	removeBtn.type = 'button'; // Prevent submission on click
	removeBtn.className = 'ui inverted red button remove-field-btn';
	removeBtn.textContent = 'X';
	removeBtn.onclick = function () {
		removeField(this);
	};
	inputsGroup.appendChild(removeBtn);

	var newNumberInput = document.createElement('input');
	newNumberInput.setAttribute('type', 'number');
	newNumberInput.setAttribute('name', numberName);
	newNumberInput.setAttribute('id', numberName);
	newNumberInput.setAttribute('min', '1');
	newNumberInput.value = '1'; // Default stoichiometry value
	inputsGroup.appendChild(newNumberInput);

	var newInput = document.createElement('input');
	newInput.setAttribute('type', 'text');
	newInput.setAttribute('name', inputName);
	newInput.setAttribute('id', inputName);
	inputsGroup.appendChild(newInput);

	if (containerId == 'substratesDiv') {
		var prod_or_subs = 'subs_comps';
	} else {
		var prod_or_subs = 'prod_comps';
	}

	var compartmentSelect = createCompartmentSelect(prod_or_subs);
	inputsGroup.appendChild(compartmentSelect);

	var selectInputType = document.createElement('select');
	selectInputType.name = inputName + '_type'; // e.g., substrates_type or products_type
	selectInputType.id = inputName + '_type';

	['VMH', 'ChEBI ID', 'ChEBI Name', 'SwissLipids', 'MDL Mol file', 'Draw', 'PubChem ID', 'Saved'].forEach(function (type) {
		var option = document.createElement('option');
		if (type === 'Saved') {
			option.text = 'My Metabolites';
		} else {
			option.text = type;
		}
		option.value = type;
		selectInputType.appendChild(option);
	});
	inputsGroup.appendChild(selectInputType);

	// Event listener to toggle file input visibility and handle metabolite type changes
	selectInputType.addEventListener('change', function () {
		toggleFileInput(inputsGroup, this.value);
		handleMetaboliteTypeChange(this);
	});
	toggleFileInput(inputsGroup, selectInputType.value); // Initial toggle based on default select value

	var statusDot = document.createElement('span');
	statusDot.className = 'status-dot';
	statusDot.style.display = 'none';
	inputsGroup.appendChild(statusDot);

	updatePlaceholder(selectInputType, newInput);

	var doneButton = document.createElement('button');
	doneButton.type = 'button';
	doneButton.className = 'done-field-btn';
	doneButton.textContent = 'verify';
	inputsGroup.appendChild(doneButton);

	// Add hidden text input field
	var hiddenInput = document.createElement('input');
	hiddenInput.type = 'text';
	hiddenInput.className = inputName + '-name';
	hiddenInput.name = inputName + '_name';
	hiddenInput.style.display = 'none';
	inputsGroup.appendChild(hiddenInput);

	document.getElementById(containerId).appendChild(inputsGroup);

	attachEventListenersToSelects();
	attachEventListenersToDoneButtons();
}

// Adds a field to the specified container with pre-filled data for substrates or products.
function addFieldWithData(container, name, schName, value, schValue, compValue, type, metab_name = null, metab_id = null) {
	return new Promise(async (resolve, reject) => {
		var inputsGroup = document.createElement('div');
		inputsGroup.className = 'inputs-group';

		// Create and attach remove button
		var removeBtn = document.createElement('button');
		removeBtn.type = 'button';
		removeBtn.className = 'ui inverted red button remove-field-btn';
		removeBtn.textContent = 'X';
		removeBtn.onclick = function () {
			removeField(this);
		};
		inputsGroup.appendChild(removeBtn);

		// Create stoichiometry input
		var newNumberInput = document.createElement('input');
		newNumberInput.type = 'number';
		newNumberInput.name = schName;
		newNumberInput.id = schName;
		newNumberInput.min = '1';
		newNumberInput.value = schValue;
		inputsGroup.appendChild(newNumberInput);

		// Create main input with pre-filled value
		var newInput = document.createElement('input');
		newInput.type = 'text';
		newInput.name = name;
		newInput.id = name;
		newInput.value = value;
		inputsGroup.appendChild(newInput);

		// Create compartment select (e.g., subs_comps or prod_comps)
		var prod_or_subs = container.id == 'substratesDiv' ? 'subs_comps' : 'prod_comps';
		var compartmentSelect = createCompartmentSelect(prod_or_subs);
		compartmentSelect.value = compValue;
		inputsGroup.appendChild(compartmentSelect);

		// Create type select
		var selectInputType = document.createElement('select');
		selectInputType.name = `${name}_type`;
		selectInputType.id = `${name}_type`;
		['VMH', 'ChEBI ID', 'ChEBI Name', 'SwissLipids', 'MDL Mol file', 'Draw', 'PubChem ID', 'Saved'].forEach(function (optionType) {
			var option = document.createElement('option');
			option.value = optionType;
			option.text = optionType === 'Saved' ? 'My Metabolites' : optionType;
			if (type === optionType) {
				option.selected = true;
			}
			selectInputType.appendChild(option);
		});
		inputsGroup.appendChild(selectInputType);

		// Attach the change event listener for toggling file input and handling type changes.
		selectInputType.addEventListener('change', function () {
			toggleFileInput(inputsGroup, this.value);
			handleMetaboliteTypeChange(this);
		});
		toggleFileInput(inputsGroup, selectInputType.value); // Initial toggle

		// Create status dot and hidden text field for the metabolite name (if needed)
		var statusDot = document.createElement('span');
		statusDot.className = 'status-dot';
		statusDot.style.display = 'none';
		inputsGroup.appendChild(statusDot);

		var hiddenInput = document.createElement('input');
		hiddenInput.type = 'text';
		var inputName = container.id == 'substratesDiv' ? 'substrates' : 'products';
		hiddenInput.className = inputName + '-name';
		hiddenInput.name = inputName + '_name';
		hiddenInput.style.display = 'none';
		if (metab_name) {
			hiddenInput.value = metab_name;
		}
		inputsGroup.appendChild(hiddenInput);

		// Append the inputsGroup to the container
		container.appendChild(inputsGroup);

		// If the field type is prefilled as 'MDL Mol file', handle async loading
		if (type === 'MDL Mol file') {
			var fileInput = inputsGroup.querySelector('input[type="file"]');
			try {
				await loadFileToInputField(value, fileInput);
				resolve();
			} catch (error) {
				reject(error);
			}
		} else {
			resolve();
		}

		updatePlaceholder(selectInputType, newInput);

		var doneButton = document.createElement('button');
		doneButton.type = 'button';
		doneButton.className = 'done-field-btn';
		doneButton.textContent = 'verify';
		inputsGroup.appendChild(doneButton);

		// If the prefilled type is 'Saved', simulate the transformation into an autocomplete widget.
		if (type === 'Saved') {
			var group = selectInputType.closest('.inputs-group');
			var autocompleteContainer = group.querySelector('.autocomplete-container');
			if (autocompleteContainer) {
				var searchInput = autocompleteContainer.querySelector('.autocomplete-input');
				var hiddenInputField = autocompleteContainer.querySelector('input[type="hidden"]');
				// Set the visible search input (display value) and the hidden input (ID).
				if (searchInput) {
					searchInput.value = metab_name;
				}
				if (hiddenInputField) {
					// Ensure you pass the correct saved metabolite ID (e.g. metab_id).
					hiddenInputField.value = value;
				}
			}
		}

		// Attach any additional event listeners.
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
	var group = selectElement.closest('.inputs-group');
	// Find the original text input that normally collects the metabolite value.
	// (Assumes the first text input in the group is the one to replace.)
	var originalInput = group.querySelector('input[type="text"]');

	if (selectElement.value === 'Saved') {
		// Create a container to hold the autocomplete elements.
		var autocompleteContainer = document.createElement('div');
		autocompleteContainer.classList.add('autocomplete-container');

		// Create the visible search input (without a name so it isn't submitted).
		var searchInput = document.createElement('input');
		searchInput.type = 'text';
		searchInput.classList.add('autocomplete-input');
		searchInput.placeholder = 'Search saved metabolites...';

		// Create a hidden input that carries the same name as the original.
		var hiddenInput = document.createElement('input');
		hiddenInput.type = 'hidden';
		hiddenInput.name = originalInput.name; // e.g., "substrates" or "products"

		// Replace the original input with the autocomplete container.
		originalInput.parentNode.replaceChild(autocompleteContainer, originalInput);
		autocompleteContainer.appendChild(searchInput);
		autocompleteContainer.appendChild(hiddenInput);

		// Create a container for the autocomplete dropdown.
		var dropdown = document.createElement('div');
		dropdown.classList.add('autocomplete-dropdown');
		autocompleteContainer.appendChild(dropdown);
		let userId = sessionStorage.getItem('userID');
		// Fetch the saved metabolites, caching them globally.
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
		// If the select value is not "Saved" and an autocomplete container exists, revert.
		var container = group.querySelector('.autocomplete-container');
		if (container) {
			var newInput = document.createElement('input');
			newInput.type = 'text';
			// Set the name from the hidden input in the autocomplete container.
			newInput.name = container.querySelector('input[type="hidden"]').name;
			newInput.required = true;
			container.parentNode.replaceChild(newInput, container);
			// Optionally update the placeholder based on the select's current value.
			updatePlaceholder(selectElement, newInput);
		}
	}
}

// Toggles between text and file input based on the selected option in the corresponding select element.
function toggleFileInput(container, selectValue) {
	let textInput = container.querySelector('input[type="text"]');
	let fileInput = container.querySelector('input[type="file"]');
	let selectElement = container.querySelector('select');

	function updateTextInputWithFileName() {
		if (fileInput.files.length > 0) {
			textInput.value = fileInput.files[0].name; // Update text input with file name
		}
	}

	if (selectValue === 'MDL Mol file') {
		if (!fileInput) {
			fileInput = document.createElement('input');
			fileInput.type = 'file';
			fileInput.name = textInput.name;
			fileInput.id = textInput.name;
			fileInput.style.display = 'none'; // Hide initially
			container.insertBefore(fileInput, selectElement); // Insert before the select element
			fileInput.addEventListener('change', updateTextInputWithFileName);
			textInput.style.display = 'none'; // Keep the text input hidden
			fileInput.style.display = 'block';
			updateTextInputWithFileName(); // Update text input in case a file is already selected
		}
		textInput.style.display = 'none';
		fileInput.style.display = 'block';
	} else {
		if (fileInput) fileInput.style.display = 'none';
		textInput.style.display = 'block';
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
	let button = event.target;
	if (button.parentNode.parentNode.id === 'substratesDiv') {
		prefix = 'subs';
		fullname = 'substrates';
	} else {
		prefix = 'prod';
		fullname = 'products';
	}
	let stoichiometryField = button.parentNode.querySelector(`input[name="${prefix}_sch"]`);
	let main_input = button.parentNode.querySelector(`input[name="${fullname}"]`);
	let compartmentField = button.parentNode.querySelector(`select[name="${prefix}_comps"]`);
	let typeField = button.parentNode.querySelector(`select[name="${fullname}_type"]`);
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
	let fileInputField = button.parentNode.querySelector('input[type="file"]');
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
				main_input.parentNode.dataset.atomCounts = JSON.stringify(data.atom_counts);
				main_input.parentNode.dataset.charge = data.charge;

				let hidden_status = document.createElement('input');
				hidden_status.style.display = 'none';
				hidden_status.value = data.found;
				hidden_status.className = 'valid-status';
				main_input.parentNode.appendChild(hidden_status);

				updateAtomChargeCounters();
			}
		})
		.finally(() => {
			button.disabled = false;
			button.textContent = originalButtonText;
		});
}

function handleDoneAllButtonClick(event) {
	let button = event.target;
	let prefix, fullname;

	if (button.parentNode.parentNode.id === 'substratesDiv') {
		prefix = 'subs';
		fullname = 'substrates';
	} else {
		prefix = 'prod';
		fullname = 'products';
	}

	let container = button.parentNode.parentNode;
	let compartmentFields = container.querySelectorAll(`select[name="${prefix}_comps"]`);

	// Check all compartment fields
	for (let i = 0; i < compartmentFields.length; i++) {
		if (compartmentFields[i].value === '-') {
			alert(`Please enter a compartment for ${fullname}`);
			return;
		}
	}

	let inputsGroupsInConrainer = container.querySelectorAll('.inputs-group');
	inputsGroupsInConrainer = Array.from(inputsGroupsInConrainer);

	for (let i = 0; i < inputsGroupsInConrainer.length; i++) {
		let doneButton = inputsGroupsInConrainer[i].querySelector('.done-field-btn');
		let validStatus = inputsGroupsInConrainer[i].querySelector('.valid-status');
		if (validStatus == null) {
			doneButton.click();
		}
	}
}

function updateNameFields(data = {}, main_input = null, typeField = null, button, showButton = false) {
	var nameField = main_input.parentNode.querySelector('[class*="name"]');
	nameField.style.display = 'block';

	if (main_input && typeField) {
		main_input.disabled = true;
		typeField.disabled = true;
		if (typeField.value === 'Saved') {
			const textChild = main_input.querySelector('input[type="text"]');
			if (textChild) {
				textChild.disabled = true;
			}
		}
	}
	button.style.display = 'none';

	let statusDot = button.parentNode.querySelector('.status-dot');
	if (statusDot) {
		statusDot.style.display = 'block';
	}

	if (data && typeof updateStatusDot === 'function') {
		updateStatusDot(statusDot, data.found, data.miriam);
		if (data.name && nameField) {
			nameField.value = data.name;
		}
	}
	if (data.name) {
		updateStatusDot(statusDot, data.found, data.miriam);
		if (data.name && nameField) {
			nameField.value = data.name;
		}
	}
	if (main_input && typeField) {
		let editDrawingButton = typeField.parentNode.querySelector('.edit-drawing-btn');
		let startDrawingButton = typeField.parentNode.querySelector('.start-drawing-btn');
		if (data.found) {
			if (typeField.value === 'Draw') {
				main_input.style.visibility = 'visible';
			}
			main_input.value = data.abbr;
			typeField.value = 'VMH';
			toggleFileInput(button.parentNode, 'VMH');
			nameField.disabled = true;
			if (editDrawingButton) {
				editDrawingButton.remove();
			}
			if (startDrawingButton) {
				startDrawingButton.remove();
			}
		} else {
			if (typeField.value === 'Saved') {
				nameField.disabled = true;
			} else {
				nameField.disabled = false;
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
