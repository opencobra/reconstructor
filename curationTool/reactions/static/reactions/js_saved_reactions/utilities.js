function formatTooltipContent(element, type, compartment) {
	// Concatenate only the 'element' part if its length exceeds 20 characters
	let formattedElement = element.length > 20 ? element.substring(0, 17) + '...' : element;
	let content = `${formattedElement} (${type}), Compartment: ${compartment}`;
	return content;
}

function setupTooltips() {
	document.querySelectorAll('.detail-item,.info-symbol').forEach((item) => {
		item.addEventListener('mouseenter', function () {
			const tooltipContent = this.getAttribute('data-tooltip-content');
			const tooltip = document.createElement('div');
			tooltip.className = 'tooltip';
			tooltip.innerHTML = tooltipContent;
			this.appendChild(tooltip);
		});
		item.addEventListener('mouseleave', function () {
			this.removeChild(this.querySelector('.tooltip'));
		});
	});
}

function validateInputs() {
	// Check for any empty non-VMH input fields
	let allInputsValid = true;
	
	// Clear previous validation states
	document.querySelectorAll('.ws-input-error').forEach(el => {
		el.classList.remove('ws-input-error');
	});
	document.querySelectorAll('.ws-row-error').forEach(el => {
		el.classList.remove('ws-row-error');
	});
	
	// Check substrate and product name inputs
	document.querySelectorAll('.sub-name-input, .prod-name-input').forEach((input) => {
		if (input.value.trim() === '') {
			allInputsValid = false;
			input.classList.add('ws-input-error');
		}
	});

	// Also check inputs in workspace detail panel
	const detailPanel = document.getElementById('wsAvailableReactionDetails');
	if (detailPanel) {
		// Check name inputs that aren't readonly (new metabolites need names)
		detailPanel.querySelectorAll('input[name="subsNameInput"]:not([readonly]), input[name="prodsNameInput"]:not([readonly])').forEach((input) => {
			if (input.value.trim() === '') {
				allInputsValid = false;
				input.classList.add('ws-input-error');
				// Highlight the row
				const row = input.closest('tr');
				if (row) row.classList.add('ws-row-error');
			}
		});
		
		// Check abbreviation inputs that aren't readonly
		detailPanel.querySelectorAll('input[name="subsAbbrInput"]:not([readonly]), input[name="prodsAbbrInput"]:not([readonly])').forEach((input) => {
			if (input.value.trim() === '') {
				allInputsValid = false;
				input.classList.add('ws-input-error');
				// Highlight the row
				const row = input.closest('tr');
				if (row) row.classList.add('ws-row-error');
			}
		});
		
		// Check reaction abbreviation
		const reactionNameInput = detailPanel.querySelector('.reaction-name-input');
		if (reactionNameInput && reactionNameInput.value.trim() === '') {
			allInputsValid = false;
			reactionNameInput.classList.add('ws-input-error');
		}

		const reactionAbbrInput = detailPanel.querySelector('.reaction-abbreviation-input');
		if (reactionAbbrInput && reactionAbbrInput.value.trim() === '') {
			allInputsValid = false;
			reactionAbbrInput.classList.add('ws-input-error');
		}

		detailPanel.querySelectorAll('.ws-availability-status[data-state="conflict"]').forEach((status) => {
			allInputsValid = false;
			const row = status.closest('tr');
			const formSection = status.closest('.ws-form-section');
			if (row) row.classList.add('ws-row-error');
			if (formSection) {
				formSection.querySelectorAll('.reaction-name-input, .reaction-abbreviation-input').forEach((input) => {
					input.classList.add('ws-input-error');
				});
			}
		});
	}

	return allInputsValid;
}

function displayValidationMessage(display, message = '') {
	let messageContainer = document.getElementById('validationMessage');
	const modalContent = document.getElementById('reactionModal').querySelector('.modal-content');
	var modal = document.getElementById('reactionModal');
	if (!messageContainer) {
		messageContainer = document.createElement('div');
		messageContainer.id = 'validationMessage';
		messageContainer.style.color = 'red';
		messageContainer.style.textAlign = 'center'; // Center the message for better visibility
		messageContainer.style.padding = '10px 0'; // Add some padding for spacing
	}

	messageContainer.textContent = message;

	if (display) {
		// Insert the message at the top of the modal content
		modal.scrollTo(0, 0);
		if (modalContent.firstChild) {
			modalContent.insertBefore(messageContainer, modalContent.firstChild);
		} else {
			modalContent.appendChild(messageContainer);
		}
	} else {
		if (messageContainer.parentNode) {
			messageContainer.parentNode.removeChild(messageContainer);
		}
	}
}
function setButtonState(isDisabled) {
	const confirmButton = document.getElementById('confirmAddToVMH');
	if (isDisabled) {
		confirmButton.classList.add('button-disabled');
		confirmButton.disabled = true;
	} else {
		confirmButton.classList.remove('button-disabled');
		confirmButton.disabled = false;
	}
}

function rgbToHex(rgb) {
	const result = rgb.match(/\d+/g);
	const r = parseInt(result[0]).toString(16).padStart(2, '0');
	const g = parseInt(result[1]).toString(16).padStart(2, '0');
	const b = parseInt(result[2]).toString(16).padStart(2, '0');
	return `#${r}${g}${b}`;
}
