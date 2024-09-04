function formatTooltipContent(element, type, compartment) {
    // Concatenate only the 'element' part if its length exceeds 20 characters
    let formattedElement = element.length > 20 ? element.substring(0, 17) + '...' : element;
    let content = `${formattedElement} (${type}), Compartment: ${compartment}`;
    return content;
}

function setupTooltips() {
    document.querySelectorAll('.detail-item,.info-symbol').forEach(item => {
        item.addEventListener('mouseenter', function() {
            const tooltipContent = this.getAttribute('data-tooltip-content');
            const tooltip = document.createElement('div');
            tooltip.className = 'tooltip';
            tooltip.innerHTML = tooltipContent;
            this.appendChild(tooltip);
        });
        item.addEventListener('mouseleave', function() {
            this.removeChild(this.querySelector('.tooltip'));
        });
    });
}

function validateInputs() {
    // Check for any empty non-VMH input fields
    let allInputsValid = true;
    document.querySelectorAll('.sub-name-input, .prod-name-input').forEach(input => {
        if (input.value.trim() === '') {
            allInputsValid = false;
        }
    });

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
