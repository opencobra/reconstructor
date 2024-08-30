

function createGeneInfoInput() {

    // Create a main container for the gene info input
    var container = document.createElement('div');
    container.id = 'geneInfoContainer';
    container.className = 'gene-info-container';

    // Create the gene input element
    var geneInput = document.createElement('div');
    geneInput.contentEditable = 'true'; // Make it editable
    geneInput.className = 'gene-info-input';
    geneInput.id = 'geneInfoInput';
    geneInput.style.border = '1px solid #ccc';
    geneInput.style.padding = '5px';
    geneInput.style.minHeight = '20px';
    geneInput.style.display = 'block';
    geneInput.style.width = '100%';
    geneInput.style.overflow = 'hidden';
    geneInput.style.boxSizing = 'border-box';
    geneInput.style.marginTop = '16px';
    geneInput.style.height = '39px';

    // Append the gene input to the main container
    container.appendChild(geneInput);

    // Optionally, create a button container if you need the AND, OR, etc., buttons
    var geneButtonContainer = document.createElement('div');
    geneButtonContainer.className = 'gene-button-container';
    ['Gene', 'AND', 'OR', '(', ')'].forEach(function (label) {
        var button = createButton(label);
        geneButtonContainer.appendChild(button);
    });

    // Append the button container to the main container if buttons are needed
    container.appendChild(geneButtonContainer);

    // Finally, append the container to the geneInfoWrapper
    document.getElementById('geneInfoWrapper').appendChild(container);
}

function createButton(label) {
    var button = document.createElement('button');
    button.textContent = label;
    button.className = 'gene-info-button';
    button.style.margin = '5px'; // Style the button (optional)

    // Add event listener to handle the button's action
    button.addEventListener('click', function () {
        handleButtonClick(label);
    });

    return button;
}
function handleButtonClick(label) {
    if (label === 'AND' || label === 'OR') {
        var cursorPosition = getCursorPosition();
        insertTextAtCursor(" " + label + " ", cursorPosition);
        placeCursorAtEnd(document.getElementById('geneInfoInput'));

    } else if (label === '(' || label === ')') {
        var cursorPosition = getCursorPosition();
        insertTextAtCursor(label, cursorPosition);
        placeCursorAtEnd(document.getElementById('geneInfoInput'));

    } else if (label === 'Gene') {
        var cursorPosition = getCursorPosition();
        createAndShowModal((inputValue, typeValue) => {
            return submitGeneInfo(inputValue, typeValue)
                .then(symbol => {
                    // Store the returned symbol in a variable
                    const returnedSymbol = symbol;
                    insertTextAtCursor(returnedSymbol, cursorPosition);
                    placeCursorAtEnd(document.getElementById('geneInfoInput'));

                    // You can now use the returnedSymbol variable as needed
                });
        });
    }
}

function getCursorPosition() {
    var selection = window.getSelection();

    if (selection.rangeCount > 0) {
        var range = selection.getRangeAt(0);
        var cursorPosition = {
            startContainer: range.startContainer,
            startOffset: range.startOffset,
            endContainer: range.endContainer,
            endOffset: range.endOffset
        };
        return cursorPosition;
    }
}

function disableKeyboardInput() {
    // Get the input element by its ID
    var inputElement = document.getElementById('geneInfoInput');

    // Function to determine if the key should be allowed
    function shouldAllowKey(event) {
        // Key codes: 8 (Backspace), 37 (Left Arrow), 39 (Right Arrow)
        var allowedKeys = [8, 37, 39];
        return allowedKeys.includes(event.keyCode);
    }

    // Add event listeners to prevent default behavior for key events, except for allowed keys
    inputElement.addEventListener('keypress', function (event) {
        if (!shouldAllowKey(event)) {
            event.preventDefault();
        }
    });

    inputElement.addEventListener('keydown', function (event) {
        if (!shouldAllowKey(event)) {
            event.preventDefault();
        }
    });

    inputElement.addEventListener('keyup', function (event) {
        if (!shouldAllowKey(event)) {
            event.preventDefault();
        }
    });
}

function insertTextAtCursor(textToInsert, cursorPosition) {
    // Select the input element
    var input = document.getElementById('geneInfoInput');
    var textSpan = document.createElement('span');
    textSpan.className = 'gene-block';
    textSpan.contentEditable = false; // Ensure the inserted text is not editable
    textSpan.innerHTML = textToInsert;

    if (cursorPosition) {
        var range = document.createRange();

        // Set the range according to the cursor position
        range.setStart(cursorPosition.startContainer, cursorPosition.startOffset);
        range.setEnd(cursorPosition.endContainer, cursorPosition.endOffset);

        if (range.endContainer === input || range.endContainer.parentNode === input) {
            range.deleteContents(); // Remove any selected content
            range.insertNode(textSpan);
            range.setStartAfter(textSpan); // Move the range start after the inserted text
            range.setEndAfter(textSpan); // Move the range end after the inserted text

            // Move the cursor to the end of the inserted content
            var selection = window.getSelection();
            selection.removeAllRanges();
            selection.addRange(range);
        } else {
            alert('Please click inside the input field to insert the gene symbol');
        }
    } else {
        console.log('No cursor position found.');
    }
}

function placeCursorAtEnd(element) {
    element.focus(); // Ensure the element is focused
    
    // Create a new range and set it to the end of the content
    var range = document.createRange();
    range.selectNodeContents(element);
    range.collapse(false); // Collapse the range to the end point

    // Apply the range to the selection
    var selection = window.getSelection();
    selection.removeAllRanges();
    selection.addRange(range);
}


function createAndShowModal(onConfirm) {
    // Check if modal already exists and remove it
    let existingModal = document.getElementById('geneModal');
    if (existingModal) {
        existingModal.remove();
    }

    // Create modal elements
    const modal = document.createElement('div');
    modal.id = 'geneModal';
    modal.className = 'modal';
    modal.style.position = 'fixed';
    modal.style.top = '50%';
    modal.style.left = '50%';
    modal.style.transform = 'translate(-50%, -50%)';
    modal.style.zIndex = '1000';

    const modalContent = document.createElement('div');
    modalContent.className = 'modal-content';
    modalContent.style.backgroundColor = '#fff';
    modalContent.style.padding = '20px';
    modalContent.style.borderRadius = '5px';
    modalContent.style.boxShadow = '0 2px 10px rgba(0, 0, 0, 0.1)';
    modalContent.style.position = 'relative'; // To position the close button absolutely within modal content

    const closeButton = document.createElement('button');
    closeButton.innerText = 'X';
    closeButton.style.background = 'none';
    closeButton.style.border = 'none';
    closeButton.style.fontSize = '16px';
    closeButton.style.fontWeight = 'bold'; // Make the 'X' bold
    closeButton.style.color = '#000'; // Ensure the color is black
    closeButton.style.cursor = 'pointer';
    closeButton.style.position = 'absolute';
    closeButton.style.top = '10px';
    closeButton.style.right = '10px';

    // Add event listener to close the modal when 'X' button is clicked
    closeButton.onclick = function () {
        modal.remove();
    };

    const modalHeader = document.createElement('h4');
    modalHeader.innerText = 'Enter Gene Information';
    modalHeader.style.marginBottom = '20px';

    const dropdown = document.createElement('select');
    dropdown.id = 'geneTypeDropdown';
    dropdown.style.marginBottom = '10px';
    dropdown.style.width = '100%';
    dropdown.style.padding = '10px';
    dropdown.style.border = '1px solid #ccc';
    dropdown.style.borderRadius = '4px';

    const option1 = document.createElement('option');
    option1.value = 'HGNC Symbol';
    option1.text = 'HGNC Symbol';
    const option2 = document.createElement('option');
    option2.value = 'Entrez ID';
    option2.text = 'Entrez ID';
    dropdown.add(option1);
    dropdown.add(option2);

    const geneInputModal = document.createElement('input');
    geneInputModal.className = 'gene-modal-input';
    geneInputModal.type = 'text';
    geneInputModal.id = 'geneInputModal' + Math.random();

    // Apply the same styles as the dropdown element
    geneInputModal.style.marginBottom = '10px';
    geneInputModal.style.width = '98%';
    geneInputModal.style.padding = '10px';
    geneInputModal.style.border = '1px solid #ccc';
    geneInputModal.style.borderRadius = '4px';

    dropdown.onchange = function () {
        geneInputModal.placeholder = 'Enter ' + dropdown.value;
    };

    const modalFooter = document.createElement('div');
    modalFooter.className = 'modal-footer';
    modalFooter.style.display = 'flex'; // Use flexbox
    modalFooter.style.justifyContent = 'center'; // Center the content
    modalFooter.style.marginTop = '20px';

    const confirmButton = document.createElement('button');
    confirmButton.id = 'modalConfirmButton';
    confirmButton.innerText = 'Confirm';
    confirmButton.style.padding = '10px 20px';
    confirmButton.style.backgroundColor = '#28a745';
    confirmButton.style.color = '#fff';
    confirmButton.style.border = 'none';
    confirmButton.style.borderRadius = '4px';
    confirmButton.style.cursor = 'pointer';

    confirmButton.onclick = function () {
        const selectedType = dropdown.value;
        const inputValue = geneInputModal.value.trim();
        const typeValue = selectedType === 'Entrez ID' ? 'Entrez ID' : 'HGNC Symbol';

        if (!inputValue) {
            console.error("Input value is undefined or empty");
            alert("Please enter a gene value");
            return;
        }

        // Call the onConfirm callback with the input values and handle the returned promise
        onConfirm(inputValue, typeValue, 'geneInfoInput', confirmButton.id);
        modal.style.display = 'none';
    };

    // Append elements to modal content
    modalContent.appendChild(closeButton); // Append close button
    modalContent.appendChild(modalHeader);
    modalContent.appendChild(dropdown);
    modalContent.appendChild(geneInputModal);
    modalFooter.appendChild(confirmButton);
    modalContent.appendChild(modalFooter);
    modal.appendChild(modalContent);
    document.body.appendChild(modal);

    // Display the modal
    modal.style.display = 'block';
}

function submitGeneInfo(inputValue, selectValue) {
    let data = new FormData();
    data.append('gene', inputValue);
    data.append('type', selectValue);

    return fetch(getGeneInfo, {
        method: 'POST',
        headers: {
            'X-Requested-With': 'XMLHttpRequest',
            'X-CSRFToken': csrfToken
        },
        body: data
    })
        .then(response => response.json())
        .then(data => {
            if (!data.error) {
                return data.symbol;  // Return the symbol here
            } else {
                console.error('Error fetching gene info:', data.message);
                alert('Error fetching gene info: ' + data.message);
                throw new Error(data.message);
            }
        })
        .catch(error => {
            console.error('Error:', error);
            alert('Failed to fetch gene info: ' + error);
            throw error;
        });
}



function clearSession() {
    fetch('/clear-session/', {
        method: 'POST',
        headers: {
            'X-CSRFToken': csrfToken,
            'Content-Type': 'application/json',
        },
    })
    .then(response => {
        if (!response.ok) {
            return response.json().then(errorData => {
                console.error('Error in session clearing:', errorData);
                throw new Error('Session clearing failed');
            });
        }
        return response.json();
    })
    .then(deletionData => {
        console.log('Session cleared successfully:', deletionData);
    })
    .catch(error => {
        console.error('Fetch error during session clearing:', error);
    });
}
