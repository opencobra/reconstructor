function hideChemdoodlestatus(show) {
    // Toggles the visibility of the ChemDoodle sketcher modal based on the 'show' parameter.
    var modal = document.getElementById('chemdoodleModal');
    if (show) {
        modal.style.display = 'block';
    } else {
        modal.style.display = 'none';
    }
}

function getIdentifierCell(selectElement) {
    return selectElement.closest('.inputs-group')?.querySelector('.cell-identifier');
}

function getIdentifierInput(selectElement) {
    return getIdentifierCell(selectElement)?.querySelector('input[type="text"]');
}

function attachEventListenersToSelects() {
    // Attaches 'change' event listeners to substrate and product type select elements.
    const substrateSelects = document.querySelectorAll('#substratesDiv select[name="substrates_type"]');
    const productSelects = document.querySelectorAll('#productsDiv select[name="products_type"]');

    substrateSelects.forEach(select => {
        select.addEventListener('change', handleSelectChange);
    });

    productSelects.forEach(select => {
        select.addEventListener('change', handleSelectChange);
    });
}

function handleSelectChange(event) {
    // Handles changes in substrate and product type select elements, adding/removing drawing buttons accordingly.
    let selectElement = event.target;
    if (selectElement.value === 'Draw') {
        addStartDrawingButton(selectElement);
        clearDrawingSketcher();
    } else {
        removeStartDrawingButton(selectElement);
        removeEditDrawingButton(selectElement);
        const textField = getIdentifierInput(selectElement);
        if (textField) {
            textField.style.visibility = 'visible';
        }
    }
}
function addStartDrawingButton(selectElement) {
    // Hide the input field and replace it with the 'Start Drawing' button
    let textField = getIdentifierInput(selectElement);
    if (textField) {
        textField.style.display = 'none'; // Hide the input field
        if (!textField.nextElementSibling || textField.nextElementSibling.className !== 'start-drawing-btn') {
            let startDrawingBtn = document.createElement('button');
            startDrawingBtn.type = 'button';
            startDrawingBtn.className = 'start-drawing-btn';
            startDrawingBtn.textContent = 'Start Drawing';
            startDrawingBtn.onclick = () => setCurrentlyDrawing(selectElement);
            textField.parentNode.insertBefore(startDrawingBtn, textField); // Insert the button in place of the input field
        }
    }
}

function removeStartDrawingButton(selectElement) {
    // Removes the 'Start Drawing' button adjacent to the specified select element, if it exists.
    let startDrawingButton = getIdentifierCell(selectElement)?.querySelector('.start-drawing-btn');
    if (startDrawingButton) {
        startDrawingButton.remove();
    }
}
function removeEditDrawingButton(selectElement) {

    let editDrawingButton = getIdentifierCell(selectElement)?.querySelector('.edit-drawing-btn');

    if (editDrawingButton) {
        editDrawingButton.remove();
    }

}
function setCurrentlyDrawing(selectElement) {
    // Marks a select element as 'currentlyDrawing' and opens the ChemDoodle sketcher modal.
    clearCurrentlyDrawing();
    selectElement.setAttribute('currentlyDrawing', true);
    hideChemdoodlestatus(true);
}

function clearCurrentlyDrawing() {
    // Clears the 'currentlyDrawing' attribute from all select elements.
    document.querySelectorAll('select[currentlyDrawing="true"]').forEach(select => {
        select.removeAttribute('currentlyDrawing');
    });
}
function saveDrawing() {
    // Close the ChemDoodle sketcher modal
    hideChemdoodlestatus(false);

    // Get the molecule data from the ChemDoodle iframe
    let iframe = document.getElementById('chemdoodleIframe');
    let molFile = iframe.contentWindow.getMoleculeData();

    // Clean up the molecule data string
    molFile = molFile.replace("from ChemDoodle Web Components", 'Drawn Molecule');
    molFile = molFile.replace("http://www.ichemlabs.com", " ");

    // Encode the molecule data for storage in the input field
    let encodedMolFile = encodeURIComponent(molFile);

    // Find the currently active select element
    let currentlyDrawingSelect = document.querySelector('select[currentlyDrawing="true"]');
    if (currentlyDrawingSelect) {
        let textField = getIdentifierInput(currentlyDrawingSelect);
        if (textField) {
            textField.value = encodedMolFile; // Save the encoded molecule data in the input field
            textField.dispatchEvent(new Event('input', { bubbles: true }));
            textField.dispatchEvent(new Event('change', { bubbles: true }));
        }

        // Remove the "Start Drawing" button and replace it with the "Edit Drawing" button
        removeStartDrawingButton(currentlyDrawingSelect);
        addEditDrawingButton(currentlyDrawingSelect);
    }

    // Clear the 'currentlyDrawing' attribute from all select elements
    clearCurrentlyDrawing();
}


function addEditDrawingButton(selectElement) {
    // Check if the input field exists
    let textField = getIdentifierInput(selectElement);
    
    if (textField) {
        // Hide the input field
        textField.style.display = 'none'; // Hide the input field

        // Ensure that an 'Edit Drawing' button is not already added next to the input field
        if (!textField.nextElementSibling || textField.nextElementSibling.className !== 'edit-drawing-btn') {
            let editDrawingBtn = document.createElement('button');
            editDrawingBtn.type = 'button';
            editDrawingBtn.className = 'edit-drawing-btn';
            editDrawingBtn.textContent = 'Edit Drawing';
            editDrawingBtn.onclick = editDrawing; // Assign the editDrawing function

            // Insert the button in place of the input field
            textField.parentNode.insertBefore(editDrawingBtn, textField);
        }
    } else {
        // If no input field is found, insert the button after the select element
        const cell = getIdentifierCell(selectElement);
        if (cell && (!cell.querySelector('.edit-drawing-btn'))) {
            let editDrawingBtn = document.createElement('button');
            editDrawingBtn.type = 'button';
            editDrawingBtn.className = 'edit-drawing-btn';
            editDrawingBtn.textContent = 'Edit Drawing';
            editDrawingBtn.onclick = editDrawing; // Assign the editDrawing function

            // Insert the button after the select element
            cell.insertBefore(editDrawingBtn, cell.firstChild);
            const identifierInput = getIdentifierInput(selectElement);
            if (identifierInput) {
                identifierInput.style.visibility = 'hidden';
            }
        }
    }
}

function editDrawing() {
    // Prepares the specified select element for editing the drawing and loads existing molecule data into the ChemDoodle sketcher.
    let selectElement = this.previousElementSibling; // Assuming the select element is right before the button
    setCurrentlyDrawing(selectElement);

    // Retrieve the molecule data from the associated text field
    let textField = getIdentifierInput(selectElement);
    let molFile = decodeURIComponent(textField.value);
    // Send the molecule data to the ChemDoodle sketcher
    let iframe = document.getElementById('chemdoodleIframe');
    iframe.contentWindow.loadMoleculeIntoSketcher(molFile);
    // Open the ChemDoodle modal
    hideChemdoodlestatus(true);
}

function clearDrawingSketcher() {
    // Clears the ChemDoodle sketcher modal and resets the currentlyDrawing attribute.
    let iframe = document.getElementById('chemdoodleIframe');
    var emptyMolFile = "Molecule%20from%20ChemDoodle%20Web%20Components%0A%0Ahttp%3A%2F%2Fwww.ichemlabs.com%0A%20%201%20%200%20%200%20%200%20%200%20%200%20%20%20%20%20%20%20%20%20%20%20%20999%20V2000%0A%20%20%20%200.0000%20%20%20%200.0000%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%0AM%20%20END";
    let molFile = decodeURIComponent(emptyMolFile);
    iframe.contentWindow.loadMoleculeIntoSketcher(molFile);
}

