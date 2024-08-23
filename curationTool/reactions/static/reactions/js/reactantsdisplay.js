document.addEventListener('DOMContentLoaded', function () {


    if (sessionStorage.getItem('userID') !== null) {
        var errorMessageContainer = document.getElementById('error-message');

        setLoggedInStatusBasedOnUrl();
        username = sessionStorage.getItem('userName');
        userID = sessionStorage.getItem('userID');
        document.getElementById('userDisplay').innerHTML = `<i class="icon user"></i> User: ${username}`;
        document.getElementById('loginButton').textContent = 'Log out';

        fetch(setSessionUser, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ 'userID': userID })
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    errorMessageContainer.style.display = 'none';
                    errorMessageContainer.innerHTML = '';

                } else {
                    errorMessageContainer.innerHTML = 'Error in setting session user: ' + data.message;
                    window.scrollTo(0, 0);
                }
            })
    }


    setupTooltips();
    createnewreaction();
    createGeneInfoInput();
    displayreactioninfo(reactionData = null);

    attachEventListenersToSelects();
    toggleStructure();


    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    const action = urlParams.get('action');

    if (reactionId || action === 'edit') {
        console.log('Starting fetch request for reactionId:', reactionId);

        fetch(getReaction + reactionId)
            .then(response => {
                console.log('Fetch request completed. Checking response status...', getReaction + reactionId);
                if (!response.ok) {
                    console.error(`HTTP error! status: ${response.status}`);
                    throw new Error(`HTTP error! status: ${response.status}`);
                }
                console.log('Response is OK. Parsing JSON...');
                return response.json();
            })
            .then(async reactionData => {
                console.log('Parsed JSON data:', reactionData);
                await updateFormFields(reactionData);
                confirmAll();
                displayDivs(reactionData);
                setLoggedInStatusBasedOnUrl();
                console.log('Data id in reactants:', reactionData);
            })
            .catch(error => {
                console.error('Error fetching reaction data:', error);
                // Optionally, handle the error by displaying a message to the user
            });
    }
    if (subsystemList.length === 0) {
        console.log('Fetching subsystems from VMH...');
        window.scrollTo(0, 0);
        fetch(getVMHsubsystems, {
            method: 'GET',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            }
        })
        .then(response => response.json())
        .then(data => {
            if (data.error) {
                showErrorModal(data.message);
            }
            else {
                subsystemList = data.subsystem_list; // Save the returned list
                hidemodal();
            }
        })
        .catch(error => console.error('Error:', error));
    }
    setupTooltips();

});


function displayreactioninfo(reactionData) {
    // Extract the reactionId from reactionData
    const reactionId = reactionData;
    var user_id = sessionStorage.getItem('userID');
    // Log the reactionId for debugging purposes
    console.log('Data id in reactants:', reactionId);

    if (reactionId === null) {
        // If reactionId is null, fetch session data
        fetch('/check-session/', {
            method: 'GET',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken, // Include CSRF token if needed
            }
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    // If session data is found, display the gene info
                    console.log('Gene Info from session:', data.gene_info);
                    displayTabContent('gene-info-content', data.gene_info, null);
                } else {
                    console.error('Error:', data.message);
                }
            })
            .catch(error => {
                console.error('Error fetching session data:', error);
            });
    } else {
        // If reactionId is not null, first fetch session data
        fetch('/check-session/', {
            method: 'GET',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken, // Include CSRF token if needed
            }
        })
        .then(sessionResponse => sessionResponse.json())
        .then(sessionData => {
            console.log('Session data:', sessionData);
        
            if (sessionData.status === 'success' && sessionData.gene_info) {
                // Process each item in the session data
                sessionData.gene_info.forEach(geneInfoItem => {
                    const dataToSubmit = {
                        userID: user_id,  // assuming you have this from your context
                        infoType: 'Gene Info',
                        infoText: geneInfoItem.info,
                        reactionId: reactionId.reaction_id  // Include the reaction ID
                    };
                    console.log('Data to submit:', dataToSubmit);
                    // Submit the data to the reaction
                    submitData(dataToSubmit);
                });
        
                // After all data has been submitted, clear the session
                clearSession();
            } else {
                console.log('No gene info in session or an error occurred.');
            }
        })
        .catch(error => {
            console.error('An error occurred:', error);
            })
            .catch(sessionError => {
                console.error('Error handling session data:', sessionError);
            })
            .then(() => {
                // After checking the session, proceed with fetching reaction details
                fetch(getReactionDetails, {
                    method: 'POST',
                    body: JSON.stringify(reactionId.reaction_id),  // Assuming the POST requires a JSON body
                    headers: {
                        'X-Requested-With': 'XMLHttpRequest',
                        'X-CSRFToken': csrfToken,
                        'Content-Type': 'application/json'  // Ensuring the content type is JSON
                    }
                })
                    .then(response => response.json())
                    .then(details => {
                        console.log('Details: func', details.comments);

                        // Display content in the appropriate tabs
                        displayTabContent('refs-content', details.references, reactionId.reaction_id);
                        displayTabContent('ext-links-content', details.external_links, reactionId.reaction_id);
                        displayTabContent('gene-info-content', details.gene_info, reactionId.reaction_id);
                        displayTabContent('comments-content', details.comments, reactionId.reaction_id);
                    })
                    .catch(error => {
                        console.error('Error fetching reaction details:', error);
                    });
            });
    }
}


function displayDivs(reactionData) {
    loadAtomMappingDiv(reactionData);
    loadChemInfoDiv(reactionData);
    loadMetaboliteInfoDiv(reactionData);
    refreshSideButtons();
    updateStatusDots('substratesDiv', reactionData.subs_found, reactionData.subs_miriams);
    updateStatusDots('productsDiv', reactionData.prod_found, reactionData.prod_miriams);
    displayReactionMessage(reactionData);
    displayreactioninfo(reactionData);
}


function updateStatusDots(containerId, foundList, miriamsList) {
    const container = document.getElementById(containerId);
    const statusDots = container.querySelectorAll('.status-dot');
    statusDots.forEach((dot, index) => {
        dot.style.display = 'block'; // Ensure the dot is displayed

        if (foundList[index]) {
            dot.className = 'status-dot found';
            dot.setAttribute('data-tooltip', 'Metabolite found in VMH');
            dot.style.backgroundColor = ''; // Reset color to default
        } else {
            dot.className = 'status-dot not-found';
            dot.setAttribute('data-tooltip', 'Metabolite not found in VMH');
            dot.style.backgroundColor = ''; // Reset color to default
        }

        if (foundList[index] && miriamsList[index]) {
            dot.onclick = () => window.open(miriamsList[index], '_blank');
            dot.style.cursor = 'pointer';
        } else {
            dot.onclick = null;
            dot.style.cursor = 'default';
        }
    });
}

function updateStatusDot(dot, found, miriam) {
    dot.style.display = 'block'; // Ensure the dot is displayed

    if (found) {
        dot.className = 'status-dot found';
        dot.setAttribute('data-tooltip', 'Metabolite found in VMH');
        dot.style.backgroundColor = ''; // Reset color to default
    } else {
        dot.className = 'status-dot not-found';
        dot.setAttribute('data-tooltip', 'Metabolite not found in VMH');
        dot.style.backgroundColor = ''; // Reset color to default
    }

    if (found && miriam) {
        dot.onclick = () => window.open(miriam, '_blank');
        dot.style.cursor = 'pointer';
    } else {
        dot.onclick = null;
        dot.style.cursor = 'default';
    }
}

function confirmAll() {
    // Get all the elements with the class 'done-field-btn-all'
    const verifyAllButtons = document.querySelectorAll('.done-field-btn-all');

    // Loop through each button and trigger a click event
    verifyAllButtons.forEach(button => button.click());
}







function setupTooltips() {
    document.querySelectorAll('.info-symbol').forEach(item => {
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


// document.getElementById("viewsavedreactions-button").addEventListener("click", function() {
//     fetch('/saved_reactions/') // URL endpoint to fetch saved reactions
//       .then(response => response.text())
//       .then(data => {
//         document.getElementById("modalContent-saverxn").innerHTML = data;
//         $('#savedReactionsModal').modal('show'); // Using Semantic UI modal
//       })
//       .catch(error => console.error('Error:', error));
//   });



function createnewreaction() {
    // Get the input element by its ID
    const resetbutton = document.getElementById('ResetButton');


    // Add an event listener to the input element to handle the click event
    resetbutton.addEventListener('click', function (event) {
        // Prevent the default form submission behavior
        event.preventDefault();

        // Redirect to the homepage
        window.location.href = window.location.origin;
    });
}




function toggleStructure() {
    const buttons = document.querySelectorAll('.toggle-button');
    buttons.forEach(button => {
        button.addEventListener('click', () => {
            const messageContainer = button.parentNode.nextElementSibling;
            if (messageContainer.style.display === 'none' || messageContainer.style.display === '') {
                messageContainer.style.display = 'block';
                button.textContent = 'Hide Message';
            } else {
                messageContainer.style.display = 'none';
                button.textContent = 'Show Message';
            }
        });
    });
}




function displayReactionMessage(reactionData) {
    // Clear the previous message
    let messageContainer = document.getElementById('reactionFoundMessage');

    messageContainer.innerHTML = '';

    if (reactionData.vmh_found) {
        const reactionFoundMessage = reactionData.vmh_found_similar ?
            'Similar Reaction found at VMH:' : 'Exact Reaction found at VMH:';

        const messageText = document.createTextNode(reactionFoundMessage);
        messageContainer.appendChild(messageText);
        messageContainer.appendChild(document.createElement('br'));

        const link = document.createElement('a');
        link.setAttribute('href', reactionData.vmh_url.trim().replace(/^"|"$/g, ''));
        link.setAttribute('target', '_blank');
        link.textContent = reactionData.vmh_url;
        messageContainer.appendChild(link);
    } else {
        const messageText = document.createTextNode('Reaction not found at VMH');
        messageContainer.appendChild(messageText);
    }

    // Display the message
    messageContainer.style.display = 'block';
}



function setLoggedInStatusBasedOnUrl() {
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');

    const action = urlParams.get('action');
    console.log(action);
    // 2. Determine the status for logged-in users
    let status = "Creating reaction"; // Default for logged-in users
    let dotClass = "dot-blue"; // Default color for "Creating reaction"

    if (action != "edit" && reactionId != null) {
        status = "Viewing reaction";
        dotClass = "dot-green"; // Green dot for "Viewing reaction"
    } else if (reactionId == null && action != "edit") {
        status = "creating reaction";
        dotClass = "dot-red"; // Orange dot for "Editing reaction"
    }
    else {
        status = "Editing reaction";
        dotClass = "dot-orange"; // Orange dot for "Editing reaction"
    }

    // 3. Display the status in the div with id="statusTitle"
    const statusElement = document.getElementById('statusTitle');
    if (statusElement) {
        statusElement.innerHTML = `Status: <span class="status-dot-top ${dotClass}"></span> ${status}`;
    } else {
        console.error("Status element not found in the DOM");
    }
}

function setLoggedOutStatusBasedOnUrl() {
    // 1. Determine the status for logged-out users
    const status = "Idle"; // Default for logged-out users
    const dotClass = "dot-grey"; // Grey dot for "Idle"

    // 2. Display the status in the div with id="statusTitle"
    const statusElement = document.getElementById('statusTitle');
    if (statusElement) {
        statusElement.innerHTML = `Status: <span class="status-dot-top ${dotClass}"></span> ${status}`;
    } else {
        console.error("Status element not found in the DOM");
    }
}



function createGeneInfoInput() {
    var container = document.createElement('div');
    container.id = 'geneInfoContainer';
    container.className = 'gene-info-container';

    // Create a container for gene input
    var geneInputContainer = document.createElement('div');
    geneInputContainer.className = 'gene-input-container';

    // Inputs container will support multiple rows for elements, but not appended yet
    var inputsContainer = document.createElement('div');
    inputsContainer.className = 'gene-inputs-container';

    // Create an initial row for inputs, this also appends an outer drop zone immediately after
    createInputRow(inputsContainer, '#000'); // Now directly modifies inputsContainer

    // Now that initial content is configured, append inputsContainer to the gene input container
    geneInputContainer.appendChild(inputsContainer);

    // Create a button container for AND, OR, and brackets
    var geneButtonContainer = document.createElement('div');
    geneButtonContainer.className = 'gene-button-container';
    ['Gene', 'AND', 'OR', '(', ')'].forEach(function (label) {
        var button = createButton(label);
        geneButtonContainer.appendChild(button);
    });

    // Append gene input container and button container to the main container
    container.appendChild(geneInputContainer);
    container.appendChild(geneButtonContainer);

    // Finally, append the container to the geneInfoWrapper
    document.getElementById('geneInfoWrapper').appendChild(container);
}

function createInputRow(inputsContainer, color = '#000') {
    var inputRow = document.createElement('div');
    inputRow.className = 'gene-input-row';

    var geneInput = document.createElement('div');
    geneInput.contentEditable = 'true';
    geneInput.className = 'gene-info-input';
    geneInput.id = 'geneInfoInput'; //+ Math.random();
    geneInput.style.border = '1px solid #ccc'; // Add some styling to make it look like an input field
    geneInput.style.padding = '5px';
    geneInput.style.minHeight = '20px'; // Ensure it has some height
    geneInput.style.display = 'inline-block';
    geneInput.style.width = '99.5%'; // Set the width to 100% to make it responsive
    geneInput.style.overflow = 'hidden'; // Hide overflow
    geneInput.style.boxSizing = 'border-box';
    geneInput.style.marginTop = '16px';
    geneInput.style.height = '39px';

    inputRow.appendChild(geneInput);

    // Append the input row to the inputs container
    inputsContainer.appendChild(inputRow);
}

// Function to create a button with an event listener
function createButton(label) {
    var button = document.createElement('button');
    button.textContent = label;
    button.className = 'gene-info-button';

    // Add event listener to the button
    button.addEventListener('click', function () {
        handleButtonClick(label);
    });

    // Inject specific styles for this button
    var style = document.createElement('style');
    style.textContent = `
        .gene-info-button {
            padding: 10px 15px;
            font-size: 13px;
            cursor: pointer;
            background-color: #007bff;
            color: white;
            border: none;
            border-radius: 4px;
            transition: background-color 0.3s;
        }
        .gene-info-button:hover {
            background-color: #0056b3;
        }
    `;
    document.head.appendChild(style); // Append the style to the document head

    return button;
}

document.querySelectorAll('button').forEach(button => {
    button.addEventListener('click', (event) => {
        handleButtonClick(event.target.textContent);
    });
});

function handleButtonClick(label) {
    if (label === 'AND' || label === 'OR') {
        var cursorPosition = getCursorPosition();
        insertTextAtCursor(" " + label + " ", cursorPosition);
    } else if (label === '(' || label === ')') {
        var cursorPosition = getCursorPosition();
        insertTextAtCursor(label, cursorPosition);
    } else if (label === 'Gene') {
        var cursorPosition = getCursorPosition();
        createAndShowModal((inputValue, typeValue) => {
            return submitGeneInfo(inputValue, typeValue)
                .then(symbol => {
                    // Store the returned symbol in a variable
                    const returnedSymbol = symbol;
                    insertTextAtCursor(returnedSymbol, cursorPosition);
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


document.getElementById('saveReactionButton').addEventListener('click', async function () {
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    var userId = sessionStorage.getItem('userID');

    // Check if reactionId is empty
    if (!reactionId) {
        alert('Reaction not created.');
        return;
    }

    // Construct the URL for the AJAX request
    var url = `/check-reaction`;

    try {
        let response = await fetch(url, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({
                user_id: userId,
                reaction_id: reactionId
            })
        });

        if (!response.ok) {
            alert('Please Log in to save the reaction.')
            return;
        }

        let data = await response.json();

        if (data.is_reaction_saved) {
            alert('Reaction is already saved.');
        } else {
            // Show the custom modal
            var modal = document.getElementById('saveReactionModal');
            modal.style.display = 'block';
            document.getElementById('modalBackground').style.display = 'block';

        }
    } catch (error) {
        console.error('Error checking reaction:', error);
        alert('An error occurred while checking the reaction.');
    }
});

// Event listener for closing the modal when the close button is clicked



document.getElementById('closeSaveReactionModal').addEventListener('click', function () {
    document.getElementById('saveReactionModal').style.display = 'none';
    document.getElementById('modalBackground').style.display = 'none';
});


document.getElementById('submitSaveReaction').addEventListener('click', function () {
    var userID = sessionStorage.getItem('userID');
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    var shortNameInput = document.getElementById('reactionNameInput');
    var shortName = shortNameInput.value;

    // Clear any previous custom validity message
    shortNameInput.setCustomValidity('');

    if (userID && reactionId) {
        if (!shortName) {
            // Alert if the short name is not provided
            alert('Please enter a short name for the reaction.');
            shortNameInput.setCustomValidity('Please enter a short name for the reaction.');
            shortNameInput.reportValidity();
            return; // Prevent form submission
        }

        var data = new FormData();
        data.append('userID', userID);
        data.append('reaction_id', reactionId);
        data.append('short_name', shortName);

        fetch(saveReaction, {  // Use the correct URL here
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken // Ensure csrfToken is correctly defined or fetched
            },
            body: data,
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    var modal = document.getElementById('saveReactionModal');
                    modal.style.display = 'none';
                    document.getElementById('modalBackground').style.display = 'none';

                    alert("Reaction saved successfully!");
                } else {
                    alert("Error: " + (data.message || "Failed to save the reaction."));
                }
            })
            .catch(error => console.error('Error:', error));
    } else if (!reactionId) {
        alert("Create the reaction first.");
    } else {
        alert("Please log in.");
    }

});




function enforceLoginRequirement() {
    // Define an array of selectors using both class and name attributes
    const selectors = [
        { class: 'content-div', name: 'references-div', header: 'References' },
        { class: 'content-div', name: 'extlinks-div', header: 'External Links' },
        { class: 'content-div', name: 'comments-div', header: 'Comments' },
        { class: 'content-div', name: 'gene_info-div', header: 'Gene Info' },
        // Add more selectors here if needed
    ];

    // Loop through each selector
    selectors.forEach(selector => {
        // Construct the CSS selector
        const query = `div.${selector.class}[name="${selector.name}"]`;

        // Select the div elements that match the query
        const targetDivs = document.querySelectorAll(query);

        // Update each targeted div
        targetDivs.forEach(div => {
            // Clear the existing content
            div.innerHTML = '';

            // Create and append the header
                const header = document.createElement('div');
                header.classList.add('div-header');
                header.textContent = selector.header;
                div.appendChild(header);

                // Create and style the message box
                const messageBox = document.createElement('div');
                messageBox.classList.add('login-message-box');

                // Create and append the message text
                const messageText = document.createElement('p');
                messageText.classList.add('login-message-text');
                messageText.textContent = 'Please log in to use this feature';

                messageBox.appendChild(messageText);
                div.appendChild(messageBox);
        });
    });
}

// Call the function to enforce the login requirement


// Function to clear the Django session
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


