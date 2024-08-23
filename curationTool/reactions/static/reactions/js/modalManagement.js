const bracketColors = ['#191919', '#F0A3FF', '#0075DC', '#993F00', '#2BCE48', '#FF5005', '#FFCC99', '#FFA8BB', '#4DB6AC', '#81C784'];
document.getElementById('loginButton').addEventListener('click', function () {
    if (sessionStorage.getItem('userID') !== null) {
        sessionStorage.removeItem('userID');
        sessionStorage.removeItem('userName');
        document.getElementById('userDisplay').textContent = 'User: None';
        document.getElementById('loginButton').textContent = 'Log in';
        document.getElementById('loginModal').style.display = 'none';
    } else {
        // Display the modal
        document.getElementById('loginModal').style.display = 'block';
    }
});

document.getElementById('loginButton').addEventListener('click', function() {
    if (this.textContent !== 'Login') {
        window.location.href = window.location.origin;
    }
});


document.getElementById('register').addEventListener('click', function () {
    document.getElementById('registerModal').style.display = 'block';

});


document.querySelector('.close-button-login').addEventListener('click', function () {
    document.getElementById('loginModal').style.display = 'none';
});

document.querySelector('.close-button').addEventListener('click', function () {
    document.getElementById('registerModal').style.display = 'none';
});

document.getElementById('loginForm').addEventListener('submit', function (event) {
    event.preventDefault(); // Prevent the form from submitting the traditional way
    var userName = document.getElementById('username').value;
    var password = document.getElementById('password').value;

    if (userName && password) {
        // Prepare the data to be sent in the POST request
        const data = new FormData();
        data.append('username', userName);
        data.append('password', password);

        // Adjust the fetch call to use POST
        fetch(`${getUser}`, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
            },
            body: data
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    sessionStorage.setItem('userID', data.userID);
                    sessionStorage.setItem('userName', data.userName);
                    document.getElementById('userDisplay').textContent = `User: ${data.userName}`;
                    document.getElementById('loginButton').textContent = 'Log out';
                    // Hide the modal after successful login
                    
                    document.getElementById('loginModal').style.display = 'none';
                    fetchAvailableReactions(data.userID);

                } else {
                    alert(data.message);
                }
            })
            .catch(error => {
                console.error('Error:', error);
                alert("Failed to fetch user details.");
            });
    }
});

async function fetchAvailableReactions(userId) {
    try {
        const response = await fetch('available_reactions', {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ user_id: userId })
        });

        if (!response.ok) {
            throw new Error('Network response was not ok');
        }

        const data = await response.json();

        if (data.error) {
            console.error(data.error);
            return;
        }

        handleLoginResponse(data);
    } catch (error) {
        console.error('There was a problem with the fetch operation:', error);
    }
}

function handleLoginResponse(data) {
    const loginModal = document.getElementById('loginModal');
    
    if (!data.available_reaction_ids || data.available_reaction_ids.length === 0) {
        loginModal.style.display = 'none';
    } else {
        const lastReactionId = data.last_index;
        window.location.href = window.location.origin + "/?reaction_id=" + lastReactionId;
    }
}


document.getElementById('registerForm').addEventListener('submit', function (event) {
    event.preventDefault(); // Prevent the form from submitting the traditional way
    var userName = document.getElementById('username-reg').value;
    var password = document.getElementById('password-reg').value;
    var orchid_id = document.getElementById('orchidid').value;
    var email = document.getElementById('email').value;

    if (userName && password && orchid_id && email) {
        // Prepare the data to be sent in the POST request
        const data = new FormData();
        data.append('username', userName);
        data.append('password', password);
        data.append('orchid_id', orchid_id);
        data.append('email', email);

        // Adjust the fetch call to use POST
        fetch(`${regUser}`, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
            },
            body: data
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    sessionStorage.setItem('userID', data.userID);
                    sessionStorage.setItem('userName', data.userName);
                    document.getElementById('userDisplay').textContent = `User: ${data.userName}`;
                    document.getElementById('loginButton').textContent = 'Log out';
                    // Hide the modal after successful login
                    document.getElementById('registerModal').style.display = 'none';
                    
                } else {
                    alert(data.message);
                }
            })
            .catch(error => {
                console.error('Error:', error);
                alert("Failed to fetch user details.");
            });
    }
});
// Opens the 'Add Information' modal and clears any existing text and response messages.
document.getElementById('addInfoButton').addEventListener('click', function () {
    userID = sessionStorage.getItem('userID');
    if (!userID) {
        alert("Please log in.");
        return;
    }
    document.getElementById('addInfoModal').style.display = 'block';
    document.getElementById('infoTextInput').textContent = '';
    document.getElementById('responseMessage').textContent = '';
});
// Closes the 'Add Information' modal when the close button is clicked.
document.getElementById('closeAddInfoModal').addEventListener('click', function () {
    document.getElementById('addInfoModal').style.display = 'none';
});

// Sets up dynamic placeholders and dropdowns for different types of information in the 'Add Information' modal.
document.addEventListener('DOMContentLoaded', function () {
    var infoTypeSelect = document.getElementById('infoTypeSelect');
    var infoTextInput = document.getElementById('infoTextInput');
    function updatePlaceholderAndDropdown() {
        var selectedOption = infoTypeSelect.value;
        switch (selectedOption) {
            case 'Reference':
                restoreInfoTextInput();
                removeAdditionalElements();
                createAdditionalDropdown(['DOI', 'PMID'], 'referenceTypeSelect');
                break;
            case 'External Link':
                restoreInfoTextInput();
                removeAdditionalElements();
                createAdditionalDropdown(['CHO Models', 'COG', 'EC Number', 'KEGG orthology', 'KEGG reaction', 'MetanetX', 'Rhea', 'SEED', 'Wikipedia'], 'externalLinkSelect');
                break;
            case 'Gene Info':
                createModal();
                hideInfoTextInput();
                removeAdditionalElements();
                createGeneInfoInput();
                disableKeyboardInput();
                focusGeneInput();
                break;
            case 'Comment':
                restoreInfoTextInput();
                infoTextInput.placeholder = 'Enter your comment here';
                removeAdditionalElements();
                break;
        }
    }

    function hideInfoTextInput() {
        infoTextInput.style.display = 'none'; // Hide the input
    }
    
    function focusGeneInput() {
        var geneInput = document.getElementById('geneInfoInput');
        geneInput.focus();
    }


    function createModal() {
        // Create the CSS for the modal
        const style = document.createElement('style');
        style.innerHTML = `
            body {
                font-family: Arial, sans-serif;
            }
    
            .modal {
                display: none; /* Hidden by default */
                position: fixed; /* Stay in place */
                z-index: 1000; /* Sit on top */
                left: 0;
                top: 0;
                width: 100%; /* Full width */
                height: 100%; /* Full height */
                overflow: auto; /* Enable scroll if needed */
                background-color: rgb(0,0,0); /* Fallback color */
                background-color: rgba(0,0,0,0.4); /* Black w/ opacity */
                padding-top: 60px;
            }
    
            .modal-content {
                background-color: #fefefe;
                margin: 15% auto; /* 15% from the top and centered */
                padding: 20px;
                border: 1px solid #888;
                width: 80%; /* Could be more or less, depending on screen size */
            }
    
            .closeBtn {
                color: #aaa;
                float: right;
                font-size: 28px;
                font-weight: bold;
            }
    
            .closeBtn:hover,
            .closeBtn:focus {
                color: black;
                text-decoration: none;
                cursor: pointer;
            }
        `;
        document.head.appendChild(style);
    
        // Create the modal HTML structure
        const modal = document.createElement('div');
        modal.id = 'myModal';
        modal.className = 'modal';
        modal.innerHTML = `
        <div class="modal-content">
        <span class="closeBtn">&times;</span>
        <h2>How to Add GPR</h2>
        <ol>
            <li>Click on the <strong>Gene</strong> button to add Gene in input field.</li>
            <li>Select whether you want to add a Entrez ID or an HGNC symbol.</li>
            <li>Add operations like <strong>AND</strong>, <strong>OR</strong>, or parentheses <strong>()</strong> by clicking on the specific button.</li>
            <li><strong>Ensure your cursor is placed in the desired input field before entering data. If the cursor is not in the input field, the value will not be seen in input field.</strong></li>
        </ol>
    </div>
    
        `;
        document.body.appendChild(modal);
    
        // Get modal element
        const modalElement = document.getElementById("myModal");
        // Get close button
        const closeBtn = document.getElementsByClassName("closeBtn")[0];
    
        // Function to open modal
        function openModal() {
            modalElement.style.display = "block";
        }
    
        // Function to close modal
        function closeModal() {
            modalElement.style.display = "none";
        }
    
        // Function to close modal if outside click
        function outsideClick(e) {
            if (e.target == modalElement) {
                modalElement.style.display = "none";
            }
        }
    
        // Listen for close click
        closeBtn.addEventListener("click", closeModal);
        // Listen for outside click
        window.addEventListener("click", outsideClick);
    
        // Show the modal immediately
        openModal();
    }

    function restoreInfoTextInput() {
        infoTextInput.style.display = ''; // Restore the input display
    }

    var bracketColors = ['#000']; // Assuming this variable exists for color options

    function createGeneInfoInput() {
        var infoTextInput = document.getElementById('infoTextInput');
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
        createInputRow(inputsContainer, bracketColors[0]); // Now directly modifies inputsContainer

        // Now that initial content is configured, append inputsContainer to the gene input container
        geneInputContainer.appendChild(inputsContainer);

        // Create a button container for AND, OR, and brackets
        var geneButtonContainer = document.createElement('div');
        geneButtonContainer.className = 'gene-button-container';
        ['Gene','AND', 'OR', '(', ')'].forEach(function (label) {
            var button = createButton(label);
            geneButtonContainer.appendChild(button);
        });

        // Append gene input container and button container to the main container
        container.appendChild(geneInputContainer);
        container.appendChild(geneButtonContainer);

        // Finally, insert the container before infoTextInput in the DOM
        infoTextInput.parentNode.insertBefore(container, infoTextInput);
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
        geneInput.style.width = '100%'; // Set the width to 100% to make it responsive
        geneInput.style.maxWidth = '500px'; // Set a max width
        geneInput.style.overflow = 'hidden'; // Hide overflow
        geneInput.style.boxSizing = 'border-box';
        geneInput.style.marginTop = '16px';
        geneInput.style.height = '39px';

        inputRow.appendChild(geneInput);

        // Append the input row to the inputs container
        inputsContainer.appendChild(inputRow);

    }

// new function i added 
// Function to create a button with an event listener
function createButton(label) {
    var button = document.createElement('button');
    button.textContent = label;
    button.className = 'gene-info-button';

    // Add event listener to the button
    button.addEventListener('click', function () {

        handleButtonClick(label);

    });

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
    inputElement.addEventListener('keypress', function(event) {
        if (!shouldAllowKey(event)) {
            event.preventDefault();
        }
    });

    inputElement.addEventListener('keydown', function(event) {
        if (!shouldAllowKey(event)) {
            event.preventDefault();
        }
    });

    inputElement.addEventListener('keyup', function(event) {
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
    closeButton.onclick = function() {
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

    dropdown.onchange = function() {
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

    confirmButton.onclick = function() {
        const selectedType = dropdown.value;
        const inputValue = geneInputModal.value.trim();
        const typeValue = selectedType === 'Entrez ID' ? 'Entrez ID' : 'HGNC Symbol';

        if (!inputValue) {
            console.error("Input value is undefined or empty");
            alert("Please enter a gene value");
            return;
        }

        // Call the onConfirm callback with the input values and handle the returned promise
        onConfirm(inputValue, typeValue, 'geneInfoInput', confirmButton.id)
            // .then(symbol => {
            //     // Handle the symbol returned from the promise
            //     console.log('Returned Symbol:', symbol);
            //     // Store the symbol in a variable or perform further actions
            //     const returnedSymbol = symbol;
            //     // Hide the modal after confirming
            //     modal.style.display = 'none';
            // })
            // .catch(error => {
            //     console.error('Error:', error);
            // });
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


function createAdditionalDropdown(options, id) {
    removeAdditionalElements(); // Ensure any existing dropdown is removed first

    var additionalSelect = document.createElement('select');
    additionalSelect.id = id;
    additionalSelect.className = 'add-info-dropdown';

    options.forEach(function (option) {
        var opt = document.createElement('option');
        opt.value = option;
        opt.innerHTML = option;
        additionalSelect.appendChild(opt);
    });

    infoTextInput.parentNode.insertBefore(additionalSelect, infoTextInput);

    // Set initial placeholder based on the first option
    updatePlaceholder(additionalSelect.options[additionalSelect.selectedIndex].value);

    additionalSelect.addEventListener('change', function () {
        updatePlaceholder(additionalSelect.options[additionalSelect.selectedIndex].value);
    });
}

function updatePlaceholder(selectedOption) {
    if (selectedOption === 'PMID') {
        infoTextInput.placeholder = "PMID Identifier(s) can add multiple PMID separated by ;";
    } else if (selectedOption === 'DOI') {
        infoTextInput.placeholder = "DOI Identifier(s) can add multiple DOI separated by ,  ";
    } else {
        infoTextInput.placeholder = selectedOption + " Identifier";
    }
}


    function removeAdditionalElements() {
        // Remove additional dropdowns
        var existingSelects = ['externalLinkSelect', 'referenceTypeSelect'];
        existingSelects.forEach(function (selectId) {
            var existingSelect = document.getElementById(selectId);
            if (existingSelect) {
                existingSelect.parentNode.removeChild(existingSelect);
            }
        });

        // Remove gene info specific elements
        var geneInputContainer = document.getElementsByClassName('gene-info-container');
        if (geneInputContainer.length > 0) {
            let arrayGeneInputContainer = Array.from(geneInputContainer);
            arrayGeneInputContainer.forEach(function (container) {
                container.parentNode.removeChild(container);
            });
        }
    }

    infoTypeSelect.addEventListener('change', updatePlaceholderAndDropdown);
    updatePlaceholderAndDropdown();
});
// Displays the 'Save Reaction' modal when the button is clicked.


// Event listener for the save reaction button
document.getElementById('saveReactionButton').addEventListener('click', async function () {
    var reactionId = document.getElementById('reactionIdInput').value;
    var userId = sessionStorage.getItem('userID');

    // Check if reactionId is empty
    if (!reactionId) {
        alert('Reaction ID cannot be empty.');
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
            throw new Error(`Error: ${response.status}`);
        }

        let data = await response.json();

        if (data.is_reaction_saved) {
            alert('Reaction is already saved.');
        } else {

            var modal = document.getElementById('saveReactionModal');
            modal.style.display = 'block';
        }
    } catch (error) {
        console.error('Error checking reaction:', error);
        alert('An error occurred while checking the reaction.');
    }
});

// Closes the 'Save Reaction' modal when the close button or overlay is clicked.
document.getElementById('closeSaveReactionModal').addEventListener('click', function () {
    var modal = document.getElementById('saveReactionModal');
    modal.style.display = 'none';
});

// Closes the 'Save Reaction' modal if a click occurs outside the modal content area.
window.addEventListener('click', function (event) {
    var modal = document.getElementById('saveReactionModal');
    if (event.target == modal) {
        modal.style.display = 'none';
    }
});




