document.addEventListener('DOMContentLoaded', function () {


    if (sessionStorage.getItem('userID') !== null) {
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
                    console.log('Session user set successfully:', data.message);
                } else {
                    var errorMessageContainer = 'Error in setting session user: ' + data.message;
                    showErrorModal(errorMessageContainer);   
                }            })
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

        fetch(getReaction + reactionId)
            .then(response => {
                if (!response.ok) {
                    console.error(`HTTP error! status: ${response.status}`);
                    throw new Error(`HTTP error! status: ${response.status}`);
                }
                return response.json();
            })
            .then(async reactionData => {
                await updateFormFields(reactionData);
                confirmAll();
                displayDivs(reactionData);
                setLoggedInStatusBasedOnUrl();
                DisplayTag(reactionData.Organs);
                })
            .catch(error => {
                console.error('Error fetching reaction data:', error);
                // Optionally, handle the error by displaying a message to the user
            });
    }
    if (subsystemList.length === 0) {
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





//TODO : Add the function to display the placeholders

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

