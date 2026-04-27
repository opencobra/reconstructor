document.getElementById('createTemplateBtn').addEventListener('click', function() {
    // check if the user is logged in
    if (!sessionStorage.getItem('userID')) {
        alert('Please log in to create a template.');
        return;
    }
    // Hide the first modal
    $('#getrxnfromtemplateModal').modal('hide');

    // Enter "Create Template" mode
    toggleCreateTemplateMode();
});

function toggleCreateTemplateMode() {
    const statusElement = document.getElementById('statusTitle');
    const dotClass = 'dot-red';
    const saveButton = document.getElementById('saveReactionButton');
    let createTemplateButton = document.getElementById('createTemplateButton');
    const resetButton = document.getElementById('ResetButton');
    const formButtonsContainer = document.getElementById('FormButtons');

    // Update status
    statusElement.innerHTML = `Status: <span class="status-dot-top ${dotClass}"></span> Creating Template`;

    // Hide save button
    saveButton.style.display = 'none';
    // Hide create reaction button
    document.getElementById('submitBtn-form').style.display = 'none';
    // Hide skip-atom-mapping
    document.getElementsByClassName('skip-atom-mapping')[0].style.display = 'none';
    // Hide getbuttoncontainer class
    document.getElementsByClassName('getbuttoncontainer')[0].style.display = 'none';

    // Create and show 'Create Template' button if it doesn't exist
    if (!createTemplateButton) {
        createTemplateButton = document.createElement('button');
        createTemplateButton.id = 'createTemplateButton';
        createTemplateButton.className = 'ui button action-btn';
        createTemplateButton.textContent = 'Create Template';

        // Copy styles from submitBtn-form
        const submitButton = document.getElementById('submitBtn-form');
        if (submitButton) {
            createTemplateButton.style.color = submitButton.style.color;
            createTemplateButton.innerHTML = submitButton.innerHTML.replace(
                'Create Reaction', 
                'Create Template'
            );
        }

        // Place 'Create Template' button before the ResetButton, respecting new layout
        const insertionTarget = resetButton && resetButton.parentElement ? resetButton.parentElement : formButtonsContainer;
        insertionTarget.insertBefore(createTemplateButton, resetButton);

        // Add event listener for creating the template
        createTemplateButton.addEventListener('click', function() {
            // Collect form data and send it to the backend
            createTemplate();
        });
    } else {
        createTemplateButton.style.display = 'inline-flex';
    }
}

async function createTemplate() {
    // Show the modal for template creation
    $('#createTemplateModal').modal({
        onApprove: async function () {
            const templateName = document.getElementById('templateNameInput').value.trim();
            const templateDescription = document.getElementById('templateDescriptionInput').value.trim();

            // Validate template name
            if (!templateName) {
                alert('Template name is required.');
                return false; // Prevent modal from closing
            }

            // Collect form data
            const reactionForm = document.getElementById('reactionForm');
            const formData = new FormData(reactionForm);
            const directionElement = document.getElementById('reactionDirection');
            formData.set('direction', directionElement ? directionElement.value : 'forward');

            // Append template name and description
            formData.append('template_name', templateName);
            formData.append('description', templateDescription);

            // Collect organ tags
            const organTagsDiv = document.getElementById('organTags');
            const organTags = Array.from(organTagsDiv.getElementsByClassName('tag'))
                .map(tag => tag.firstChild.textContent.trim());
            formData.append('organs', JSON.stringify(organTags));

            // Validate subsystem field
            if (subsystemList.length === 0) {
                subsystemList = await updateSubsystems();
            }
            const subsystemField = document.getElementById('subsystemField').value.trim();
            if (subsystemField) {
                const isValidSubsystem = subsystemList.some(
                    subsystem => subsystem.toLowerCase() === subsystemField.toLowerCase()
                );
                if (!isValidSubsystem) {
                    const userConfirmed = confirm(`Are you sure you want to add a new subsystem "${subsystemField}"?`);
                    if (!userConfirmed) {
                        alert('The subsystem entered is not valid.');
                        return false; // Prevent modal from closing
                    } else if (!sessionStorage.getItem('userID')) {
                        alert('Please log in to add a new subsystem.');
                        return false; // Prevent modal from closing
                    } else {
                        subsystemList.push(subsystemField); // Add to local list
                        persistSubsystemList(subsystemField); // Save subsystem if logged in
                    }
                }
                formData.append('subsystem', subsystemField);
            } else {
                formData.append('subsystem', '');
            }

            // Add the user ID if logged in
            const userID = sessionStorage.getItem('userID');
            if (userID) {
                formData.append('userID', userID);
            }

            // Send the template creation request
            try {
                const response = await fetch('/create_template/', {
                    method: 'POST',
                    body: formData,
                    headers: {
                        'X-Requested-With': 'XMLHttpRequest',
                        'X-CSRFToken': document.querySelector('[name=csrfmiddlewaretoken]').value,
                    },
                });

                if (response.ok) {
                    const data = await response.json();
                    alert('Template created successfully!');
                    location.reload(); // Refresh the page
                    return true; // Allow modal to close
                } else {
                    const errorData = await response.json();
                    alert(`Error creating template: ${errorData.message}`);
                    return false; // Prevent modal from closing
                }
            } catch (error) {
                console.error('Error creating template:', error);
                alert('An error occurred while creating the template.');
                return false; // Prevent modal from closing
            }
        },
        onDeny: function () {
            // Clear the form fields if the user cancels
            document.getElementById('templateNameInput').value = '';
            document.getElementById('templateDescriptionInput').value = '';
            return true; // Allow modal to close
        },
    }).modal('show');
}