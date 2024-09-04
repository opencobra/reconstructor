document.getElementById('submitBtn-form').addEventListener('click', function(event) {
    event.preventDefault(); // Prevent the default form submission

    // Trigger the form's submit event
    document.getElementById('reactionForm').requestSubmit();
});

function hidemodal(){

    document.getElementById('error-modal').style.display = 'none';

}

document.getElementById('close-button').onclick = function() {
    document.getElementById('error-modal').style.display = 'none';
};


document.getElementById('reactionForm').addEventListener('submit', function(e) {
    var loadingIndicator = document.getElementById('loadingIndicator');


    e.preventDefault(); // Prevent the default form submission

    var submitBtn = document.getElementById('submitBtn-form');
    currentUrl = window.location.href;
    if (currentUrl.includes('edit')) {
        // pop up that says "are you sure you want to update this reaction?"
        var userConfirmed = confirm('Are you sure you want to update this reaction?');
        if (!userConfirmed) {
            return; // Exit the function and do not submit form
        }
    }

    var loadingIndicator = document.getElementById('loadingIndicator');



    var divElement = document.getElementById("organTags");
    var organTags = Array.from(divElement.getElementsByClassName('tag'))
                            .map(tag => tag.firstChild.textContent.trim());


    // Check if the subsystem field is filled
    var subsystemField = document.getElementById('subsystemField').value;
    if (!subsystemField.trim()) {
        alert('Please enter a subsystem.');
        return; // Exit the function and do not submit form
    }

    var inputsGroups = document.querySelectorAll('.inputs-group');
    inputsGroups.forEach(function(group) {
        var statusDot = group.querySelector('.status-dot');
        if (!statusDot){
            console.error('Status dot not found in inputs group:', group);
            alert('Please confirm all metabolites before submitting.');
            return;
        }
    });


    var isValidSubsystem = subsystemList.some(subsystem => subsystem.toLowerCase() === subsystemField.toLowerCase());

    if (!isValidSubsystem) {
        var userConfirmed = confirm(`Are you sure you want to add a new subsystem "${subsystemField}"?`);
        if (!userConfirmed) {
            var errorMessage = 'The subsystem entered is not valid.';
            showErrorModal(errorMessage);
            window.scrollTo(0, 0);
            return; // Exit the function and do not submit form
        } else {
            if (sessionStorage.getItem('userID') !== null) {
            // Add the new subsystem to the list
            subsystemList.push(subsystemField);
            persistSubsystemList(subsystemField); // Persist the new subsystem
            }
            else {
                var errorMessage = 'Please login to add a new subsystem.';
                showErrorModal(errorMessage);
                window.scrollTo(0, 0);
                return; // Exit the function and do not submit form
            }
        }
    }

    // Continue with form submission

    const disabledInputs = this.querySelectorAll('input:disabled, select:disabled');
    disabledInputs.forEach(input => input.disabled = false);
    var formData = new FormData(this);

    disabledInputs.forEach(input => input.disabled = true);
    var nameData = {};
    var metaboliteFields = document.querySelectorAll('.substrates-name, .products-name');
    var allNamesEntered = true;
    metaboliteFields.forEach(function(input, index) {
        var key = input.name + (index + 1);
        var value = input.value;
        nameData[key] = value;

        if (input.value === '') {
            allNamesEntered = false;
        }
    });

    if (!allNamesEntered) {
        var errorMessage = 'Enter all metabolite names before creating reaction.';
        showErrorModal(errorMessage);
        window.scrollTo(0, 0);
        return; // Exit the function and do not submit form
    }
    submitBtn.disabled = true;
    loadingIndicator.style.display = 'flex';

    var skipAtomMapping = document.getElementById('skipAtomMapping').checked;
    formData.append('skipAtomMapping', skipAtomMapping);
    formData.append('nameData', JSON.stringify(nameData));
    formData.append('organs', JSON.stringify(organTags));
    const url = window.location.href;
    const parsedUrl = new URL(url);
    
    // Extract query parameters
    const params = new URLSearchParams(parsedUrl.search);
    const reactionId = params.get('reaction_id');
    // Get the value of reaction_id
    const action = params.get('action');
    formData.append('action', action);

    if (action === 'edit') {
        formData.append('action', action);

        formData.append('reaction_id', reactionId);
    }   
    formData.append('userID', sessionStorage.getItem('userID'));

    fetch(inputReactionUrl, {
        method: 'POST',
        body: formData,
        headers: {
            'X-Requested-With': 'XMLHttpRequest',
            'X-CSRFToken': csrfToken
        }
    })
    .then(response => response.json())
    .then(async (data) => {
        if (data.status === 'error') {
            showErrorModal(data.message);
            window.scrollTo(0, 0);
            submitBtn.disabled = false;
            loadingIndicator.style.display = 'none';
            return;
        }
        else if (action !== 'edit' && data.reaction_id) { // reactionId is taken from the response
            const userID = sessionStorage.getItem('userID');
            const reactionId = data.reaction_id; // Ensure reactionId is obtained from the response
            
            // Send additional request to save CreatedReaction
            fetch('create-reaction/', {
                method: 'POST',
                body: JSON.stringify({
                    user_id: userID,
                    reaction_id: reactionId
                }),
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                }
            })
            .then(response => response.json())
            .then(result => {
                if (result.success) {
                    console.log('CreatedReaction saved successfully.');
                } else {
                    console.error('Failed to save CreatedReaction:', result.error);
                }
            })
            .catch(error => {
                console.error('Error saving CreatedReaction:', error);
            });
        }
        
        setTimeout(function() {
            const redirectUrl = window.location.origin + "/?reaction_id=" + data.reaction_id;
            window.location.href = action === 'edit' ? redirectUrl + "&action=edit" : redirectUrl;
        }, 10); 
    })
    .catch(error => {
        console.error('Error:', error);
        const errorMessage = 'An unexpected error occurred.';
        showErrorModal(errorMessage);
        window.scrollTo(0, 0);
        submitBtn.disabled = false;
        loadingIndicator.style.display = 'none';
    });
});    


function persistSubsystemList(newSubsystem) {
    fetch('/update_subsystems/', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrfToken
        },
        body: JSON.stringify({ subsystems: [newSubsystem] })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            console.error('Error in updating subsystems:', data.message);
        } else {
            console.log('Subsystem list updated successfully.');
        }
    })
    .catch(error => console.error('Error:', error));
}


function showLoader() {
    document.getElementById('loadingIndicator').style.display = 'flex';
}

// Function to hide the loading indicator
function hideLoader() {
    document.getElementById('loadingIndicator').style.display = 'none';
}





function showErrorModal(message) {
    var errorMessageElement = document.getElementById('error-message');
    
    // Debugging output

    errorMessageElement.innerText = message;
    document.getElementById('error-message').style.display = 'block';
    document.getElementById('error-modal').style.display = 'block';
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