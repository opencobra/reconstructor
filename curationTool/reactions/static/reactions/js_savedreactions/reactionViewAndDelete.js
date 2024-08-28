
// Attach event listener to buttons with the class 'view-btn'
document.querySelectorAll('.view-btn').forEach(function(button, index) {
    button.addEventListener('click', function() {
        const reactionData = reactions[index];
        const reactionId = reactionData.pk;
        window.location.href = `/?reaction_id=${reactionId}&action=edit`;
    });
});

function cloneReaction(reactionId) {
    var data = new FormData();
    data.append('reaction_id', reactionId);
    data.append('userID', sessionStorage.getItem('userID'));
    
    var name = prompt('Enter a name for the cloned reaction:')
    if (!name) {
        return;
    }
    data.append('name', name);

    fetch("reactions/clone/", {
        method: 'POST',
        headers: {
            'X-Requested-With': 'XMLHttpRequest',
            'X-CSRFToken': csrfToken,
        },
        body: data  // Use FormData object here
    })
    .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok');
        }
        return response.json();
    })
    .then(data => {
        if (data.status === 'success') {
            alert('Reaction cloned successfully!');
            location.reload();

        } else {
            alert(`Error cloning reaction: ${data.error}`);
        }
    })
    .catch(error => {
        console.error('There was a problem with the fetch operation:', error);
        alert('An error occurred while cloning the reaction.');
    });
}


// get user flag 
// add user flag

function savedReactionModal(userID, reactionId, saveReaction) {
    var shortNameInput = document.getElementById('reactionNameInputSavedPage');
    var shortName = shortNameInput.value;
    var flagSelect = document.getElementById('flagSelectSavedPage');
    var flag_name = flagSelect.options[flagSelect.selectedIndex].text;
    var flag_color = flagSelect.options[flagSelect.selectedIndex].getAttribute('data-color');

    // Clear any previous custom validity message
    shortNameInput.setCustomValidity('');

    if (userID && reactionId) {
        // if (!shortName) {
        //     // Alert if the short name is not provided
        //     alert('Please enter a short name for the reaction.');
        //     shortNameInput.setCustomValidity('Please enter a short name for the reaction.');
        //     shortNameInput.reportValidity();
        //     return; // Prevent form submission
        // }
        
        var data = new FormData();
        data.append('userID', userID);
        data.append('reaction_id', reactionId);
        data.append('short_name', shortName);
        data.append('flag_name', flag_name);
        data.append('flag_color', flag_color);
        console.log('Saving reaction:', reactionId);
        console.log('Short name:', shortName);  
        console.log('Flag name:', flag_name);
        console.log('Flag color:', flag_color);

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
                alert("Reaction saved successfully!");
                var modal = document.getElementById('savedReactionsModal');
                modal.style.display = 'none';
                document.getElementById('modalBackground').style.display = 'none';
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
}




    // ask user for new reaction name
    // make fetch call to clone reaction (makes identical copy of reaction with new name and adds to users saved reactions)
    // refresh page
function confirmDelete(reactionName) {
    return confirm('Are you sure you want to delete reaction ' + reactionName + '?');
}
