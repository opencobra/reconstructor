
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




document.getElementById('closeSaveReactionModal').addEventListener('click', function () {
    document.getElementById('saveReactionModal').style.display = 'none';
    document.getElementById('modalBackground').style.display = 'none';
});

function rgbToHex(rgb) {
    const result = rgb.match(/\d+/g);
    const r = parseInt(result[0]).toString(16).padStart(2, '0');
    const g = parseInt(result[1]).toString(16).padStart(2, '0');
    const b = parseInt(result[2]).toString(16).padStart(2, '0');
    return `#${r}${g}${b}`;
}


document.getElementById('submitSaveReaction').addEventListener('click', function () {
    var userID = sessionStorage.getItem('userID');
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    var shortNameInput = document.getElementById('reactionNameInput');
    var shortName = shortNameInput.value;
    var flag_name_element = document.getElementById('selectedOption');
    flag_name = flag_name_element.textContent;
    flag_name = flag_name.trim();   
    const flagIcon = flag_name_element.querySelector('i');
    var flag_color = flagIcon.style.color;
    flag_color = rgbToHex(flag_color);
    // var flag_color = document.getElementById('flagSelectCustom').options[document.getElementById('flagSelectCustom').selectedIndex].getAttribute('data-color');

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
        data.append('flag_name', flag_name);
        data.append('flag_color', flag_color);

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
                    var modal = document.getElementById('saveReactionModal');
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

});
