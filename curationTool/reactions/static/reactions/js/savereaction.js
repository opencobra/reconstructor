// saved_reactions.js

document.addEventListener('DOMContentLoaded', function () {
    const saveReactionButton = document.getElementById('saveReactionButton');

    async function checkAndSaveReaction() {
        const urlParams = new URLSearchParams(window.location.search);
        const reactionId = urlParams.get('reaction_id');
        const userId = sessionStorage.getItem('userID');

        if (!reactionId) {
            alert('Reaction not created.');
            return;
        }

        const url = `/check-reaction`;

        try {
            let response = await fetch(url, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken
                },
                body: JSON.stringify({ user_id: userId, reaction_id: reactionId })
            });

            const data = await response.json();

            if (data.is_reaction_saved) {
                alert('Reaction is already saved.');
            } else {
                var modal = document.getElementById('saveReactionModal');
                modal.style.display = 'block';
                document.getElementById('modalBackground').style.display = 'block';

                if (!flagsLoaded) {
                    loadFlags();
                }
            }
        } catch (error) {
            console.error('Error checking reaction:', error);
        }
    }
    document.getElementById('submitSaveReaction').addEventListener('click', function () {
        const userID = sessionStorage.getItem('userID');
        const urlParams = new URLSearchParams(window.location.search);
        const reactionId = urlParams.get('reaction_id');
        const shortNameInput = document.getElementById('reactionNameInput');
        const shortName = shortNameInput.value;
        const flagNameElement = document.getElementById('selectedOption');
        let flagName = flagNameElement.textContent.trim();
        const flagIcon = flagNameElement.querySelector('i');
        let flagColor = flagIcon ? flagIcon.style.color : '';
    
        // Ensure the rgbToHex function is defined
        function rgbToHex(rgb) {
            const result = rgb.match(/\d+/g);
            const r = parseInt(result[0]).toString(16).padStart(2, '0');
            const g = parseInt(result[1]).toString(16).padStart(2, '0');
            const b = parseInt(result[2]).toString(16).padStart(2, '0');
            return `#${r}${g}${b}`;
        }
    
        // Convert the color if available
        flagColor = flagColor ? rgbToHex(flagColor) : '';
    
        shortNameInput.setCustomValidity(''); // Clear any previous custom validity message
    
        if (userID && reactionId) {
            if (!shortName) {
                alert('Please enter a short name for the reaction.');
                shortNameInput.setCustomValidity('Please enter a short name for the reaction.');
                shortNameInput.reportValidity();
                return; // Prevent form submission
            }
    
            const data = new FormData();
            data.append('userID', userID);
            data.append('reaction_id', reactionId);
            data.append('short_name', shortName);
            data.append('flag_name', flagName);
            data.append('flag_color', flagColor);
    
            fetch(saveReaction, { // Ensure the correct URL endpoint here
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
                    document.getElementById('saveReactionModal').style.display = 'none';
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
    
    saveReactionButton.addEventListener('click', checkAndSaveReaction);
});
