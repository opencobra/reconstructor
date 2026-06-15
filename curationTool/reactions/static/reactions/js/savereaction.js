/**
 * Convert a CSS color into a #rrggbb hex string.
 *
 * Flags are stored as hex (the color picker emits e.g. "#ff0000"), but reading
 * a flag icon's `style.color` back from the DOM yields "rgb(255, 0, 0)". This
 * normalizes it back to hex so the saved flag_color matches the stored flag.
 * Already-hex values are returned unchanged.
 */
function rgbToHex(color) {
    if (!color) {
        return '';
    }
    color = color.trim();
    if (color.charAt(0) === '#') {
        return color;
    }
    const match = color.match(/^rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i);
    if (!match) {
        return color; // Unrecognized format; pass through unchanged.
    }
    const toHex = (n) => parseInt(n, 10).toString(16).padStart(2, '0');
    return `#${toHex(match[1])}${toHex(match[2])}${toHex(match[3])}`;
}

document.getElementById("confidenceScoreInput").addEventListener("input", function () {
    let value = this.value.trim();

    if (value === "") {
        this.value = ""; // Allow blank (null)
    } else {
        let numValue = parseInt(value, 10);
        if (isNaN(numValue) || numValue < 1 || numValue > 4) {
            this.value = ""; // Reset invalid inputs
        }
    }
});
async function checkAndSaveReaction() {
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    const userId = sessionStorage.getItem('userID');

    if (!reactionId) {
        Notify.error('Reaction not created.');
        return;
    }
    try {
        let response = await fetch(alreadySavedURL, {
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
            Notify.info('Reaction is already saved.');
        } else {
            var modal = document.getElementById('saveReactionModal');
            modal.style.display = 'block';
            document.getElementById('reactionDescriptionInput').value = sessionStorage.getItem('lastTemplateDescription') || '';
            document.getElementById('modalBackground').style.display = 'block';

            if (!flagsLoaded) {
                loadFlags();
            }
        }
    } catch (error) {
        console.error('Error checking reaction:', error);
    }
}
async function reactionNameExists(shortName, userId) {
    try {
        const response = await fetch(reactionNameExistsURL, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ short_name: shortName, user_id: userId })
        });

        const data = await response.json();
        return data.is_name_saved;
    } catch (error) {
        console.error('Error checking reaction name:', error);
        return false;
    }
}
document.addEventListener('DOMContentLoaded', function () {
    const saveReactionButton = document.getElementById('saveReactionButton');

    document.getElementById('submitSaveReaction').addEventListener('click', async function () {
        const userID = sessionStorage.getItem('userID');
        const urlParams = new URLSearchParams(window.location.search);
        const reactionId = urlParams.get('reaction_id');
        const shortNameInput = document.getElementById('reactionNameInput');
        const shortName = shortNameInput.value;
        const descriptionInput = document.getElementById('reactionDescriptionInput'); // New description input
        const reactionDescription = descriptionInput ? descriptionInput.value : ''; // Get description
        const flagNameElement = document.getElementById('selectedOption');
        let flagName = flagNameElement.textContent.trim();
        const flagIcon = flagNameElement.querySelector('i');
        let flagColor = flagIcon ? flagIcon.style.color : '';
    
        // Convert the color if available
        flagColor = flagColor ? rgbToHex(flagColor) : '';

        // Get Confidence Score
        const confidenceScoreInput = document.getElementById("confidenceScoreInput");
        let confidenceScore = confidenceScoreInput.value.trim() === "" ? null : parseInt(confidenceScoreInput.value);

        // Ensure only valid numbers are sent (1-4) OR null
        if (confidenceScore !== null && (confidenceScore < 1 || confidenceScore > 4)) {
            Notify.error("Confidence Score must be between 1 and 4.");
            return;
        }

        shortNameInput.setCustomValidity(''); // Clear any previous custom validity message
    
        if (userID && reactionId) {
            if (!shortName) {
                Notify.error('Please enter a short name for the reaction.');
                shortNameInput.setCustomValidity('Please enter a short name for the reaction.');
                shortNameInput.reportValidity();
                return; // Prevent form submission
            }
    
            const data = new FormData();
            data.append('userID', userID);
            data.append('reaction_id', reactionId);
            data.append('short_name', shortName);
            data.append('description', reactionDescription); // Include description
            data.append('flag_name', flagName);
            data.append('flag_color', flagColor);
            if (confidenceScore !== null) {
                data.append("confidence_score", confidenceScore); // Include confidence score
            }
            const nameExists = await reactionNameExists(shortName, userID);
            if (nameExists) {
                const userConfirmation = await Notify.confirm({
                    title: 'Name already exists',
                    message: 'A reaction with this name already exists. Continue anyway?\nYou will have two reactions with the same name.',
                    confirmText: 'Continue',
                    cancelText: 'Cancel',
                });
                if (!userConfirmation) {
                    return; // Prevent form submission and keep the modal open
                }
            }            
            fetch(saveReaction, {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken
                },
                body: data,
            })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    Notify.success("Reaction saved successfully!");
                    document.getElementById('saveReactionModal').style.display = 'none';
                    document.getElementById('modalBackground').style.display = 'none';
                } else {
                    Notify.error("Error: " + (data.message || "Failed to save the reaction."));
                }
            })
            .catch(error => console.error('Error:', error));
        } else if (!reactionId) {
            Notify.error("Create the reaction first.");
        } else {
            Notify.error("Please log in.");
        }
    });

    saveReactionButton.addEventListener('click', checkAndSaveReaction);
});
