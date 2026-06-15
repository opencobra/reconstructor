
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




// ---------------------------------------------------------------------------
// Delete reaction(s) — custom confirmation modal (intentionally hard to misclick)
//
// Both the per-row "Delete" button and the bulk "Delete" toolbar button funnel
// through the same modal. The confirm button stays disabled until the user types
// "delete"; pressing Enter in the field confirms when the text matches.
// ---------------------------------------------------------------------------
(function () {
    const modal = document.getElementById('deleteReactionModal');
    if (!modal) return;

    const messageEl = document.getElementById('deleteReactionMessage');
    const input = document.getElementById('deleteConfirmInput');
    const confirmBtn = document.getElementById('confirmDeleteReactionBtn');
    const cancelBtn = document.getElementById('cancelDeleteReactionBtn');

    let pendingIds = [];

    function openModal(ids, message) {
        pendingIds = ids;
        messageEl.textContent = message;
        input.value = '';
        confirmBtn.disabled = true;
        modal.style.display = 'flex';
        // Focus the field so the user can immediately type the confirmation word.
        setTimeout(() => input.focus(), 0);
    }

    function closeModal() {
        modal.style.display = 'none';
        pendingIds = [];
        input.value = '';
        confirmBtn.disabled = true;
    }

    function isConfirmed() {
        return input.value.trim().toLowerCase() === 'delete';
    }

    function performDelete() {
        if (!isConfirmed() || pendingIds.length === 0) return;
        confirmBtn.disabled = true;

        fetch(deleteReactionsUrl, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-CSRFToken': csrfToken,
            },
            body: JSON.stringify({
                userID: userID,
                reactionIds: pendingIds,
            }),
        })
            .then((response) => response.json())
            .then((data) => {
                if (data.status === 'success') {
                    location.reload();
                } else {
                    alert('Error deleting reaction(s): ' + (data.message || 'Unknown error.'));
                    closeModal();
                }
            })
            .catch((error) => {
                console.error('Error deleting reaction(s):', error);
                alert('An error occurred while deleting the reaction(s).');
                closeModal();
            });
    }

    // Enable the confirm button only once the user types the exact word.
    input.addEventListener('input', function () {
        confirmBtn.disabled = !isConfirmed();
    });

    // Enter confirms (only when the typed word matches).
    input.addEventListener('keydown', function (e) {
        if (e.key === 'Enter') {
            e.preventDefault();
            performDelete();
        }
    });

    confirmBtn.addEventListener('click', performDelete);
    cancelBtn.addEventListener('click', closeModal);

    // Click on the dimmed backdrop (outside the card) cancels.
    modal.addEventListener('click', function (e) {
        if (e.target === modal) closeModal();
    });

    // Per-row delete buttons (event delegation so it survives table rebuilds).
    document.addEventListener('click', function (e) {
        const btn = e.target.closest('.delete-btn');
        if (!btn) return;
        const reactionId = btn.getAttribute('data-reaction-id');
        if (!reactionId) return;
        const reactionName = btn.getAttribute('data-reaction-name') || 'this reaction';
        openModal([reactionId], `You are about to delete reaction "${reactionName}".`);
    });

    // Bulk "Delete" toolbar button — operates on the selected reactions.
    const deleteSelectedBtn = document.getElementById('deleteSelected');
    if (deleteSelectedBtn) {
        deleteSelectedBtn.addEventListener('click', function () {
            const ids = (window.checkedReactions || []).slice();
            if (ids.length === 0) {
                if (typeof showToast === 'function') {
                    showToast('Please select at least one reaction', '#f44336');
                } else {
                    alert('Please select at least one reaction.');
                }
                return;
            }
            const noun = ids.length === 1 ? 'reaction' : 'reactions';
            openModal(ids, `You are about to delete ${ids.length} selected ${noun}.`);
        });
    }
})();
