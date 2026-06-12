document.addEventListener('DOMContentLoaded', function() {
    const reactionSelect = document.getElementById('reactionSelect');
    const reactionInput = document.getElementById('reactionAbbreviation');
    const submitButton = document.getElementById('submitButton');

    $('#getrxnfromvmhrhea').click(function() {
        $('.ui.modal.getrxnfromvmhrhea').modal('show');
    });

    reactionSelect.addEventListener('change', function() {
        handleReactionSelectChange();
    });

    submitButton.addEventListener('click', function() {
        // Disable the button immediately after click
        submitButton.disabled = true;

        // Handle the actions based on the selection
        if (reactionSelect.value === 'VMH') {
            handleVmhFetch();
        } else if (reactionSelect.value === 'RHEA') {
            handleRheaFetch();
        }
    });

    function handleReactionSelectChange() {
        const selectedValue = reactionSelect.value;

        if (selectedValue === 'VMH') {
            reactionInput.disabled = false;
            reactionInput.placeholder = 'Enter VMH Abbreviation';
            submitButton.style.display = 'block';
        } else if (selectedValue === 'RHEA') {
            reactionInput.disabled = false;
            reactionInput.placeholder = 'Enter RHEA ID';
            submitButton.style.display = 'block';
        } else {
            reactionInput.disabled = true;
            reactionInput.placeholder = '';
            submitButton.style.display = 'none';
        }
    }

    function handleVmhFetch() {
        if (reactionInput.value.trim()) {
            fetch('/get_from_vmh/', {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ reactionAbbreviation: reactionInput.value.trim() })
            })
            .then(response => response.json())
            .then(data => {
                if (data.error) {
                    Notify.error('Error: ' + data.error);
                } else {
                    $('.ui.modal.getrxnfromvmhrhea').modal('hide');
                    updateFormFields(data);
                }
            })
            .catch(error => console.error('Error:', error.message))
            .finally(() => {
                // Re-enable the button after the fetch completes
                submitButton.disabled = false;
            });
        } else {
            // Re-enable the button if no input is provided
            submitButton.disabled = false;
        }
    }

    function handleRheaFetch() {
        if (reactionInput.value.trim()) {
            fetch('fetch_rhea_rxn', {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ reactionAbbreviation: reactionInput.value.trim() })
            })
            .then(response => response.json())
            .then(data => {
                if (data.error) {
                    Notify.error('Error: ' + data.error);
                } else {
                    $('.ui.modal.getrxnfromvmhrhea').modal('hide');
                    updateFormFields(data);
                }
            })
            .catch(error => console.error('Error:', error))
            .finally(() => {
                // Re-enable the button after the fetch completes
                submitButton.disabled = false;
            });
        } else {
            // Re-enable the button if no input is provided
            submitButton.disabled = false;
        }
    }

    // Initial call to handleReactionSelectChange to set initial state
    handleReactionSelectChange();
});
