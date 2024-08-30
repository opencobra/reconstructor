document.addEventListener('DOMContentLoaded', function() {
    const reactionSelect = document.getElementById('reactionSelect');
    const reactionInput = document.getElementById('reactionAbbreviation');
    const submitButton = document.getElementById('submitButton');
    const cancelButton = document.getElementById('cancel-getrxnfromvmhrhea-button');
    $('#getrxnfromvmhrhea').click(function() {
        $('.ui.modal.getrxnfromvmhrhea').modal('show');
    });

    reactionSelect.addEventListener('change', function() {
        handleReactionSelectChange();
    });

    submitButton.addEventListener('click', function() {
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
                    alert('Error: ' + data.error);
                } else {
                    $('.ui.modal.getrxnfromvmhrhea').modal('hide');

                    updateFormFields(data);
                }
            })
            .catch(error => console.error('Error:', error.message));
        }
    }

    function handleRheaFetch() {
        if (reactionInput.value.trim()) {
            fetch('reaction_view', {
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
                    alert('Error: ' + data.error);
                } else {
+
                    $('.ui.modal.getrxnfromvmhrhea').modal('hide');

                    updateFormFields(data);
                }
            })
            .catch(error => console.error('Error:', error));
        }
    }

    // Initial call to handleReactionSelectChange to set initial state
    handleReactionSelectChange();
});
