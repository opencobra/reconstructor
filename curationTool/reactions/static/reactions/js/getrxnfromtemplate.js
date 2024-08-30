document.addEventListener('DOMContentLoaded', function() {
    // Initialize dropdown
    var reactionSelect = document.getElementById('reactionSelectemp');
    
    // Show modal on button click
    document.getElementById('getrxnfromtemplate').addEventListener('click', function() {
            $('.ui.modal.getrxntemplate').modal('show');

    });

    // Enable apply button when a valid option is selected
    reactionSelect.addEventListener('change', function() {
        var selectedValue = reactionSelect.value;
        var applyTemplateButton = document.getElementById('applyTemplate');
        if (selectedValue) {
            applyTemplateButton.disabled = false;
            applyTemplateButton.style = 'block';
        } else {
            applyTemplateButton.disabled = true;
            applyTemplateButton.style = 'none';
        }
    });

    document.getElementById('applyTemplate').addEventListener('click', async function() {
        var selectedValue = reactionSelect.value;
        await handleTemplateButtonClick(selectedValue);
        document.querySelector('.ui.modal.getrxntemplate').classList.remove('show');
    });

    async function handleTemplateButtonClick(option) {
        try {
            const response = await fetch('/get_rxn_template/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ reaction_type: option })
            });

            if (response.ok) {
                const data = await response.json();
                $('.ui.modal.getrxntemplate').modal('hide');

                await updateFormFields(data);
            } else {
                console.error('Failed to fetch reaction template');
            }
        } catch (error) {
            console.error('Error:', error);
        }
    }
});
