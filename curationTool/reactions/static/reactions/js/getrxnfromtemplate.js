document.addEventListener('DOMContentLoaded', function () {
    const modalButton = document.getElementById('getrxnfromtemplateBtn');
    const reactionField = document.getElementById('reactionField'); 
    const reactionDropdown = document.getElementById('reactionDropdown');
    let templateList = [];
    let templatesFetched = { value: false }; // Mutable tracking

    // Populate dropdown with available templates
    function populateDropdown(list) {
        reactionDropdown.innerHTML = ''; // Clear existing dropdown content
        list.forEach(template => {
            const element = document.createElement('div');
            element.textContent = template;
            element.addEventListener('click', function () {
                reactionField.value = this.textContent; // Set selected value in the input field
                reactionDropdown.style.display = 'none'; // Hide dropdown
            });
            reactionDropdown.appendChild(element);
        });
        reactionDropdown.style.display = list.length > 0 ? 'block' : 'none';
    }

    // Show dropdown on input click
    reactionField.addEventListener('click', function () {
        populateDropdown(templateList); // Show the full list on click
    });

    // Filter dropdown based on user input
    reactionField.addEventListener('keyup', function () {
        const inputVal = this.value.toLowerCase();
        const matches = templateList.filter(template => template.toLowerCase().includes(inputVal));
        populateDropdown(matches);
    });

    // Hide dropdown when clicking outside
    document.addEventListener('click', function (event) {
        if (!reactionDropdown.contains(event.target) && !reactionField.contains(event.target)) {
            reactionDropdown.style.display = 'none';
        }
    });

    modalButton.addEventListener('click', async function () {
        $('#getrxnfromtemplateModal').modal('show');
        await window.ReactionUtils.fetchTemplates(templateList, templatesFetched);
    });

    // Apply selected template
    document.getElementById('applyTemplate').addEventListener('click', async function () {
        const selectedValue = reactionField.value;
        if (!templateList.includes(selectedValue)) {
            Notify.error('Please select a valid template from the dropdown.');
            return;
        }

        try {
            const userId = sessionStorage.getItem('userID');
            const response = await fetch('/get_rxn_template/', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ reaction_type: selectedValue, userID: userId })
            });
            if (response.ok) {
                const data = await response.json();
                sessionStorage.setItem('lastTemplateDescription', data.description || '');
                await updateFormFields(data);
                $('.ui.modal.getrxntemplate').modal('hide');
            } else {
                console.error('Failed to fetch template');
            }
        } catch (error) {
            console.error('Error:', error);
        }
    });
});
