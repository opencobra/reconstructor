document.addEventListener('DOMContentLoaded', function () {
    const userId = sessionStorage.getItem('userID');
    const createFlagModal = document.getElementById('createFlagModalCustom');
    const modalOverlay = document.getElementById('modalOverlayCustom');
    const createFlagButton = document.getElementById('createFlagButtonCustom');
    const closeCreateFlagModal = document.getElementById('closeCreateFlagModalCustom');
    const saveReactionModal = document.getElementById('saveReactionModal');
    const submitCreateFlagButton = document.getElementById('submitCreateFlagCustom');
    const flagColorDisplay = document.getElementById('flagColorDisplayCustom'); // Reference to color display box

    // Function to open the create flag modal and close the save reaction modal
    function openCreateFlagModal() {
        createFlagModal.classList.add('active');
        modalOverlay.classList.add('active');
        saveReactionModal.style.display = 'none';
        document.getElementById('modalBackground').style.display = 'none';
    }

    // Function to close the create flag modal and open the save reaction modal
    function closeCreateFlagModalAndReopenSaveReaction() {
        createFlagModal.classList.remove('active');
        modalOverlay.classList.remove('active');
        document.getElementById('modalBackground').style.display = 'block';
        saveReactionModal.style.display = 'block';
    }

    // Open the create flag modal when the create flag button is clicked
    createFlagButton.addEventListener('click', openCreateFlagModal);

    // Close the create flag modal when the close button or overlay is clicked
    closeCreateFlagModal.addEventListener('click', closeCreateFlagModalAndReopenSaveReaction);
    modalOverlay.addEventListener('click', closeCreateFlagModalAndReopenSaveReaction);

    // Function to update the color display box based on selected flag
    function updateColorDisplay() {
        const flagSelect = document.getElementById('flagSelectCustom');
        const selectedOption = flagSelect.options[flagSelect.selectedIndex];
        const color = selectedOption.getAttribute('data-color');

        if (color) {
            flagColorDisplay.style.backgroundColor = color;
        } else {
            flagColorDisplay.style.backgroundColor = 'transparent'; // Reset if no color
        }
    }

    // Load existing flags for the user
    function loadFlags() {
        fetch(`/flags/${userId}/`)
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    const flagSelect = document.getElementById('flagSelectCustom');
                    flagSelect.innerHTML = '<option value="">None</option>';
                    data.flags.forEach(flag => {
                        const option = document.createElement('option');
                        option.value = flag.id;
                        option.textContent = flag.name_flag;

                        // Attach flag color as a data attribute
                        option.setAttribute('data-color', flag.color);

                        flagSelect.appendChild(option);
                    });
                    updateColorDisplay(); // Update color display after loading flags
                } else {
                    console.error('Failed to load flags:', data.message);
                }
            })
            .catch(error => console.error('Error loading flags:', error));
    }

    // Submit the new flag and reload flags
    submitCreateFlagButton.addEventListener('click', function () {
        const flagName = document.getElementById('newFlagNameCustom').value;
        const flagColor = document.getElementById('newFlagColorCustom').value;

        fetch('/add_flag/', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ user_id: userId, name_flag: flagName, color: flagColor })
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    loadFlags(); // Reload the flags after adding a new one
                    closeCreateFlagModalAndReopenSaveReaction(); // Close flag modal and reopen save reaction modal
                } else {
                    alert(data.message);
                }
            })
            .catch(error => console.error('Error creating flag:', error));
    });

    // Event listener to update color display when a flag is selected
    document.getElementById('flagSelectCustom').addEventListener('change', updateColorDisplay);

    var user = sessionStorage.getItem('userID');

    if (user !== null) {
        loadFlags(); // Load flags if user is logged in
    }
});
