document.addEventListener('DOMContentLoaded', function () {
    const userId = sessionStorage.getItem('userID');
    const createFlagModal = document.getElementById('createFlagModalCustom');
    const modalOverlay = document.getElementById('modalOverlayCustom');
    const createFlagButton = document.getElementById('createFlagButtonCustom');
    const closeCreateFlagModal = document.getElementById('closeCreateFlagModalCustom');
    const saveReactionModal = document.getElementById('saveReactionModal');
    const submitCreateFlagButton = document.getElementById('submitCreateFlagCustom');

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

    // Load existing flags for the user
    function loadFlags() {
        fetch(`/flags/${userId}/`)
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    const dropdownMenu = document.getElementById('dropdownMenu');
                    const selectedOption = document.getElementById('selectedOption');
                
                    data.flags.forEach(flag => {
                        const option = document.createElement('div');
                        option.style.padding = '10px';
                        option.style.cursor = 'pointer';
                        option.style.display = 'flex';
                        option.style.alignItems = 'center';
                        
                        option.innerHTML = `
                        <span>${flag.name_flag}</span>
                        <i class="fas fa-flag" style="color: ${flag.color}; margin-right: 10px;"></i>
                        
                    `;
                    
                    option.addEventListener('click', function() {
                        selectedOption.innerHTML = `
                        ${flag.name_flag}    
                        <i class="fas fa-flag" style="color: ${flag.color}; margin-right: 10px; display: inline-block;"></i>
                            
                        `;
                        dropdownMenu.style.display = 'none';
                    });
                    
                
                        dropdownMenu.appendChild(option);
                    });
                
                    selectedOption.addEventListener('click', function() {
                        dropdownMenu.style.display = dropdownMenu.style.display === 'none' ? 'block' : 'none';
                    });
                
                    document.addEventListener('click', function(event) {
                        if (!event.target.closest('#customDropdown')) {
                            dropdownMenu.style.display = 'none';
                        }
                    });
                
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

    var user = sessionStorage.getItem('userID');

    if (user !== null) {
        loadFlags(); // Load flags if user is logged in
    }
});
