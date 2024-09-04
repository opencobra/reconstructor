// flags.js
flagsLoaded = false;
document.addEventListener('DOMContentLoaded', function () {
    const userId = sessionStorage.getItem('userID');
    // Utility function to fetch flags for a user
    async function fetchFlags(userId) {
        try {
            const response = await fetch(`flags/${userId}/`);
            const data = await response.json();
            if (data.status === 'success') {
                return data.flags;
            } else {
                console.error('Failed to load flags:', data.message);
                return [];
            }
        } catch (error) {
            console.error('Error loading flags:', error);
            return [];
        }
    }

    // Utility function to populate dropdown with flags
    function populateFlagDropdown(flags, dropdownMenu, selectedOption) {
        dropdownMenu.innerHTML = ''; // Clear existing flags

        flags.forEach(flag => {
            const option = document.createElement('div');
            option.style.padding = '10px';
            option.style.cursor = 'pointer';
            option.style.display = 'flex';
            option.style.alignItems = 'center';
            option.innerHTML = `<i class="fas fa-flag" style="color: ${flag.color}; margin-right: 10px;"></i><span>${flag.name_flag}</span>`;
            
            option.addEventListener('click', function() {
                selectedOption.innerHTML = `${flag.name_flag}<i class="fas fa-flag" style="color: ${flag.color}; margin-right: 10px; display: inline-block;"></i>`;
            });

            dropdownMenu.appendChild(option);
        });
    }

    // Load flags and populate dropdown
    async function loadFlags() {
        if (!flagsLoaded) {
            const flags = await fetchFlags(userId);
            const dropdownMenu = document.getElementById('dropdownMenu');
            const selectedOption = document.getElementById('selectedOption');
            populateFlagDropdown(flags, dropdownMenu, selectedOption);
            flagsLoaded = true;
        }
    }

    // Event listeners for dropdown interaction
    document.getElementById('customDropdown').addEventListener('click', function() {
        const dropdownMenu = document.getElementById('dropdownMenu');
        dropdownMenu.style.display = dropdownMenu.style.display === 'none' ? 'block' : 'none';
    });

    document.addEventListener('click', function(event) { 
        if (!event.target.closest('#customDropdown')) {
            document.getElementById('dropdownMenu').style.display = 'none';
        }
    });

    // Add new flag
    document.getElementById('submitCreateFlagCustom').addEventListener('click', async function () {
        const flagName = document.getElementById('newFlagNameCustom').value;
        const flagColor = document.getElementById('newFlagColorCustom').value;

        try {
            const response = await fetch('/add_flag/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                },
                body: JSON.stringify({ user_id: userId, name_flag: flagName, color: flagColor })
            });

            const data = await response.json();
            if (data.status === 'success') {
                flagsLoaded = false; // Reset flagsLoaded to reload flags
                loadFlags();
                document.getElementById('closeCreateFlagModalCustom').click();
            } else {
                alert(data.message);
            }
        } catch (error) {
            console.error('Error creating flag:', error);
        }
    });

    // Initial load if user is logged in
    if (userId !== null) {
        loadFlags();
    }
});
