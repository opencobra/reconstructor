document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('addFlag').addEventListener('click', function() {
        document.getElementById('flagModal').style.display = 'block';
        fetchFlags(userID); 
    });

    document.querySelector('.custom-dropdown .dropdown-trigger').addEventListener('click', function() {
        document.querySelector('.custom-dropdown').classList.toggle('active');
    });

    document.querySelector('.dropdown-menu').addEventListener('click', function(event) {
        if (event.target.tagName === 'DIV') {
            const selectedOption = event.target;
            const color = selectedOption.getAttribute('data-color');
            const textContent = selectedOption.textContent;
    
            // Create the icon HTML with the color
            const iconHTML = `<i class="fas fa-flag" style="color: ${color}; margin-right: 10px;"></i>`;
    
            // Set the combined content into the dropdown trigger
            document.querySelector('.custom-dropdown .dropdown-trigger').innerHTML = textContent + iconHTML;
    
            document.querySelector('.custom-dropdown').classList.remove('active');
        }
    });

    function rgbToHex(rgb) {
        const result = rgb.match(/\d+/g);
        const r = parseInt(result[0]).toString(16).padStart(2, '0');
        const g = parseInt(result[1]).toString(16).padStart(2, '0');
        const b = parseInt(result[2]).toString(16).padStart(2, '0');
        return `#${r}${g}${b}`;
    }

    document.getElementById('AddFlagtosavedreaction').addEventListener('click', function() {
        const selectedReactionIds = getSelectedReactionIds();  // Implement this function to retrieve selected reactions
        const selectedFlagElement = document.querySelector('.dropdown-trigger');
        var selectedFlagName = selectedFlagElement.textContent.trim();  
        const flagIcon = selectedFlagElement.querySelector('i');
        var selectedFlagColor = rgbToHex(flagIcon.style.color);

        if (!selectedReactionIds.length) {
            console.error('No reactions selected');
            return;
        }

        if (selectedFlagName === 'None' || !selectedFlagColor) {
            console.error('No flag selected or flag color missing');
            return;
        }

        const data = {
            userID: userID,  
            reaction_ids: selectedReactionIds,
            flag_name: selectedFlagName,
            flag_color: selectedFlagColor
        };

        fetch('/saved_reactions/save_flags_in_saved_reactions/', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/x-www-form-urlencoded',
                'X-CSRFToken': csrfToken  
            },
            body: JSON.stringify(data) 
        })
        .then(response => response.json())
        .then(data => {
            if (data.status === 'success') {
                console.log('Flag added to selected reactions successfully');
                location.reload();
            } else {
                console.error(data.message);
            }
        })
        .catch(error => console.error('Error:', error));
    });

    function getSelectedReactionIds() {
        let selected = [];
        selected = checkedReactions;
        return selected;
    }

    document.getElementById('closeflagmodal').addEventListener('click', function() {
        document.getElementById('flagModal').style.display = 'none';
    });

    document.getElementById('createFlagButtonCustom').addEventListener('click', function() {
        const newFlagFields = document.getElementById('newFlagFields');
        newFlagFields.style.display = (newFlagFields.style.display === 'block') ? 'none' : 'block';
    });

    window.addEventListener('click', function(event) {
        if (event.target == document.getElementById('flagModal')) {
            document.getElementById('flagModal').style.display = 'none';
        }
    });

    document.getElementById('saveFlagButton').addEventListener('click', function() {
        const flagNameInput = document.getElementById('flagNameInput').value;
        const flagColorInput = document.getElementById('flagColorInput').value;

        if (flagNameInput && flagColorInput) {
            const data = {
                user_id: userID,  
                name_flag: flagNameInput,
                color: flagColorInput
            };

            fetch('/saved_reactions/add_flag/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken  
                },
                body: JSON.stringify(data)
            })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    resetDropdown();
                    document.getElementById('newFlagFields').style.display = 'none';
                } else {
                    console.error(data.message);
                }
            })
            .catch(error => console.error('Error:', error));
        } else {
            alert('Flag name and color are required');
        }
    });

    function fetchFlags(userId) {
        fetch(`/saved_reactions/flags/${userId}/`)
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    populateDropdown(data.flags);
                } else {
                    console.error(data.message);
                }
            })
            .catch(error => console.error('Error:', error));
    }

    function populateDropdown(flags) {
        const dropdownMenu = document.querySelector('.custom-dropdown .dropdown-menu');
        dropdownMenu.innerHTML = '';  // Clear the dropdown before populating it
    
        flags.forEach(flag => {
            const div = document.createElement('div');
            div.innerHTML = `
                <i class="fas fa-flag" style="color: ${flag.color}; margin-right: 10px;"></i>
                <span>${flag.name_flag}</span>`;
            div.setAttribute('data-color', flag.color);
            dropdownMenu.appendChild(div);
        });
    }

    function resetDropdown() {
        fetchFlags(userID);  // Re-fetch flags and repopulate dropdown
    }

});
