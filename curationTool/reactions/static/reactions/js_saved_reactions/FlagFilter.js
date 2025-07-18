document.addEventListener('DOMContentLoaded', function() {
    const filterButton = document.getElementById('FilterFlag');
    const dropdown = document.getElementById('Filteredflag');
    let selectedFlags = []; // Array to hold selected flags

    filterButton.addEventListener('click', function() {
        const userId = userID; // Replace with the actual user ID

        if (dropdown.style.display === 'inline-block') {
            // If dropdown is already open, close it
            toggleDropdown();
            return;
        }

        dropdown.style.display = 'inline-block';

        // Fetch flags and populate dropdown
        fetchFlags(userId);
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
        // Clear existing items in the dropdown container
        dropdown.innerHTML = '';
    
        // Create the "None" option with an arrow
        const noneOption = document.createElement('div');
        noneOption.className = 'dropdown-item';
        noneOption.innerHTML = 'None';
        noneOption.setAttribute('data-value', '');
        noneOption.addEventListener('click', () => {
            selectedFlags = []; // Clear selected flags
            filterTable(); // Show all rows
            resetFlagColors(); // Reset all flag colors
        });
        const arrow = document.createElement('button');
        arrow.className = 'fas fa-chevron-up';
        arrow.id = 'arrow';
    
        arrow.addEventListener('click', toggleDropdown);
        dropdown.appendChild(noneOption);
        dropdown.appendChild(arrow);
    
        // Create dropdown items for each flag
        flags.forEach(flag => {
            const item = document.createElement('div');
            item.className = 'dropdown-item';
    
            const icon = document.createElement('i');
            icon.className = 'fas fa-flag';
            icon.style.color = flag.color;
            icon.style.marginRight = '10px';
    
            item.appendChild(icon);
            item.appendChild(document.createTextNode(flag.name_flag));
            item.setAttribute('data-color', flag.color.toLowerCase());
    
            item.addEventListener('click', () => {
                const color = flag.color.toLowerCase();
                
                if (selectedFlags.includes(color)) {
                    // If flag is already selected, remove it and restore original color
                    selectedFlags = selectedFlags.filter(c => c !== color);
                    item.classList.remove('selected'); // Remove highlight class
                } else {
                    // Otherwise, add it to the selected flags and apply the highlight
                    selectedFlags.push(color);
                    item.classList.add('selected'); // Add highlight class
                }
    
                filterTable(); // Apply filter based on selected flags
            });
    
            dropdown.appendChild(item);
        });

        // Reapply selected flags' highlights after populating the dropdown
        reapplySelectedFlagHighlights();
    }

    function reapplySelectedFlagHighlights() {
        const items = dropdown.querySelectorAll('.dropdown-item');
        items.forEach(item => {
            const color = item.getAttribute('data-color');
            if (selectedFlags.includes(color)) {
                item.classList.add('selected'); // Apply consistent highlight
            } else {
                item.classList.remove('selected'); // Remove highlight
            }
        });
    }
    function hexToRgb(hex) {
        hex = hex.replace(/^#/, '');
        
        let bigint = parseInt(hex, 16);
        let r = (bigint >> 16) & 255;
        let g = (bigint >> 8) & 255;
        let b = bigint & 255;
        
        return `rgb(${r}, ${g}, ${b})`;
    }

    function filterTable() {
        const reactionListBody = document.querySelector('#reactionList tbody');
        const rows = reactionListBody.querySelectorAll('tr'); // Only select rows in tbody
        
        for (const row of rows) {
            // ignore first row (header)
            if (row.rowIndex === 0) {
                continue;
            }
            const flagIcons = row.querySelectorAll('.fas.fa-flag');
            let rowMatches = false;

            for (const icon of flagIcons) {
                const iconColor = icon.style.color;
                if (selectedFlags.some(color => hexToRgb(color) === iconColor)) {
                    rowMatches = true;
                    break;
                }
            }

            row.style.display = (selectedFlags.length === 0 || rowMatches) ? '' : 'none';
        }
    }

    function resetFlagColors() {
        const items = dropdown.querySelectorAll('.dropdown-item');
        items.forEach(item => {
            item.classList.remove('selected'); // Reset all to non-highlighted
        });
    }

    function toggleDropdown() {
        if (dropdown.style.display === 'none') {
            // Show the dropdown
            dropdown.style.display = 'inline-block';
            
            // Reapply the selected flags' highlights
            reapplySelectedFlagHighlights();
        } else {
            // Hide the dropdown
            dropdown.style.display = 'none';
        }
        
        // Toggle the filter button visibility
        filterButton.style.display = (dropdown.style.display === 'none') ? 'inline-block' : 'none';
    }
});