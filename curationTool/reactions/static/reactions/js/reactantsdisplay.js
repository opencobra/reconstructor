document.addEventListener('DOMContentLoaded', async function () {
    $('#savedMetabolitesModal').modal();

    // Check if we should open the saved metabolites modal (from redirect)
    if (sessionStorage.getItem('openSavedMetabolites') === 'true') {
        sessionStorage.removeItem('openSavedMetabolites');
        // Small delay to ensure everything is loaded
        setTimeout(() => {
            loadSavedMetabolites();
            $('#savedMetabolitesModal').modal('show');
        }, 500);
    }

    // Custom dropdown functionality
    const userDropdown = document.getElementById('userDropdown');
    const userDropdownTrigger = document.getElementById('userDisplay');
    
    if (userDropdownTrigger && userDropdown) {
        userDropdownTrigger.addEventListener('click', function(e) {
            e.stopPropagation();
            userDropdown.classList.toggle('active');
        });
        
        // Close dropdown when clicking outside
        document.addEventListener('click', function(e) {
            if (!userDropdown.contains(e.target)) {
                userDropdown.classList.remove('active');
            }
        });
        
        // Handle dropdown item clicks
        document.getElementById('dropdown-saved-reactions').addEventListener('click', function() {
            const userId = sessionStorage.getItem('userID');
            if (userId) {
                window.location.href = '/saved_reactions';
            } else {
                var errorMessage = 'Please login to view saved reactions.';
                showErrorModal(errorMessage);
            }
            userDropdown.classList.remove('active');
        });
        
        document.getElementById('dropdown-saved-metabolites').addEventListener('click', function() {
            loadSavedMetabolites();
            $('#savedMetabolitesModal').modal('show');
            userDropdown.classList.remove('active');
        });
    }

    // Add event listeners for nav links
    const aboutNav = document.getElementById('about-nav');
    if (aboutNav) {
        aboutNav.addEventListener('click', function(e) {
            e.preventDefault();
            window.location.href = aboutUrl;
        });
    }

    if (sessionStorage.getItem('userID') !== null) {
        username = sessionStorage.getItem('userName');
        userID = sessionStorage.getItem('userID');
        const userDisplayEl = document.getElementById('userDisplay');
        if (userDisplayEl) {
            userDisplayEl.innerHTML = `<i class="fas fa-user"></i><span>User: ${username}</span><i class="fas fa-chevron-down dropdown-arrow"></i>`;
        }
        document.getElementById('loginButton').textContent = 'Log out';
        setLoggedInStatusBasedOnUrl('');
        fetch(setSessionUser, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            },
            body: JSON.stringify({ 'userID': userID })
        })
            .then(response => response.json())
            .then(data => {
                if (data.status === 'success') {
                    console.log('Session user set successfully:', data.message);
                } else {
                    var errorMessageContainer = 'Error in setting session user: ' + data.message;
                    showErrorModal(errorMessageContainer);   
                }            })

    }

    setupTooltips();
    createnewreaction();
    createGeneInfoInput();
    setupdate();
    displayreactioninfo(reactionData = null);
    attachEventListenersToSelects();
    toggleStructure();
    updateAtomChargeCounters();
    initResponsiveHeaders(); // Initialize responsive column headers
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    const action = urlParams.get('action');
    if (reactionId || action === 'edit') {
        fetch(getReaction + reactionId)
            .then(response => {
                if (!response.ok) {
                    console.error(`HTTP error! status: ${response.status}`);
                    throw new Error(`HTTP error! status: ${response.status}`);
                }
                return response.json();
            })
            .then(async reactionData => {
                await updateFormFields(reactionData);
                confirmAll();
                displayDivs(reactionData);
                if (reactionData.short_name) {
                    let reactionDataName = reactionData.short_name;
                    if(reactionData.short_name.length>20){
                        reactionDataName = reactionData.short_name.substring(0, 30) + "...";
                    }
                    let reactionDataDescription = reactionData.description;
                    reactionStatusInfo = {
                        'name': reactionDataName,
                        'description': reactionDataDescription
                      };
                    setLoggedInStatusBasedOnUrl(reactionStatusInfo);
                    DisplayTag(reactionData.Organs);
                }
                if (window.ReactantsFormDirty) {
                    ReactantsFormDirty.captureBaseline();
                }
            })
            .catch(error => {
                console.error('Error fetching reaction data:', error);
                // Optionally, handle the error by displaying a message to the user
            });
    }
    subsystemList = await updateSubsystems();
    setupTooltips();

});
async function updateSubsystems() {
    if (subsystemList.length === 0) {
        try {
            window.scrollTo(0, 0);

            const response = await fetch(getVMHsubsystems, {
                method: 'GET',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken
                }
            });

            const data = await response.json();

            if (data.error) {
                showErrorModal(data.message);
                return []; // Return an empty list if there's an error
            } else {
                hidemodal();
                return data.subsystem_list; // Return the fetched list
            }
        } catch (error) {
            console.error('Error fetching subsystems:', error);
            return []; // Return an empty list in case of failure
        }
    }
    else {
        return subsystemList;
    }
}
function createnewreaction() {
    // Get the input element by its ID
    const resetbutton = document.getElementById('ResetButton');
    // Add an event listener to the input element to handle the click event
    resetbutton.addEventListener('click', function (event) {
        // Prevent the default form submission behavior
        event.preventDefault();

        // Redirect to the homepage
        window.location.href = window.location.origin;
    });
}

function setupdate(){
    const urlParams = new URLSearchParams(window.location.search);
    const action = urlParams.get('action');
    if (action === 'edit') {
        document.getElementById('submitBtn-form').childNodes[2].nodeValue = 'Update Reaction';
    }
    // Hide "Save Reaction" in edit mode (redundant for an already-saved reaction).
    if (typeof updateEditModeUI === 'function') {
        updateEditModeUI();
    }
}

/**
 * Responsive column header labels
 * Shrinks font or switches to abbreviations when columns are narrow.
 */
function initResponsiveHeaders() {
    const headers = document.querySelectorAll('.reactant-grid-head');
    if (!headers.length) return;

    const updateHeaders = () => {
        headers.forEach(header => {
            const cols = header.querySelectorAll('.col[data-full]');
            cols.forEach(col => {
                const full = col.dataset.full;
                const mid = col.dataset.mid;
                const short = col.dataset.short;
                const colWidth = col.offsetWidth;

                // Determine which label to show based on available width
                let label = full;
                if (colWidth < 90) {
                    label = short;
                } else if (colWidth < 130) {
                    label = mid;
                }

                // For compartment column, update the inner .col-label span
                const labelEl = col.querySelector('.col-label');
                if (labelEl) {
                    labelEl.textContent = label;
                } else {
                    // Simple text column
                    col.textContent = label;
                }
            });
        });
    };

    // Run once on load and on resize
    updateHeaders();
    window.addEventListener('resize', updateHeaders);
    // Also update when panel is resized (MutationObserver fallback)
    const container = document.getElementById('workspacePanels');
    if (container && window.MutationObserver) {
        const observer = new MutationObserver(() => {
            requestAnimationFrame(updateHeaders);
        });
        observer.observe(container, { attributes: true, subtree: true, attributeFilter: ['style'] });
    }
}


