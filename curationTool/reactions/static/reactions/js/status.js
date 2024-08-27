
function setLoggedInStatusBasedOnUrl() {
    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');

    const action = urlParams.get('action');
    // 2. Determine the status for logged-in users
    let status = "Creating reaction"; // Default for logged-in users
    let dotClass = "dot-blue"; // Default color for "Creating reaction"

    if (action != "edit" && reactionId != null) {
        status = "Viewing reaction";
        dotClass = "dot-green"; // Green dot for "Viewing reaction"
    } else if (reactionId == null && action != "edit") {
        status = "creating reaction";
        dotClass = "dot-red"; // Orange dot for "Editing reaction"
    }
    else {
        status = "Editing reaction";
        dotClass = "dot-orange"; // Orange dot for "Editing reaction"
    }

    // 3. Display the status in the div with id="statusTitle"
    const statusElement = document.getElementById('statusTitle');
    if (statusElement) {
        statusElement.innerHTML = `Status: <span class="status-dot-top ${dotClass}"></span> ${status}`;
    } else {
        console.error("Status element not found in the DOM");
    }
}

function setLoggedOutStatusBasedOnUrl() {
    // 1. Determine the status for logged-out users
    const status = "Idle"; // Default for logged-out users
    const dotClass = "dot-grey"; // Grey dot for "Idle"

    // 2. Display the status in the div with id="statusTitle"
    const statusElement = document.getElementById('statusTitle');
    if (statusElement) {
        statusElement.innerHTML = `Status: <span class="status-dot-top ${dotClass}"></span> ${status}`;
    } else {
        console.error("Status element not found in the DOM");
    }
}

