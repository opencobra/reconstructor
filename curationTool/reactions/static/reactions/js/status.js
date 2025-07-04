function setLoggedInStatusBasedOnUrl(reactionData) {
    let status = '';
    let dotClass = '';

    if (reactionData === '') {
        status = "Creating Reaction";
        dotClass = "dot-red"; 
    } else {
        const urlParams = new URLSearchParams(window.location.search);
        const reactionId = urlParams.get('reaction_id');
        const action = urlParams.get('action');
        let reactionName = reactionData.name;
        let reactionDescription = reactionData.description;

        if (action !== "edit" && reactionId !== null) {
            status = `Viewing reaction <br>
                      <span class="reaction-name" data-tooltip-this="${reactionDescription}">${reactionName}</span>`;
            dotClass = "dot-green";
        } else if (reactionId === null && action !== "edit") {
            status = "Creating Reaction";
            dotClass = "dot-red"; 
        } else {
            status = `Editing reaction 
                      <span class="reaction-name" data-tooltip-this="${reactionDescription}">${reactionName}</span>`;
            dotClass = "dot-orange";
        }
    }

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

