// Function to create and show the modal
function showEditWarningModal(reactionId) {
    // Create the modal
    let modal = document.createElement("div");
    modal.style.position = "fixed";
    modal.style.zIndex = "1000";
    modal.style.left = "50%";
    modal.style.top = "50%";
    modal.style.transform = "translate(-50%, -50%)";
    modal.style.width = "300px";
    modal.style.padding = "20px";
    modal.style.backgroundColor = "white";
    modal.style.boxShadow = "0 4px 8px rgba(0, 0, 0, 0.2)";
    modal.style.borderRadius = "8px";
    modal.style.textAlign = "center";

    // Modal content
    let modalContent = document.createElement("p");
    modalContent.innerText = "⚠️ Warning: You are currently editing a reaction";
    modal.appendChild(modalContent);

    // Close button
    let closeButton = document.createElement("button");
    closeButton.innerText = "Close";
    closeButton.style.marginTop = "10px";
    closeButton.style.padding = "5px 10px";
    closeButton.style.border = "none";
    closeButton.style.backgroundColor = "#007bff";
    closeButton.style.color = "white";
    closeButton.style.borderRadius = "4px";
    closeButton.style.cursor = "pointer";
    closeButton.onclick = function() {
        document.body.removeChild(modal);
    };
    modal.appendChild(closeButton);

    // Append modal to body
    document.body.appendChild(modal);

    // Redirect after showing modal
    setTimeout(function() {
        window.location.href = `/?reaction_id=${reactionId}&action=edit`;
    }, 3000);  // Adjust the delay as needed
}

// Attach event listener to buttons with the class 'view-btn'
document.querySelectorAll('.view-btn').forEach(function(button, index) {
    button.addEventListener('click', function() {
        const reactionData = reactions[index];
        const reactionId = reactionData.pk;
        showEditWarningModal(reactionId);
    });
});


function confirmDelete(reactionName) {
    return confirm('Are you sure you want to delete reaction ' + reactionName + '?');
}
