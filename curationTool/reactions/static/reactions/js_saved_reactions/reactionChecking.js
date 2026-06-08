
var checkedReactions = []; // Store IDs of checked reactions
window.checkedReactions = checkedReactions; // Expose to window for graph visualization

// Use event delegation so this survives table rebuilds (e.g. after sort)
document.addEventListener('change', function (e) {
    if (!e.target.classList.contains('reaction-checkbox')) return;
    const reactionId = e.target.getAttribute('data-reaction-id');
    if (e.target.checked) {
        if (!checkedReactions.includes(reactionId)) {
            checkedReactions.push(reactionId);
        }
    } else {
        const index = checkedReactions.indexOf(reactionId);
        if (index > -1) {
            checkedReactions.splice(index, 1);
        }
    }
    window.checkedReactions = checkedReactions;
});

document.getElementById("applyMassConfidence").addEventListener("click", async function () {
    let selectedReactions = document.querySelectorAll(".reaction-checkbox:checked");
    let confidenceScore = document.getElementById("massEditConfidence").value;

    if (!confidenceScore) {
        alert("Please select a confidence score.");
        return;
    }

    if (selectedReactions.length === 0) {
        alert("No reactions selected.");
        return;
    }

    let reactionIds = [];
    selectedReactions.forEach(checkbox => {
        reactionIds.push(checkbox.dataset.reactionId);
    });

    try {
        let response = await fetch(updateConfidenceScoresURL, {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
                "X-CSRFToken": csrfToken
            },
            body: JSON.stringify({
                user_id: userID,
                reaction_ids: reactionIds,
                confidence_score: confidenceScore
            })
        });

        let data = await response.json();
        if (data.status === "success") {
            alert("Confidence scores updated successfully!");
        } else {
            alert("Error updating confidence scores: " + data.message);
        }
    } catch (error) {
        console.error("Error:", error);
        alert("An error occurred while updating confidence scores.");
    }
});
