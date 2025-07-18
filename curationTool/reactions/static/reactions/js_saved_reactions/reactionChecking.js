
var checkedReactions = []; // Store IDs of checked reactions

document.addEventListener('DOMContentLoaded', function () {
    const checkboxes = document.querySelectorAll('.reaction-checkbox');
    checkboxes.forEach(function(checkbox) {
        checkbox.addEventListener('change', function() {
            const reactionId = this.getAttribute('data-reaction-id');
            if (this.checked) {
                // Add the reaction ID to the array
                if (!checkedReactions.includes(reactionId)) {
                    checkedReactions.push(reactionId);
                }
            } else {
                // Remove the reaction ID from the array
                const index = checkedReactions.indexOf(reactionId);
                if (index > -1) {
                    checkedReactions.splice(index, 1);
                }
            }
        });
    });
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
