
var checkedReactions = []; // Store IDs of checked reactions
window.checkedReactions = checkedReactions; // Expose to window for graph visualization

// Keep the toolbar's selection count in sync so it's clear that the bulk
// actions apply to whatever is currently selected.
function updateSelectionCount() {
    const countEl = document.getElementById('selectionCount');
    if (!countEl) return;
    const n = checkedReactions.length;
    countEl.dataset.count = String(n);
    countEl.textContent =
        n === 0 ? 'No reactions selected' : `${n} reaction${n === 1 ? '' : 's'} selected`;
}

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
    updateSelectionCount();
});

document.addEventListener('DOMContentLoaded', updateSelectionCount);

// Build the CS badge markup for a (server-confirmed) confidence score. Mirrors
// the server-side template so updated rows match a freshly rendered page.
function renderConfidenceBadge(score) {
    const valid = ['1', '2', '3', '4'];
    if (score && valid.includes(String(score))) {
        return `<span class="cs-badge cs-badge--${score}">${score}</span>`;
    }
    return '<span class="cs-badge cs-badge--none">–</span>';
}

// Refresh just the CS cell of each updated row, in place, keeping row order and
// everything else untouched. `updated` maps reaction_id -> normalized score.
function applyConfidenceUpdates(updated) {
    Object.entries(updated).forEach(([reactionId, score]) => {
        const checkbox = document.querySelector(
            `.reaction-checkbox[data-reaction-id="${reactionId}"]`);
        const row = checkbox && checkbox.closest('tr');
        const cell = row && row.querySelector('.cs-col');
        if (cell) {
            cell.innerHTML = renderConfidenceBadge(score);
        }
    });
}

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
            // Reflect the server's persisted values in the affected rows only.
            applyConfidenceUpdates(data.updated || {});
            const n = Object.keys(data.updated || {}).length;
            if (typeof showToast === 'function') {
                showToast(`Confidence score updated for ${n} reaction${n === 1 ? '' : 's'}`);
            }
        } else {
            alert("Error updating confidence scores: " + data.message);
        }
    } catch (error) {
        console.error("Error:", error);
        alert("An error occurred while updating confidence scores.");
    }
});
