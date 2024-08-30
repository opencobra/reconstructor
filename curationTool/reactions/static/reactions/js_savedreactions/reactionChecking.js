
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