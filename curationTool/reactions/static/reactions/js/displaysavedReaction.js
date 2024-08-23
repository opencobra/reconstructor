document.addEventListener('DOMContentLoaded', function () {
    // Function to redirect to the saved reactions page
    if (sessionStorage.getItem('userID') !== null) {
    document.getElementById('viewsavedreactions-button').addEventListener('click', function () {
        window.location.href = '/saved_reactions';
    });
}
else {
    document.getElementById('viewsavedreactions-button').addEventListener('click', function () {
        var errorMessage = 'Please login to view saved reactions.';
        showErrorModal(errorMessage);    });
}
    
});

// Get the modal
// Get the modal
// Get the modal
// var savedReactionsModal = document.getElementById("savedReactionsModal");

// // Get the button that opens the modal
// var savedReactionsBtn = document.getElementById("viewsavedreactions-button");

// // Get the <span> element that closes the modal
// var closeSavedReactions = document.getElementsByClassName("close-saved-reactions")[0];

// // When the user clicks the button, open the modal and load the content
// savedReactionsBtn.onclick = function() {
//     savedReactionsModal.style.display = "block";

//     // Load content dynamically via AJAX or Fetch API
//     fetch('saved_reactions')
//         .then(response => response.text())
//         .then(data => {
//             document.getElementById("savedReactionsModalBody").innerHTML = data;

//             // Option 1: Dynamically load the required JavaScript files using static paths
//             loadScript("{% static 'reactions/js_savedreactions/reactionChecking.js' %}");
//             loadScript("{% static 'reactions/js_savedreactions/dynamicEventListeners.js' %}");
//             loadScript("{% static 'reactions/js_savedreactions/reactionViewAndDelete.js' %}");
//             loadScript("{% static 'reactions/js_savedreactions/queryai.js' %}");
//             loadScript("{% static 'reactions/js_savedreactions/utilities.js' %}");
//             loadScript("{% static 'reactions/js_savedreactions/handleAdd2VMH.js' %}");

//             // Option 2: Dynamically load the required JavaScript files using data attributes
//             loadScript(savedReactionsBtn.getAttribute('data-reaction-checking'));
//             loadScript(savedReactionsBtn.getAttribute('data-dynamic-listeners'));
//             loadScript(savedReactionsBtn.getAttribute('data-view-delete'));
//             loadScript(savedReactionsBtn.getAttribute('data-query-ai'));
//             loadScript(savedReactionsBtn.getAttribute('data-utilities'));
//             loadScript(savedReactionsBtn.getAttribute('data-handle-vmh'));
//         })
//         .catch(error => console.error('Error loading the modal content:', error));
// }

// // When the user clicks on <span> (x), close the modal
// closeSavedReactions.onclick = function() {
//     savedReactionsModal.style.display = "none";
// }

// // When the user clicks anywhere outside of the modal, close it
// window.onclick = function(event) {
//     if (event.target == savedReactionsModal) {
//         savedReactionsModal.style.display = "none";
//     }
// }

// // Function to dynamically load scripts
// function loadScript(src) {
//     return new Promise(function(resolve, reject) {
//         var script = document.createElement('script');
//         script.src = src;
//         script.onload = function() {
//             resolve();
//         };
//         script.onerror = function() {
//             reject(new Error("Script load error for " + src));
//         };
//         document.head.append(script);
//     });
// }
