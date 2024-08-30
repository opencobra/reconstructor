document.addEventListener('DOMContentLoaded', function () {

    document.getElementById('viewsavedreactions-button').addEventListener('click', function () {
        userId = sessionStorage.getItem('userID');
        if (userId) {
            window.location.href = '/saved_reactions';
        }
        else{
            var errorMessage = 'Please login to view saved reactions.';
            showErrorModal(errorMessage);
        }
    });

});