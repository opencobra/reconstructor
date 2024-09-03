var button_to_div = {
    'reactants-button': 'reactants',
    'atommapping-button': 'atommapping',
    'cheminfo-button': 'cheminfo',
    'metinfo-button': 'metaboliteinfo',
    'references-button': 'references',
    'extlinks-button': 'extlinks',
    'reactinfo-button': 'reactioninfo',
    'geneinfo-button': 'gene_info',
    'comments-button': 'comments',
    'reactiontemps-button': 'reaction_temps'
    // 'viewsavedreactions-button' is excluded
};

document.addEventListener('DOMContentLoaded', function() {
    var buttons = document.querySelectorAll('.dynamic-button-side-button');
    buttons.forEach(function(button) {
        // Check if the button id is in the button_to_div object
        if (button.id in button_to_div) {
            button.addEventListener('click', function() {
                const divName = button_to_div[button.id];
                const div = document.getElementsByName(divName + '-div')[0];

                // Toggle the visibility of the clicked div
                if (div.style.display === 'block') {
                    div.style.display = 'none';
                    button.classList.remove('active');
                } else {
                    div.style.display = 'block';
                    button.classList.add('active');
                }
            });
        }
    });
});

function refreshSideButtons() {
    for (const [button, div] of Object.entries(button_to_div)) {
        const buttonElement = document.getElementById(button);
        const divElement = document.getElementsByName(div + '-div')[0];
        if (divElement) {
            if (divElement.style.display === 'block') {
                buttonElement.classList.add('active');
            } else {
                buttonElement.classList.remove('active');
            }
        }
    }
}
document.addEventListener('DOMContentLoaded', function() {
    refreshSideButtons();
    
    // Add event listener to item about-item to go to the about page (in new tab)
    var aboutButton = document.getElementById('item about-item');
    
    // Check if the element exists to prevent errors

    aboutButton.addEventListener('click', function() {  
        // Open the About page in a new tab
        window.open(aboutUrl, '_blank');
    });
});