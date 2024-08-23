function loadAtomMappingDiv(reactionData) {
    let messageContainer = document.getElementById('reactionFoundMessage');
    let contentDiv = document.querySelector('.content-div[name="atommapping-div"]');
    messageContainer.innerHTML = ''; // Clear previous messages

    if (reactionData.error) {
        let formattedErrorMessage = reactionData.error.replace(/\n/g, '<br>');
        let errorMessage = formattedErrorMessage;
        $('#error-message').html(errorMessage);
        $('#error-modal').modal('show');
        window.scrollTo(0, 0);
    } else {
        let reactionImage = contentDiv.getElementsByClassName('reaction-image')[0];
        let imageClass = reactionImage.className;
        let reactionImageName = reactionImage.name;
        let reactionImageId = reactionImage.id;
        
        // Remove the element and ensure it's done before continuing
        reactionImage.remove();
        
        // Create a new image element with the previous attributes
        let newReactionImage = document.createElement('img');
        newReactionImage.className = imageClass;
        newReactionImage.id = reactionImageId;
        newReactionImage.name = reactionImageName;  // Set the name directly on the new element
    
        // Generate cache buster and set the new image source
        let cacheBuster = new Date().getTime() + "_" + Math.random();
        newReactionImage.src = MEDIA_URL + reactionData.visualization[0] + '?v=' + cacheBuster; // Assuming data.visualization[0] contains the image path
    
        // Append the new image to the contentDiv
        contentDiv.appendChild(newReactionImage);
    
        // Apply the blowup effect to the newly created image
        newReactionImage.addEventListener('load', function() {
            $("#" + reactionImageId).blowup({
                "width": 300,            // Custom width of the lens
                "height": 300,            // Custom height of the lens
                "border" : "6px solid #f2711c",
                "scale" : 2.3 
            });
        });
    
        // Ensure the contentDiv is visible
        contentDiv.style.display = 'block';
    }
}
