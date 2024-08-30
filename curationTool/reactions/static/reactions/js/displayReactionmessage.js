function displayReactionMessage(reactionData) {
    // Clear the previous message
    let messageContainer = document.getElementById('reactionFoundMessage');

    messageContainer.innerHTML = '';

    if (reactionData.vmh_found) {
        const reactionFoundMessage = reactionData.vmh_found_similar ?
            'Similar Reaction found at VMH:' : 'Exact Reaction found at VMH:';

        const messageText = document.createTextNode(reactionFoundMessage);
        messageContainer.appendChild(messageText);
        messageContainer.appendChild(document.createElement('br'));

        const link = document.createElement('a');
        link.setAttribute('href', reactionData.vmh_url.trim().replace(/^"|"$/g, ''));
        link.setAttribute('target', '_blank');
        link.textContent = reactionData.vmh_url;
        messageContainer.appendChild(link);
    } else {
        const messageText = document.createTextNode('Reaction not found at VMH');
        messageContainer.appendChild(messageText);
    }

    // Display the message
    messageContainer.style.display = 'block';
}
