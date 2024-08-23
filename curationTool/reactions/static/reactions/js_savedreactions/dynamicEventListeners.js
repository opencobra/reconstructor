function attachDynamicEventListeners() {
    // Attach a single event listener to the parent container that holds all your dynamic content
    const modalList = document.getElementById('modalReactionsList');

    modalList.addEventListener('click', function(e) {
        // Handle the "Add" button clicks using event delegation
        if (e.target && e.target.matches('.add-reference, .add-ext-link, .add-comment')) {
            const type = e.target.classList.contains('add-reference') ? 'reference' :
                         e.target.classList.contains('add-ext-link') ? 'ext-link' : 'comment';
            const parentSection = e.target.parentNode;
            const reactionId = e.target.getAttribute('data-reaction-id');
            const newItem = document.createElement('div');
            newItem.className = `${type}-item`;
            const items = parentSection.querySelectorAll(`.${type}-item`);
            const newIndex = items.length; // Calculate new index based on existing items
            newItem.setAttribute('data-reaction-id', reactionId);
            newItem.setAttribute('data-index', newIndex);
            newItem.innerHTML = type === 'ext-link' ? createExtLinkSelect({}, reactionId, newIndex) : '';
            newItem.innerHTML += type === 'reference' ? createRefSelect({}, reactionId, newIndex) : '';
            newItem.innerHTML += `
                <input type="text" class="${type}-input" placeholder="Enter ${type}" value="" data-reaction-id="${reactionId}" data-index="${newIndex}">
                <button class="remove-${type}" data-reaction-id="${reactionId}" data-index="${newIndex}">Remove</button>
            `;
            parentSection.insertBefore(newItem, e.target);
        }

        // Handle the "Remove" button clicks using event delegation
        if (e.target && e.target.matches('.remove-reference, .remove-ext-link, .remove-comment,.remove-gene-info')) {
            e.target.parentElement.remove();
        }
    });
}  