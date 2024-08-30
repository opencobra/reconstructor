document.addEventListener('DOMContentLoaded', function () {
    const organTagsContainer = document.getElementById('organTags');
    const organDropdown = document.getElementById('organDropdown');

    const organList_old = [
        "Adipocytes", "Agland", "Brain", "Brain", "Breast", 
        "Brain", "Brain", "Cervix", "Brain", "Colon", "sIEC", 
        "Uterus", "Testis", "Lung", "Uterus", "Gall", "Heart", 
        "Brain", "Brain", "Kidney", "Liver", "Lung", 
        "Brain", "Ovary", "Pancreas", "Pthyroidgland", 
        "Prostate", "Colon", "Retina", "Testis", "Muscle", 
        "Skin", "sIEC", "Scord", "Spleen", "Stomach", 
        "Testis", "Thyroidgland", "Urinarybladder", "Cervix"
    ];
    
    const organList = [...new Set(organList_old)];

    // Placeholder text for the div
    organTagsContainer.setAttribute('data-placeholder', 'Type or select organs...');

    // Show dropdown when the input container is clicked
    organTagsContainer.addEventListener('click', function() {
        renderDropdown(organList, organDropdown);
        organDropdown.style.display = 'block'; // Show the dropdown
    });

    // Handle typing within the div
    organTagsContainer.addEventListener('keydown', function(e) {
        if (e.key === 'Enter' || e.key === ',') {
            e.preventDefault(); // Prevent the default action (adding the character)
            const tagText = organTagsContainer.textContent.trim();
            if (tagText) {
                addTag(tagText);
            }
            organTagsContainer.textContent = ''; // Clear the input area
        } else if (e.key === 'Backspace' && !organTagsContainer.textContent.trim()) {
            // If Backspace is pressed and the input area is empty, delete the last tag
            const tags = organTagsContainer.querySelectorAll('.tag');
            if (tags.length > 0) {
                const lastTag = tags[tags.length - 1];
                organTagsContainer.removeChild(lastTag);
            }
        }
    });

    // Add event listener to hide the dropdown when clicking outside
    document.addEventListener('click', function(event) {
        if (!organDropdown.contains(event.target) && !organTagsContainer.contains(event.target)) {
            organDropdown.style.display = 'none';
        }
    });

    // Add tag within the tags container
    function addTag(tagText) {
        if ([...organTagsContainer.children].some(tag => tag.textContent.replace('×', '').trim() === tagText)) return; // Prevent duplicate tags

        const tag = document.createElement('span');
        tag.className = 'tag';
        tag.contentEditable = false; // Make the tag non-editable
        tag.textContent = tagText;

        const closeBtn = document.createElement('span');
        closeBtn.className = 'close-btn';
        closeBtn.textContent = '×';
        closeBtn.addEventListener('click', function() {
            organTagsContainer.removeChild(tag); // Remove the tag when close button is clicked
        });

        tag.appendChild(closeBtn);
        organTagsContainer.appendChild(tag); // Insert the tag inside the tags container

        // Add a space and an empty span to ensure the cursor moves outside the tag
        const spaceNode = document.createTextNode(' ');
        const caretSpan = document.createElement('span');
        caretSpan.innerHTML = '&nbsp;'; // Invisible space to ensure caret positioning

        organTagsContainer.appendChild(spaceNode);
        organTagsContainer.appendChild(caretSpan);

        // Move the cursor to the end of the contenteditable div
        const range = document.createRange();
        const selection = window.getSelection();
        range.setStartAfter(caretSpan);
        range.collapse(true);
        selection.removeAllRanges();
        selection.addRange(range);

        organTagsContainer.focus(); // Refocus the input area
    }

    // Render the dropdown options
    function renderDropdown(items, dropdown) {
        dropdown.innerHTML = ''; // Clear the dropdown content
        items.forEach(item => {
            const element = document.createElement('div');
            element.textContent = item;
            element.addEventListener('click', function() {
                addTag(item); // Add the selected item as a tag within the input
                dropdown.style.display = 'none'; // Hide the dropdown after selection
            });
            dropdown.appendChild(element);
        });
    }
});




function DisplayTag(tags) {
    if (!tags) {
        return;
    }
    const organTagsContainer = document.getElementById('organTags');

    let str = tags;
    let result = str
    .slice(1, -1)  // Remove the surrounding brackets
    .split(',')     // Split the string by commas
    .map(item => item.replace(/[^\w\s]/g, ''));
    if (result.length === 1 && result[0] === '') {
        return; // Prevent empty tags
    }
    result.forEach(tagText => {
        if ([...organTagsContainer.children].some(tag => tag.textContent.replace('×', '').trim() === tagText)) {
            return; // Prevent duplicate tags
        }

        const tag = document.createElement('span');
        tag.className = 'tag';
        tag.textContent = tagText;

        const closeBtn = document.createElement('span');
        closeBtn.className = 'close-btn';
        closeBtn.textContent = '×';
        closeBtn.style.marginLeft = '5px';
        closeBtn.style.cursor = 'pointer';

        closeBtn.addEventListener('click', function() {
            organTagsContainer.removeChild(tag);
        });

        tag.appendChild(closeBtn);
        organTagsContainer.appendChild(tag);
    });
}

