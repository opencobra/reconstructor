function displayReactionMessage(reactionData) {
    const messageContainer = document.getElementById('reactionFoundMessage');
    if (!messageContainer) {
        return;
    }

    messageContainer.innerHTML = '';

    const hasMatch = Boolean(reactionData && reactionData.vmh_found);
    const isSimilar = Boolean(reactionData && reactionData.vmh_found_similar);
    const addedViaConstructor = Boolean(
        reactionData && reactionData.vmh_added_via_constructor);
    const rawUrl = reactionData && reactionData.vmh_url ? reactionData.vmh_url : '';
    const cleanedUrl = rawUrl.trim().replace(/^"|"$/g, '');

    if (!hasMatch) {
        const emptyBadge = document.createElement('span');
        emptyBadge.className = 'vmh-badge vmh-badge--none';
        emptyBadge.innerHTML =
            '<i class="fas fa-search" aria-hidden="true"></i><span>No VMH match found</span>';
        messageContainer.appendChild(emptyBadge);
        return;
    }

    const reactionId = cleanedUrl.split('/').filter(Boolean).pop() || 'Open in VMH';
    const hasLink = cleanedUrl !== '';

    const badge = hasLink ? document.createElement('a') : document.createElement('span');
    let modifierClass = 'vmh-badge--exact';
    if (addedViaConstructor) {
        modifierClass = 'vmh-badge--added';
    } else if (isSimilar) {
        modifierClass = 'vmh-badge--similar';
    }
    badge.className = `vmh-badge ${modifierClass}`;

    if (hasLink) {
        badge.href = cleanedUrl;
        badge.target = '_blank';
        badge.rel = 'noopener noreferrer';
    }

    const icon = document.createElement('i');
    if (addedViaConstructor) {
        icon.className = 'fas fa-plus-circle';
    } else {
        icon.className = isSimilar ? 'fas fa-code-branch' : 'fas fa-check-circle';
    }
    icon.setAttribute('aria-hidden', 'true');

    const label = document.createElement('span');
    if (addedViaConstructor) {
        label.textContent = 'Reaction added to VMH via constructor';
    } else if (isSimilar) {
        label.textContent = `VMH similar reaction: ${reactionId}`;
    } else {
        label.textContent = `VMH exact reaction: ${reactionId}`;
    }

    badge.appendChild(icon);
    badge.appendChild(label);
    messageContainer.appendChild(badge);
}
