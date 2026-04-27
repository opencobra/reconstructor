function displayReactionMessage(reactionData) {
    const messageContainer = document.getElementById('reactionFoundMessage');
    if (!messageContainer) {
        return;
    }

    messageContainer.innerHTML = '';

    const hasMatch = Boolean(reactionData && reactionData.vmh_found);
    const isSimilar = Boolean(reactionData && reactionData.vmh_found_similar);
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
    badge.className = isSimilar ? 'vmh-badge vmh-badge--similar' : 'vmh-badge vmh-badge--exact';

    if (hasLink) {
        badge.href = cleanedUrl;
        badge.target = '_blank';
        badge.rel = 'noopener noreferrer';
    }

    const icon = document.createElement('i');
    icon.className = isSimilar ? 'fas fa-code-branch' : 'fas fa-check-circle';
    icon.setAttribute('aria-hidden', 'true');

    const label = document.createElement('span');
    label.textContent = isSimilar
        ? `VMH similar reaction: ${reactionId}`
        : `VMH exact reaction: ${reactionId}`;

    badge.appendChild(icon);
    badge.appendChild(label);
    messageContainer.appendChild(badge);
}
