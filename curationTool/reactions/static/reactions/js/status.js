function buildReactionNameElement(name, description) {
    if (!name) {
        return null;
    }

    const reactionName = document.createElement('span');
    reactionName.className = 'reaction-name';
    reactionName.textContent = name;
    if (description) {
        reactionName.setAttribute('data-tooltip-this', description);
    }

    return reactionName;
}

function renderStatusBadge(options) {
    const {
        dotClass,
        iconClass,
        label,
        reactionName,
        reactionDescription,
    } = options;

    const statusElement = document.getElementById('statusTitle');
    if (!statusElement) {
        console.error('Status element not found in the DOM');
        return;
    }

    const badge = document.createElement('span');
    badge.className = `status-badge ${dotClass}`;

    const icon = document.createElement('i');
    icon.className = iconClass;
    icon.setAttribute('aria-hidden', 'true');
    badge.appendChild(icon);

    const labelSpan = document.createElement('span');
    labelSpan.className = 'status-label';
    labelSpan.textContent = label;
    badge.appendChild(labelSpan);

    const reactionNameElement = buildReactionNameElement(reactionName, reactionDescription);
    if (reactionNameElement) {
        const separator = document.createElement('span');
        separator.className = 'status-separator';
        separator.textContent = '·';
        badge.appendChild(separator);
        badge.appendChild(reactionNameElement);
    }

    statusElement.replaceChildren(badge);
}

function setLoggedInStatusBasedOnUrl(reactionData) {
    if (!reactionData || reactionData === '') {
        renderStatusBadge({
            dotClass: 'dot-red',
            iconClass: 'fas fa-plus-circle',
            label: 'Creating reaction',
        });
        return;
    }

    const urlParams = new URLSearchParams(window.location.search);
    const reactionId = urlParams.get('reaction_id');
    const action = urlParams.get('action');

    const reactionName = reactionData.name || 'Unnamed reaction';
    const reactionDescription = reactionData.description || '';

    if (action !== 'edit' && reactionId !== null) {
        renderStatusBadge({
            dotClass: 'dot-green',
            iconClass: 'fas fa-eye',
            label: 'Viewing reaction',
            reactionName,
            reactionDescription,
        });
    } else if (reactionId === null && action !== 'edit') {
        renderStatusBadge({
            dotClass: 'dot-red',
            iconClass: 'fas fa-plus-circle',
            label: 'Creating reaction',
        });
    } else {
        renderStatusBadge({
            dotClass: 'dot-orange',
            iconClass: 'fas fa-pencil-alt',
            label: 'Editing reaction',
            reactionName,
            reactionDescription,
        });
    }
}

function setLoggedOutStatusBasedOnUrl() {
    renderStatusBadge({
        dotClass: 'dot-grey',
        iconClass: 'fas fa-clock',
        label: 'Idle',
    });
}

