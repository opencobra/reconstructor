document.getElementById('submitBtn-form').addEventListener('click', function(event) {
    event.preventDefault(); // Prevent the default form submission

    // Trigger the form's submit event
    document.getElementById('reactionForm').requestSubmit();
});

// Keep loader outside hidden modal containers so it is always visible.
(function ensureGlobalLoadingIndicator() {
    const loadingIndicator = document.getElementById('loadingIndicator');
    if (loadingIndicator && loadingIndicator.parentElement !== document.body) {
        document.body.appendChild(loadingIndicator);
    }
})();

var ReactantsFormDirty = (function () {
    const CLEAN_TITLE = 'No reaction changes to update';
    const EMPTY_LABEL = 'empty';
    let baseline = null;
    let currentDirty = false;
    let evaluateTimer = null;
    let observer = null;

    function params() {
        return new URLSearchParams(window.location.search);
    }

    function isEditMode() {
        return params().get('action') === 'edit';
    }

    function text(value) {
        return String(value == null ? '' : value).trim();
    }

    function normalizeStoich(value) {
        const raw = text(value);
        if (!raw) return '';
        const numberValue = Number(raw);
        return Number.isFinite(numberValue) ? String(numberValue) : raw;
    }

    function normalizeDirectionValue(value) {
        if (typeof normalizeReactionDirection === 'function') {
            return normalizeReactionDirection(value);
        }
        return text(value).toLowerCase() === 'bidirectional' ? 'bidirectional' : 'forward';
    }

    function getIdentifierValue(group, fieldName) {
        const autocomplete = group.querySelector('.autocomplete-container');
        if (autocomplete) {
            const hidden = autocomplete.querySelector('input[type="hidden"]');
            const visible = autocomplete.querySelector('.autocomplete-input');
            return text((hidden && hidden.value) || (visible && visible.value));
        }

        const fileInput = group.querySelector('.cell-identifier input[type="file"]');
        if (fileInput && fileInput.files && fileInput.files.length > 0) {
            return text(fileInput.files[0].name);
        }

        const input = group.querySelector('.cell-identifier input[name="' + fieldName + '"]');
        return text(input && input.value);
    }

    function readRows(containerId, side) {
        const container = document.getElementById(containerId);
        if (!container) return [];
        return Array.from(container.querySelectorAll('.inputs-group')).map(function (group) {
            const stoichInput = group.querySelector('input[name="' + side.stoich + '"]');
            const compartmentSelect = group.querySelector('select[name="' + side.compartment + '"]');
            const typeSelect = group.querySelector('select[name="' + side.type + '"]');
            return {
                identifier: getIdentifierValue(group, side.field),
                stoichiometry: normalizeStoich(stoichInput && stoichInput.value),
                compartment: text(compartmentSelect && compartmentSelect.value),
                inputType: text(typeSelect && typeSelect.value),
            };
        });
    }

    function readOrgans() {
        const container = document.getElementById('organTags');
        if (!container) return [];
        return Array.from(container.querySelectorAll('.tag'))
            .map(function (tag) {
                return text(tag.firstChild ? tag.firstChild.textContent : tag.textContent.replace('×', ''));
            })
            .filter(Boolean);
    }

    function collectState() {
        const directionEl = document.getElementById('reactionDirection');
        const subsystemEl = document.getElementById('subsystemField');
        return {
            substrates: readRows('substratesDiv', {
                field: 'substrates',
                stoich: 'subs_sch',
                compartment: 'subs_comps',
                type: 'substrates_type',
            }),
            products: readRows('productsDiv', {
                field: 'products',
                stoich: 'prod_sch',
                compartment: 'prod_comps',
                type: 'products_type',
            }),
            direction: normalizeDirectionValue(directionEl ? directionEl.value : 'forward'),
            subsystem: text(subsystemEl && subsystemEl.value),
            organs: readOrgans(),
        };
    }

    function stateKey(state) {
        return JSON.stringify(state || null);
    }

    function setOpenReactionDirty(value) {
        if (window.OpenReactionState) {
            OpenReactionState.setDirty(undefined, !!value);
        }
    }

    function updateButton() {
        const submitBtn = document.getElementById('submitBtn-form');
        if (!submitBtn) return;

        if (submitBtn.dataset.submitting === 'true') {
            submitBtn.disabled = true;
            submitBtn.classList.remove('reactants-update-disabled');
            submitBtn.removeAttribute('aria-disabled');
            submitBtn.removeAttribute('title');
            return;
        }

        if (!isEditMode()) {
            submitBtn.disabled = false;
            submitBtn.classList.remove('reactants-update-disabled');
            submitBtn.removeAttribute('aria-disabled');
            submitBtn.removeAttribute('title');
            return;
        }

        if (currentDirty) {
            submitBtn.disabled = false;
            submitBtn.classList.remove('reactants-update-disabled');
            submitBtn.setAttribute('aria-disabled', 'false');
            submitBtn.removeAttribute('title');
        } else {
            submitBtn.disabled = true;
            submitBtn.classList.add('reactants-update-disabled');
            submitBtn.setAttribute('aria-disabled', 'true');
            submitBtn.setAttribute('title', CLEAN_TITLE);
        }
    }

    function evaluateNow() {
        evaluateTimer = null;
        if (!isEditMode()) {
            baseline = null;
            currentDirty = false;
            updateButton();
            return false;
        }
        if (!baseline) {
            currentDirty = false;
            updateButton();
            return false;
        }
        currentDirty = stateKey(collectState()) !== stateKey(baseline);
        setOpenReactionDirty(currentDirty);
        updateButton();
        return currentDirty;
    }

    function scheduleEvaluate() {
        if (evaluateTimer) {
            clearTimeout(evaluateTimer);
        }
        evaluateTimer = setTimeout(evaluateNow, 0);
    }

    function captureBaseline() {
        if (!isEditMode()) {
            baseline = null;
            currentDirty = false;
            setOpenReactionDirty(false);
            updateButton();
            return;
        }
        baseline = collectState();
        currentDirty = false;
        setOpenReactionDirty(false);
        updateButton();
    }

    function labelValue(value) {
        const normalized = text(value);
        return normalized || EMPTY_LABEL;
    }

    function rowTerm(row) {
        if (!row) return EMPTY_LABEL;
        const stoich = row.stoichiometry && row.stoichiometry !== '1'
            ? row.stoichiometry + ' '
            : '';
        const identifier = labelValue(row.identifier);
        const compartment = row.compartment ? '[' + row.compartment + ']' : '';
        return stoich + identifier + compartment;
    }

    function equationFor(state) {
        const substrates = state.substrates.length ? state.substrates.map(rowTerm).join(' + ') : EMPTY_LABEL;
        const products = state.products.length ? state.products.map(rowTerm).join(' + ') : EMPTY_LABEL;
        const arrow = state.direction === 'bidirectional' ? '<=>' : '->';
        return substrates + ' ' + arrow + ' ' + products;
    }

    function fieldLabel(key) {
        const labels = {
            identifier: 'metabolite identifier',
            stoichiometry: 'stoichiometry',
            compartment: 'compartment',
            inputType: 'input type',
        };
        return labels[key] || key;
    }

    function pushDiff(diffs, label, beforeValue, afterValue) {
        diffs.push({
            label: label,
            before: labelValue(beforeValue),
            after: labelValue(afterValue),
        });
    }

    function rowKey(row) {
        return ['identifier', 'stoichiometry', 'compartment', 'inputType']
            .map(function (key) { return row ? row[key] : ''; })
            .join('\u0001');
    }

    function rowsEqual(beforeRow, afterRow) {
        return rowKey(beforeRow) === rowKey(afterRow);
    }

    function alignRows(beforeRows, afterRows) {
        const removeCost = 1;
        const addCost = 1;
        const changeCost = 1.5;
        const m = beforeRows.length;
        const n = afterRows.length;
        const dp = Array.from({ length: m + 1 }, function () {
            return Array(n + 1).fill(0);
        });
        const op = Array.from({ length: m + 1 }, function () {
            return Array(n + 1).fill(null);
        });

        for (let i = m - 1; i >= 0; i -= 1) {
            dp[i][n] = removeCost + dp[i + 1][n];
            op[i][n] = 'remove';
        }
        for (let j = n - 1; j >= 0; j -= 1) {
            dp[m][j] = addCost + dp[m][j + 1];
            op[m][j] = 'add';
        }

        for (let i = m - 1; i >= 0; i -= 1) {
            for (let j = n - 1; j >= 0; j -= 1) {
                if (rowsEqual(beforeRows[i], afterRows[j])) {
                    dp[i][j] = dp[i + 1][j + 1];
                    op[i][j] = 'match';
                    continue;
                }

                const remove = removeCost + dp[i + 1][j];
                const add = addCost + dp[i][j + 1];
                const change = changeCost + dp[i + 1][j + 1];
                let bestCost = change;
                let bestOp = 'change';

                if (remove < bestCost) {
                    bestCost = remove;
                    bestOp = 'remove';
                }
                if (add < bestCost) {
                    bestCost = add;
                    bestOp = 'add';
                }

                dp[i][j] = bestCost;
                op[i][j] = bestOp;
            }
        }

        const operations = [];
        let i = 0;
        let j = 0;
        while (i < m || j < n) {
            const action = op[i][j];
            if (action === 'match') {
                operations.push({ type: 'match', beforeIndex: i, afterIndex: j });
                i += 1;
                j += 1;
            } else if (action === 'change') {
                operations.push({ type: 'change', beforeIndex: i, afterIndex: j });
                i += 1;
                j += 1;
            } else if (action === 'remove') {
                operations.push({ type: 'remove', beforeIndex: i });
                i += 1;
            } else if (action === 'add') {
                operations.push({ type: 'add', afterIndex: j });
                j += 1;
            } else {
                break;
            }
        }
        return operations;
    }

    function pushRowFieldDiffs(diffs, beforeRow, afterRow, rowLabel) {
        ['identifier', 'stoichiometry', 'compartment', 'inputType'].forEach(function (key) {
            if (beforeRow[key] !== afterRow[key]) {
                pushDiff(diffs, rowLabel + ' ' + fieldLabel(key), beforeRow[key], afterRow[key]);
            }
        });
    }

    function compareRows(diffs, beforeRows, afterRows, singularLabel) {
        alignRows(beforeRows, afterRows).forEach(function (operation) {
            if (operation.type === 'match') {
                return;
            }

            if (operation.type === 'remove') {
                const beforeRow = beforeRows[operation.beforeIndex];
                pushDiff(diffs, singularLabel + ' ' + (operation.beforeIndex + 1), rowTerm(beforeRow), 'removed');
                return;
            }

            if (operation.type === 'add') {
                const afterRow = afterRows[operation.afterIndex];
                pushDiff(diffs, singularLabel + ' ' + (operation.afterIndex + 1), 'not present', rowTerm(afterRow));
                return;
            }

            if (operation.type === 'change') {
                const beforeRow = beforeRows[operation.beforeIndex];
                const afterRow = afterRows[operation.afterIndex];
                const rowLabel = singularLabel + ' ' + (operation.beforeIndex + 1);
                pushRowFieldDiffs(diffs, beforeRow, afterRow, rowLabel);
            }
        });
    }

    function diffStates(beforeState, afterState) {
        const diffs = [];
        compareRows(diffs, beforeState.substrates, afterState.substrates, 'Substrate');
        compareRows(diffs, beforeState.products, afterState.products, 'Product');
        if (beforeState.direction !== afterState.direction) {
            pushDiff(diffs, 'Direction', beforeState.direction, afterState.direction);
        }
        if (beforeState.subsystem !== afterState.subsystem) {
            pushDiff(diffs, 'Subsystem', beforeState.subsystem, afterState.subsystem);
        }
        if (stateKey(beforeState.organs) !== stateKey(afterState.organs)) {
            pushDiff(diffs, 'Organs', beforeState.organs.join(', '), afterState.organs.join(', '));
        }
        return diffs;
    }

    function escapeHtml(value) {
        const div = document.createElement('div');
        div.textContent = value == null ? '' : String(value);
        return div.innerHTML;
    }

    function confirmationHtml() {
        const beforeState = baseline || collectState();
        const afterState = collectState();
        const diffs = diffStates(beforeState, afterState);
        const diffItems = diffs.length
            ? diffs.map(function (diff) {
                return '<li><strong>' + escapeHtml(diff.label) + ':</strong> ' +
                    escapeHtml(diff.before) + ' → ' + escapeHtml(diff.after) + '</li>';
            }).join('')
            : '<li>No Reactants field changes detected.</li>';

        return '' +
            '<div class="reaction-update-summary">' +
                '<p class="reaction-update-warning">This updates the saved reaction and cannot be undone.</p>' +
                '<div class="reaction-update-equations">' +
                    '<div class="reaction-update-label">Before</div>' +
                    '<div class="reaction-update-equation">' + escapeHtml(equationFor(beforeState)) + '</div>' +
                    '<div class="reaction-update-label">After</div>' +
                    '<div class="reaction-update-equation">' + escapeHtml(equationFor(afterState)) + '</div>' +
                '</div>' +
                '<div class="reaction-update-diff-title">Changed fields</div>' +
                '<ul class="reaction-update-diff-list">' + diffItems + '</ul>' +
            '</div>';
    }

    function bind() {
        const container = document.getElementById('reactionFormContainer');
        if (container) {
            container.addEventListener('input', scheduleEvaluate, true);
            container.addEventListener('change', scheduleEvaluate, true);
        }

        const observed = ['substratesDiv', 'productsDiv', 'organTags']
            .map(function (id) { return document.getElementById(id); })
            .filter(Boolean);
        if (window.MutationObserver && observed.length) {
            observer = new MutationObserver(scheduleEvaluate);
            observed.forEach(function (el) {
                observer.observe(el, { childList: true, subtree: true });
            });
        }
        updateButton();
    }

    document.addEventListener('DOMContentLoaded', bind);

    return {
        captureBaseline: captureBaseline,
        collectState: collectState,
        confirmationHtml: confirmationHtml,
        hasChanges: function () { return evaluateNow(); },
        isEditMode: isEditMode,
        refreshButton: updateButton,
        scheduleEvaluate: scheduleEvaluate,
    };
})();

window.ReactantsFormDirty = ReactantsFormDirty;

function normalizeReactionDirection(direction) {
    const value = (direction || '').toString().trim().toLowerCase();
    if (value === 'forward' || value === '->' || value === '=>' || value === '=') {
        return 'forward';
    }
    if (
        value === 'bidirectional' ||
        value === 'reversible' ||
        value === '<=>' ||
        value === '<->' ||
        value === '<-->' ||
        value === '↔' ||
        value === '⇌' ||
        value === 'reverse' ||
        value === 'backward' ||
        value === '<=' ||
        value === '<-'
    ) {
        return 'bidirectional';
    }
    return 'forward';
}

function getMetaboliteIdentifierInput(group) {
    if (!group) return null;
    return group.querySelector(
        '.cell-identifier input[name="substrates"], .cell-identifier input[name="products"]'
    );
}

function reactionRowHasMetabolite(group) {
    if (!group) return false;
    const fileInput = group.querySelector('.cell-identifier input[type="file"]');
    if (fileInput && fileInput.files && fileInput.files.length > 0) {
        return true;
    }
    const identifierInput = getMetaboliteIdentifierInput(group);
    return Boolean(identifierInput && identifierInput.value && identifierInput.value.trim());
}

function getActiveReactionRows(root = document) {
    return Array.from(root.querySelectorAll('.inputs-group')).filter(reactionRowHasMetabolite);
}

function hidemodal(){

    document.getElementById('error-modal').style.display = 'none';

}

document.getElementById('close-button').onclick = function() {
    document.getElementById('error-modal').style.display = 'none';
};

function showIdenticalReactionModal(matches, submitBtn) {
    return new Promise((resolve) => {
        const modal = $('#identicalReactionModal');

        // 1. Build a dynamic message
        let messageText = '';
        if (matches.length > 1) {
            messageText += `Multiple identical reactions have already been saved in your reactions: <br><br>`;
        } else {
            messageText += `An identical reaction has already been saved in your reactions: <br><br>`;
        }

        // 2. Generate a "View" button for each match
        matches.forEach(({ reaction_id, reaction_name }) => {
            messageText += `
                <button 
                  class="ui button dynamic-view-btn" 
                  data-reaction-id="${reaction_id}" 
                  style="margin-bottom: 0.5em; margin-right: 0.5em;"
                >
                  View "${reaction_name}"
                </button>
            `;
        });

        // 3. Prompt user they can create new
        messageText += `<br><br>Click "Create" to proceed with creating a new reaction.`;

        // Insert the HTML into the modal
        document.getElementById('identicalReactionMessage').innerHTML = messageText;

        // 4. Attach handlers for each generated "View" button
        setTimeout(() => {
            const viewButtons = document.querySelectorAll('.dynamic-view-btn');
            viewButtons.forEach(btn => {
                btn.addEventListener('click', function() {
                    const rxnId = this.getAttribute('data-reaction-id');
                    resolve(`view:${rxnId}`);
                    modal.modal('hide');
                });
            });
        }, 0);

        // 5. Keep the existing "Create" button
        document.getElementById('createReactionButton').onclick = function () {
            resolve('create');
            modal.modal('hide');
        };

        // Re-enable submit button if the user closes modal
        modal.modal({
            onHide: function() {
                submitBtn.disabled = false;
            }
        });

        // Show the modal
        modal.modal('show');
    });
}


async function checkIdenticalReaction(formData, loadingIndicator, submitBtn, logTiming) {
    try {
        if (logTiming) logTiming('identicalReaction fetch:start');
        const response = await fetch(identicalReactionUrl, {
            method: 'POST',
            body: formData,
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken
            }
        });
        if (logTiming) logTiming('identicalReaction response headers', { status: response.status });

        const data = await response.json();
        if (logTiming) logTiming('identicalReaction response json parsed', { exists: !!data.exists, status: data.status });

        if (data.status === 'error') {
            showErrorModal(data.message);
            window.scrollTo(0, 0);
            submitBtn.disabled = false;
            loadingIndicator.style.display = 'none';
            return false;
        }

        if (data.exists) {
            loadingIndicator.style.display = 'none';
            if (logTiming) logTiming('identicalReaction duplicate modal shown', { matches: data.matches ? data.matches.length : 0 });
            const userDecision = await showIdenticalReactionModal(data.matches, submitBtn);
            if (logTiming) logTiming('identicalReaction duplicate modal resolved', { decision: userDecision });
            if (userDecision.startsWith('view:')) {
                const rxnId = userDecision.split(':')[1];
                window.location.href = `/?reaction_id=${rxnId}`;
                return false;
            }
        }
        loadingIndicator.style.display = 'flex';
        return true; // Proceed with reaction creation
    } catch (error) {
        console.error('Error checking for identical reaction:', error);
        return false;
    }
}

document.getElementById('reactionForm').addEventListener('submit', async function(e) {
    var loadingIndicator = document.getElementById('loadingIndicator');

    e.preventDefault(); // Prevent the default form submission

    var submitBtn = document.getElementById('submitBtn-form');
    const timingStart = performance.now();
    let timingLast = timingStart;
    const logTiming = function (label, details) {
        const now = performance.now();
        const delta = (now - timingLast).toFixed(1);
        const total = (now - timingStart).toFixed(1);
        console.log(`[TEMP ReactionTiming client submit] ${label}: +${delta}ms total=${total}ms`, details || '');
        timingLast = now;
    };
    const stopLoadingState = function () {
        delete submitBtn.dataset.submitting;
        if (window.ReactantsFormDirty && ReactantsFormDirty.isEditMode()) {
            ReactantsFormDirty.refreshButton();
        } else {
            submitBtn.disabled = false;
        }
        hideLoader();
    };

    const parsedUrl = new URL(window.location.href);
    const params = new URLSearchParams(parsedUrl.search);
    const reactionId = params.get('reaction_id');
    const action = params.get('action');
    const editing = action === 'edit';
    logTiming('submit:start', { action: action || 'create', reactionId: reactionId || null });

    if (editing) {
        if (window.ReactantsFormDirty && !ReactantsFormDirty.hasChanges()) {
            ReactantsFormDirty.refreshButton();
            logTiming('edit submit stopped: no Reactants changes');
            return;
        }
        logTiming('edit dirty check complete');
        // Context-aware confirm: saving a clone vs. updating an existing reaction.
        var inClone = !!(window.CloneReaction && CloneReaction.cloneId);
        var confirmOpts = inClone
            ? {
                title: 'Save clone',
                messageHtml: ReactantsFormDirty
                    ? ReactantsFormDirty.confirmationHtml().replace(
                        'This updates the saved reaction and cannot be undone.',
                        'Save these changes to this clone? The original reaction is not affected.'
                    )
                    : undefined,
                message: ReactantsFormDirty ? undefined : 'Save your changes as this new clone? The original reaction is not affected.',
                confirmText: 'Save clone',
                cancelText: 'Keep editing',
            }
            : {
                title: 'Update reaction',
                messageHtml: ReactantsFormDirty ? ReactantsFormDirty.confirmationHtml() : undefined,
                message: ReactantsFormDirty ? undefined : 'This updates the saved reaction and cannot be undone.',
                confirmText: 'Update reaction',
                cancelText: 'Cancel',
                danger: true,
            };
        logTiming('update confirm modal opening');
        var userConfirmed = await Notify.confirm(confirmOpts);
        logTiming('update confirm modal resolved', { confirmed: userConfirmed });
        if (!userConfirmed) {
            return; // Exit the function and do not submit form
        }
    }

    submitBtn.dataset.submitting = 'true';
    submitBtn.disabled = true;
    showLoader();
    logTiming('loader shown');

    var divElement = document.getElementById("organTags");
    var organTags = Array.from(divElement.getElementsByClassName('tag'))
                            .map(tag => tag.firstChild.textContent.trim());
    logTiming('organ tags read', { count: organTags.length });

    // Check if the subsystem field is filled
    var subsystemField = document.getElementById('subsystemField').value;
    if (!subsystemField.trim()) {
        Notify.error('Please enter a subsystem.');
        stopLoadingState();
        return; // Exit the function and do not submit form
    }

    var inputsGroups = getActiveReactionRows();
    if (inputsGroups.length === 0) {
        var noMetabolitesMessage = 'Enter at least one substrate or one product before creating a reaction.';
        showErrorModal(noMetabolitesMessage);
        window.scrollTo(0, 0);
        stopLoadingState();
        return;
    }
    
    for (var i = 0; i < inputsGroups.length; i++) {
        var group = inputsGroups[i];
        var statusDot = group.querySelector('.status-dot');
        if (statusDot.style.display === 'none') {
            var errorMessage = 'Verify all metabolites before creating reaction.';
            showErrorModal(errorMessage);
            window.scrollTo(0, 0);
            stopLoadingState();
            return
        }
    }
    logTiming('client validation complete', { activeRows: inputsGroups.length });
    // Continue with form submission

    const disabledInputs = this.querySelectorAll('input:disabled, select:disabled');
    disabledInputs.forEach(input => input.disabled = false);
    var formData = new FormData(this);

    disabledInputs.forEach(input => input.disabled = true);
    logTiming('FormData built', { disabledInputs: disabledInputs.length });
    const directionEl = document.getElementById('reactionDirection');
    formData.set('direction', normalizeReactionDirection(directionEl ? directionEl.value : 'forward'));
    var nameData = {};
    var metaboliteFields = document.querySelectorAll('.substrates-name, .products-name');
    var allNamesEntered = true;
    metaboliteFields.forEach(function(input, index) {
        var group = input.closest('.inputs-group');
        var key = input.name + (index + 1);
        var value = input.value;
        nameData[key] = value;

        if (reactionRowHasMetabolite(group) && input.value === '') {
            allNamesEntered = false;
        }
    });

    if (!allNamesEntered) {
        var errorMessage = 'Enter all metabolite names before creating reaction.';
        showErrorModal(errorMessage);
        window.scrollTo(0, 0);
        stopLoadingState();
        return; // Exit the function and do not submit form
    }
    loadingIndicator.style.display = 'flex';
    logTiming('metabolite names collected', { metaboliteNameFields: metaboliteFields.length });
    if (subsystemList.length === 0) {
        logTiming('updateSubsystems fetch:start');
        subsystemList = await updateSubsystems();
        logTiming('updateSubsystems fetch:complete', { count: subsystemList.length });
    }
    var isValidSubsystem = subsystemList.some(subsystem => subsystem.toLowerCase() === subsystemField.toLowerCase());
    logTiming('subsystem validation complete', { valid: isValidSubsystem });

    if (!isValidSubsystem) {
        logTiming('new subsystem confirm opening', { subsystem: subsystemField });
        var userConfirmed = await Notify.confirm({ title: 'New subsystem', message: `Are you sure you want to add a new subsystem "${subsystemField}"?`, confirmText: 'Add subsystem' });
        logTiming('new subsystem confirm resolved', { confirmed: userConfirmed });
        if (!userConfirmed) {
            var errorMessage = 'The subsystem entered is not valid.';
            showErrorModal(errorMessage);
            window.scrollTo(0, 0);
            stopLoadingState();
            return; // Exit the function and do not submit form
        } else {
            if (sessionStorage.getItem('userID') !== null) {
            // Add the new subsystem to the list
            subsystemList.push(subsystemField);
            persistSubsystemList(subsystemField); // Persist the new subsystem
            }
            else {
                var errorMessage = 'Please login to add a new subsystem.';
                showErrorModal(errorMessage);
                window.scrollTo(0, 0);
                stopLoadingState();
                return; // Exit the function and do not submit form
            }
        }
    }
    var skipAtomMapping = document.getElementById('skipAtomMapping').checked;
    formData.append('skipAtomMapping', skipAtomMapping);
    formData.append('nameData', JSON.stringify(nameData));
    formData.append('organs', JSON.stringify(organTags));
    formData.append('action', action);
    logTiming('form payload finalized', { skipAtomMapping: skipAtomMapping, action: action || 'create' });

    if (action === 'edit') {
        formData.append('action', action);

        formData.append('reaction_id', reactionId);
    }   
    formData.append('userID', sessionStorage.getItem('userID'));


    // Check if identical reaction exists
    if (!editing){
        const shouldProceed = await checkIdenticalReaction(formData, loadingIndicator, submitBtn, logTiming);
        if (!shouldProceed) {
            logTiming('create stopped after identicalReaction check');
            stopLoadingState();
            return;
        }
    }
    // Create a new reaction
    logTiming('inputReaction fetch:start');
    fetch(inputReactionUrl, {
        method: 'POST',
        body: formData,
        headers: {
            'X-Requested-With': 'XMLHttpRequest',
            'X-CSRFToken': csrfToken
        }
    })
    .then(response => {
        logTiming('inputReaction response headers', { status: response.status });
        return response.json();
    })
    .then(async (data) => {
        logTiming('inputReaction response json parsed', {
            status: data.status,
            reactionId: data.reaction_id || null,
        });
        if (data.status === 'error') {
            showErrorModal(data.message);
            window.scrollTo(0, 0);
            stopLoadingState();
            return;
        }
        else if (action !== 'edit' && data.reaction_id) { // reactionId is taken from the response
            const userID = sessionStorage.getItem('userID');
            const reactionId = data.reaction_id; // Ensure reactionId is obtained from the response
            
            // Send additional request to save CreatedReaction (to keep track of which user created which reaction)
            logTiming('create-reaction link fetch:start');
            fetch('create-reaction/', {
                method: 'POST',
                body: JSON.stringify({
                    user_id: userID,
                    reaction_id: reactionId
                }),
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                }
            })
            .then(response => {
                logTiming('create-reaction link response headers', { status: response.status });
                return response.json();
            })
            .then(result => {
                logTiming('create-reaction link response json parsed', { success: result.success });
                if (result.success) {
                    console.log('CreatedReaction saved successfully.');
                } else {
                    console.error('Failed to save CreatedReaction:', result.error);
                }
            })
            .catch(error => {
                console.error('Error saving CreatedReaction:', error);
            });
        }
        
        setTimeout(function() {
            const redirectUrl = window.location.origin + "/?reaction_id=" + data.reaction_id;
            logTiming('redirecting after successful submit', { url: action === 'edit' ? redirectUrl + "&action=edit" : redirectUrl });
            window.location.href = action === 'edit' ? redirectUrl + "&action=edit" : redirectUrl;
        }, 10); 
    })
    .catch(error => {
        logTiming('inputReaction fetch failed');
        console.error('Error:', error);
        const errorMessage = 'An unexpected error occurred.';
        showErrorModal(errorMessage);
        window.scrollTo(0, 0);
        stopLoadingState();
    });
});    


function persistSubsystemList(newSubsystem) {
    fetch('/update_subsystems/', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrfToken
        },
        body: JSON.stringify({ subsystems: [newSubsystem] })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            console.error('Error in updating subsystems:', data.message);
        } else {
            console.log('Subsystem list updated successfully.');
        }
    })
    .catch(error => console.error('Error:', error));
}


function showLoader() {
    document.getElementById('loadingIndicator').style.display = 'flex';
}

// Function to hide the loading indicator
function hideLoader() {
    document.getElementById('loadingIndicator').style.display = 'none';
}





function showErrorModal(message) {
    var errorMessageElement = document.getElementById('error-message');
    
    // Debugging output

    errorMessageElement.innerText = message;
    document.getElementById('error-message').style.display = 'block';
    document.getElementById('error-modal').style.display = 'block';
}




function updateStatusDots(containerId, foundList, miriamsList) {
    const container = document.getElementById(containerId);
    const statusDots = container.querySelectorAll('.status-dot');
    statusDots.forEach((dot, index) => {
        dot.style.display = 'block'; // Ensure the dot is displayed

        if (foundList[index]) {
            dot.className = 'status-dot found';
            dot.setAttribute('data-tooltip', 'Metabolite found in VMH');
            dot.style.backgroundColor = ''; // Reset color to default
        } else {
            dot.className = 'status-dot not-found';
            dot.setAttribute('data-tooltip', 'Metabolite not found in VMH');
            dot.style.backgroundColor = ''; // Reset color to default
        }

        if (foundList[index] && miriamsList[index]) {
            dot.onclick = () => window.open(miriamsList[index], '_blank');
            dot.style.cursor = 'pointer';
        } else {
            dot.onclick = null;
            dot.style.cursor = 'default';
        }
    });
}

function updateStatusDot(dot, found, miriam) {
    dot.style.display = 'block'; // Ensure the dot is displayed

    if (found) {
        dot.className = 'status-dot found';
        dot.setAttribute('data-tooltip', 'Metabolite found in VMH');
        dot.style.backgroundColor = ''; // Reset color to default
    } else {
        dot.className = 'status-dot not-found';
        dot.setAttribute('data-tooltip', 'Metabolite not found in VMH');
        dot.style.backgroundColor = ''; // Reset color to default
    }

    if (found && miriam) {
        dot.onclick = () => window.open(miriam, '_blank');
        dot.style.cursor = 'pointer';
    } else {
        dot.onclick = null;
        dot.style.cursor = 'default';
    }
}

function confirmAll() {
    // Get all the elements with the class 'done-field-btn-all'
    const verifyAllButtons = document.querySelectorAll('.done-field-btn-all');

    // Loop through each button and trigger a click event
    verifyAllButtons.forEach(button => button.click());
}


function setupTooltips() {
    const OFFSET = -20;

    const positionTooltip = (trigger, tooltip) => {
        if (!tooltip) return;
        const rect = trigger.getBoundingClientRect();
        const tooltipRect = tooltip.getBoundingClientRect();

        let top = rect.top - tooltipRect.height - OFFSET;
        let left = rect.right + OFFSET;

        const maxLeft = Math.max(0, window.innerWidth - tooltipRect.width - 12);
        const minTop = 8;

        if (left > maxLeft) {
            left = maxLeft;
        }
        if (top < minTop) {
            top = minTop;
        }

        tooltip.style.top = `${top}px`;
        tooltip.style.left = `${left}px`;
    };

    document.querySelectorAll('.info-symbol').forEach((item) => {
        let activeTooltip = null;

        const removeTooltip = () => {
            if (activeTooltip && activeTooltip.parentNode) {
                activeTooltip.parentNode.removeChild(activeTooltip);
            }
            activeTooltip = null;
        };

        item.addEventListener('mouseenter', function () {
            const tooltipContent = this.getAttribute('data-tooltip-content');
            if (!tooltipContent) return;

            const tooltip = document.createElement('div');
            tooltip.className = 'info-tooltip';
            tooltip.innerHTML = tooltipContent.replace(/\n/g, '<br />');
            document.body.appendChild(tooltip);
            activeTooltip = tooltip;
            positionTooltip(this, tooltip);
        });

        item.addEventListener('mousemove', function () {
            if (activeTooltip) {
                positionTooltip(this, activeTooltip);
            }
        });

        item.addEventListener('mouseleave', () => {
            removeTooltip();
        });

        item.addEventListener('blur', removeTooltip);
    });
}
