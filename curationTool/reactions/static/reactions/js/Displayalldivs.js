/**
 * Displayalldivs.js — orchestrates rendering a reaction into the detail panels.
 *
 * Multi-reaction refactor (001-multi-reaction-tabbed, FR-003 / SC-006):
 * `displayDivs` now funnels the reaction object through `OpenReactionState` so the
 * renderers draw from the active open-reaction entry rather than an ad-hoc global.
 * Each open reaction therefore has isolated state keyed by its id; switching the
 * selection (see `loadReactionById`) swaps which entry the panels show.
 */
function displayDivs(reactionData) {
    // Record this reaction as the active open-reaction state (SC-006). The id is
    // the DB reaction_id when known (URL param), else a temp id (e.g. a clone).
    if (window.OpenReactionState) {
        var urlParams = new URLSearchParams(window.location.search);
        var id = urlParams.get('reaction_id') || OpenReactionState.getActiveId();
        if (id === null || id === undefined) {
            id = OpenReactionState.newTempId();
        }
        OpenReactionState.setReaction(id, reactionData, { active: true });
        reactionData = OpenReactionState.getActiveReaction();
    }
    loadAtomMappingDiv(reactionData);
    loadChemInfoDiv(reactionData);
    loadMetaboliteInfoDiv(reactionData);
    refreshSideButtons();
    updateStatusDots('substratesDiv', reactionData.subs_found, reactionData.subs_miriams);
    updateStatusDots('productsDiv', reactionData.prod_found, reactionData.prod_miriams);
    displayReactionMessage(reactionData);
    displayreactioninfo(reactionData);
}

/**
 * Load a saved reaction into the detail panels in-page (no full reload) and make
 * it the active open reaction. Reuses the proven load path used on page load
 * (updateFormFields → confirmAll → displayDivs) and rewrites the URL's
 * reaction_id via history.replaceState so the existing save / add-info handlers —
 * which read reaction_id from the URL — target the newly selected reaction.
 *
 * Used by groupView.js when the curator selects a reaction card (FR-002).
 * Returns a Promise that resolves when the panels have been populated.
 */
function loadReactionById(reactionId) {
    if (reactionId === null || reactionId === undefined) {
        return Promise.reject(new Error('No reaction id'));
    }
    // Clear any pending session gene-info bucket BEFORE loading so it cannot be
    // flushed into the reaction we are about to open (FR-003 isolation). Without
    // this, a gene/GPR staged while creating a reaction would leak into every
    // reaction subsequently selected.
    var clearGuard = (typeof clearGeneSessionUrl !== 'undefined')
        ? fetch(clearGeneSessionUrl, {
            method: 'POST',
            headers: { 'X-Requested-With': 'XMLHttpRequest', 'X-CSRFToken': csrfToken },
        }).catch(function () {})
        : Promise.resolve();

    return clearGuard
        .then(function () { return fetch(getReaction + reactionId); })
        .then(function (response) {
            if (!response.ok) {
                throw new Error('HTTP error! status: ' + response.status);
            }
            return response.json();
        })
        .then(async function (reactionData) {
            await updateFormFields(reactionData);
            confirmAll();
            // Selecting a group member opens it for editing: put the page into
            // edit mode (?reaction_id=<id>&action=edit) so save/add-info target
            // this reaction and the primary button reads "Update Reaction".
            history.replaceState({}, '',
                window.location.pathname + '?reaction_id=' + reactionId + '&action=edit');
            var submitBtn = document.getElementById('submitBtn-form');
            if (submitBtn && submitBtn.childNodes.length > 2) {
                submitBtn.childNodes[2].nodeValue = 'Update Reaction';
            }
            updateEditModeUI();
            displayDivs(reactionData);
            if (reactionData.short_name) {
                var reactionDataName = reactionData.short_name;
                if (reactionData.short_name.length > 20) {
                    reactionDataName = reactionData.short_name.substring(0, 30) + '...';
                }
                var reactionStatusInfo = {
                    name: reactionDataName,
                    description: reactionData.description,
                };
                setLoggedInStatusBasedOnUrl(reactionStatusInfo);
                DisplayTag(reactionData.Organs);
            }
            if (window.OpenReactionState) {
                OpenReactionState.setDirty(reactionId, false);
            }
            return reactionData;
        });
}

window.loadReactionById = loadReactionById;

/**
 * Show/hide action buttons that only make sense outside edit mode. In edit mode
 * (?action=edit) the reaction is an existing saved row, so "Save Reaction" —
 * which only adds a brand-new reaction to the saved list — is redundant and is
 * hidden. Called on initial load, on card selection, and when leaving clone mode.
 */
function updateEditModeUI() {
    var action = new URLSearchParams(window.location.search).get('action');
    var saveBtn = document.getElementById('saveReactionButton');
    if (saveBtn) {
        saveBtn.style.display = (action === 'edit') ? 'none' : '';
    }
}

window.updateEditModeUI = updateEditModeUI;
