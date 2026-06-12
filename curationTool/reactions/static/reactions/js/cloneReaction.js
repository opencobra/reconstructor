/**
 * cloneReaction.js — clone-with-edits flow (001-multi-reaction-tabbed, US4;
 * FR-013, FR-014, FR-015).
 *
 * Cloning is done SERVER-SIDE: the active reaction is copied into a brand-new
 * Reaction row (its own reaction_id) which is added to the current group and
 * opened in edit mode. Every edit, GPR change, reference and save then binds to
 * the clone's id — the original row is never referenced, so it can never change
 * (FR-014). Discarding deletes the clone row outright, persisting nothing
 * (FR-015).
 *
 * This replaces an earlier client-side "create mode" approach that reused the
 * original's data object (which carries reaction_id) and so leaked edits back
 * into the original — that is the exact bug this rewrite eliminates.
 */
(function (global) {
    'use strict';

    function userId() { return sessionStorage.getItem('userID'); }
    function csrf() { return (typeof csrfToken !== 'undefined') ? csrfToken : ''; }

    function postJSON(url, body) {
        return fetch(url, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrf(),
            },
            body: JSON.stringify(body),
        }).then(function (r) { return r.json(); });
    }

    function currentReactionId() {
        return new URLSearchParams(window.location.search).get('reaction_id');
    }

    var CloneReaction = {
        cloneId: null,

        /** Clone the currently open saved reaction into a new editable variant. */
        clone: function () {
            if (!userId()) {
                Notify.error('Log in to clone reactions.');
                return;
            }
            var reactionId = currentReactionId();
            if (!reactionId) {
                Notify.error('Open a saved reaction first, then clone it.');
                return;
            }
            var data = window.OpenReactionState ? OpenReactionState.getActiveReaction() : null;
            var baseName = (data && data.short_name) ? data.short_name : 'reaction';
            var self = this;

            var proceed = global.GroupView ? GroupView.guardUnsaved() : Promise.resolve(true);
            proceed.then(function (ok) {
                if (!ok) return;
                if (typeof showLoader === 'function') showLoader();

                // Server-side deep copy → new row with its own id.
                var form = new FormData();
                form.append('reaction_id', reactionId);
                form.append('userID', userId());
                form.append('name', baseName + ' (clone)');
                fetch(global.cloneReactionUrl, {
                    method: 'POST',
                    headers: { 'X-Requested-With': 'XMLHttpRequest', 'X-CSRFToken': csrf() },
                    body: form,
                })
                    .then(function (r) { return r.json(); })
                    .then(function (res) {
                        if (res.status !== 'success' || !res.new_reaction_id) {
                            if (typeof hideLoader === 'function') hideLoader();
                            Notify.error('Could not clone: ' + (res.message || 'unknown error'));
                            return;
                        }
                        self.cloneId = res.new_reaction_id;
                        var groupId = global.GroupView ? GroupView.activeGroupId : null;
                        var addStep = groupId
                            ? postJSON(global.groupAddUrl, {
                                userID: userId(), groupId: groupId,
                                reactionIds: [self.cloneId],
                            })
                            : Promise.resolve();
                        addStep.then(function () { self.openClone(self.cloneId, baseName); });
                    })
                    .catch(function (err) {
                        if (typeof hideLoader === 'function') hideLoader();
                        console.error('Clone failed', err);
                        Notify.error('Could not clone the reaction.');
                    });
            });
        },

        /** Open the freshly created clone in edit mode and switch into clone UI. */
        openClone: function (cloneId, sourceName) {
            var self = this;
            var loaded = (typeof window.loadReactionById === 'function')
                ? window.loadReactionById(cloneId)
                : Promise.reject();
            loaded.then(function () {
                if (typeof hideLoader === 'function') hideLoader();
                self.enterCloneUI(sourceName);
                if (global.GroupView) GroupView.loadGroup(GroupView.activeGroupId);
                if (global.TabsController) TabsController.refreshCounts();
            }).catch(function () {
                // Fallback: navigate to the clone in edit mode.
                window.location.href = '/?reaction_id=' + cloneId + '&action=edit';
            });
        },

        /** Switch the action bar into a focused, unambiguous clone-editing mode. */
        enterCloneUI: function (sourceName) {
            this.toggleButtons(true);
            var banner = document.getElementById('cloneModeBanner');
            if (banner) {
                var label = banner.querySelector('.clone-banner-label');
                if (label) {
                    label.textContent = sourceName
                        ? 'Editing a new clone of "' + sourceName + '". Make your changes, then Save clone. The original is untouched.'
                        : 'Editing a new clone. Make your changes, then Save clone. The original is untouched.';
                }
                banner.style.display = 'flex';
            }
        },

        /**
         * Hide the standalone Clone button while already editing a clone (avoids
         * clone-of-clone ambiguity). The "Save Reaction" button is owned by the
         * edit-mode logic (updateEditModeUI) since a clone always opens in edit
         * mode, so it is hidden there automatically.
         */
        toggleButtons: function (inClone) {
            var clone = document.getElementById('cloneReactionButton');
            if (clone) clone.style.display = inClone ? 'none' : '';
        },

        /** Persist the clone's edits via the normal Update path, then it becomes
         *  an ordinary saved reaction (banner clears on the resulting reload). */
        saveClone: function () {
            var submitBtn = document.getElementById('submitBtn-form');
            if (submitBtn) submitBtn.click();
        },

        /** Discard the clone — delete the row so nothing is persisted (FR-015). */
        discard: function () {
            var id = this.cloneId || currentReactionId();
            if (!id) { this.exitCloneUI(); return; }
            var self = this;
            Notify.confirm({
                title: 'Discard clone',
                message: 'Discard this clone? It will be deleted and nothing will be saved.',
                confirmText: 'Discard',
                cancelText: 'Keep editing',
                danger: true,
            }).then(function (ok) {
                if (!ok) return;
                if (window.OpenReactionState) OpenReactionState.setDirty(undefined, false);
                postJSON(global.groupDiscardCloneUrl, { userID: userId(), reactionId: id })
                    .then(function () {
                        self.cloneId = null;
                        // Return to a clean workspace; the clone no longer exists.
                        window.location.href = window.location.origin;
                    })
                    .catch(function (err) {
                        console.error('Discard failed', err);
                        Notify.error('Could not discard the clone.');
                    });
            });
        },

        exitCloneUI: function () {
            // Leaving the clone-editing context: forget the clone id so a later
            // normal update is not mistaken for a clone save (context-aware confirm).
            this.cloneId = null;
            this.toggleButtons(false);
            var banner = document.getElementById('cloneModeBanner');
            if (banner) banner.style.display = 'none';
            // Re-apply edit-mode button visibility for the reaction now open.
            if (typeof window.updateEditModeUI === 'function') updateEditModeUI();
        },
    };

    global.CloneReaction = CloneReaction;

    document.addEventListener('DOMContentLoaded', function () {
        var cloneBtn = document.getElementById('cloneReactionButton');
        if (cloneBtn) cloneBtn.addEventListener('click', function () { CloneReaction.clone(); });
        var saveBtn = document.getElementById('saveCloneButton');
        if (saveBtn) saveBtn.addEventListener('click', function () { CloneReaction.saveClone(); });
        var discardBtn = document.getElementById('discardCloneButton');
        if (discardBtn) discardBtn.addEventListener('click', function () { CloneReaction.discard(); });
    });
})(window);
