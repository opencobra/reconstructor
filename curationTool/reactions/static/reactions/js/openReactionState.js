/**
 * openReactionState.js — per-open-reaction state store
 * (feature 001-multi-reaction-tabbed, FR-003 / SC-006).
 *
 * Today a reaction is loaded by Displayalldivs.js and passed as a `reactionData`
 * object into the display*.js renderers, while the live *editing* state lives in
 * the DOM panels. To support several open reactions at once without one bleeding
 * into another, this module keeps a store keyed by an `openReactionId`:
 *
 *   - For a saved group member, the openReactionId is the Reaction's DB id.
 *   - For a clone (US4), it is a temporary client id (e.g. "tmp-<n>") until saved.
 *
 * Each entry mirrors the shape the current `reactionData` carries plus a `dirty`
 * flag (unsaved-edit guard, FR-017) and a `temp` marker (unsaved clone, FR-015).
 *
 * This is an additive module: it introduces no behaviour change on its own. The
 * display*.js renderers are migrated to read `getActiveReaction()` in a separate,
 * incremental step (T007) so the single-reaction page keeps working throughout.
 */
(function (global) {
    'use strict';

    const store = new Map();   // openReactionId -> { data, dirty, temp }
    let activeId = null;
    let tempCounter = 0;

    function newTempId() {
        tempCounter += 1;
        return 'tmp-' + tempCounter;
    }

    /**
     * Insert or replace the state for an open reaction.
     * @param {string|number} id      openReactionId (DB id or temp id)
     * @param {object} data           the reactionData-shaped object
     * @param {object} [opts]         { temp: bool, active: bool }
     * @returns {string|number} the id used
     */
    function setReaction(id, data, opts) {
        opts = opts || {};
        if (id === undefined || id === null) {
            id = newTempId();
        }
        store.set(id, {
            data: data || {},
            dirty: false,
            temp: !!opts.temp,
        });
        if (opts.active || activeId === null) {
            activeId = id;
        }
        return id;
    }

    /** Mark which open reaction the detail panels currently show. */
    function setActiveReaction(id) {
        if (store.has(id)) {
            activeId = id;
            return true;
        }
        return false;
    }

    function getActiveId() {
        return activeId;
    }

    /** @returns {object|null} the active reaction's reactionData-shaped object. */
    function getActiveReaction() {
        const entry = activeId === null ? null : store.get(activeId);
        return entry ? entry.data : null;
    }

    function getReaction(id) {
        const entry = store.get(id);
        return entry ? entry.data : null;
    }

    function hasReaction(id) {
        return store.has(id);
    }

    /** Replace just the data for an existing entry (e.g. after a re-fetch). */
    function updateReaction(id, data) {
        const entry = store.get(id);
        if (entry) {
            entry.data = data || {};
        }
    }

    function setDirty(id, value) {
        const entry = store.get(id);
        if (entry) {
            entry.dirty = value === undefined ? true : !!value;
        }
    }

    function isDirty(id) {
        const entry = store.get(id === undefined ? activeId : id);
        return !!(entry && entry.dirty);
    }

    function isTemp(id) {
        const entry = store.get(id === undefined ? activeId : id);
        return !!(entry && entry.temp);
    }

    /** Remove an open reaction (e.g. discarding an unsaved clone, FR-015). */
    function removeReaction(id) {
        const existed = store.delete(id);
        if (activeId === id) {
            activeId = store.size ? store.keys().next().value : null;
        }
        return existed;
    }

    function clear() {
        store.clear();
        activeId = null;
    }

    function ids() {
        return Array.from(store.keys());
    }

    /**
     * Create a deep-copied temp entry from an existing reaction's state, used by
     * the clone-with-edits flow (US4 / FR-013). Returns the new temp id.
     */
    function cloneReactionState(sourceId) {
        const source = getReaction(sourceId);
        if (!source) {
            return null;
        }
        const copy = JSON.parse(JSON.stringify(source));
        const id = newTempId();
        store.set(id, { data: copy, dirty: false, temp: true });
        return id;
    }

    global.OpenReactionState = {
        setReaction,
        setActiveReaction,
        getActiveReaction,
        getActiveId,
        getReaction,
        hasReaction,
        updateReaction,
        setDirty,
        isDirty,
        isTemp,
        removeReaction,
        cloneReactionState,
        clear,
        ids,
        newTempId,
    };
})(window);
