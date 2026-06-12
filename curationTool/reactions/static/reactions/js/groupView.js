/**
 * groupView.js — renders the active group's reactions as a rail of cards and
 * drives selecting / adding / removing members (001-multi-reaction-tabbed,
 * US1 / US3; FR-001, FR-002, FR-010, FR-011, FR-012).
 *
 * Selecting a card loads that reaction into the existing detail panels in-page
 * via window.loadReactionById (no full reload). Unsaved edits to the current
 * reaction are guarded before switching (FR-017). Removing a card drops the
 * reaction from the group only — the saved reaction is untouched (FR-011).
 */
(function (global) {
    'use strict';

    var RAIL_ID = 'groupRail';
    var EMPTY_ID = 'groupEmptyState';

    function userId() {
        return sessionStorage.getItem('userID');
    }

    function csrf() {
        return (typeof csrfToken !== 'undefined') ? csrfToken : '';
    }

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
        var params = new URLSearchParams(window.location.search);
        return params.get('reaction_id');
    }

    var GroupView = {
        activeGroupId: null,
        pickerReactions: [],
        pendingSelection: false,
        loadingReactionId: null,

        rail: function () { return document.getElementById(RAIL_ID); },

        /** Fetch a group's contents and render them (FR-001). */
        loadGroup: function (groupId) {
            var self = this;
            if (!userId()) {
                this.renderEmpty(true);
                return Promise.resolve(null);
            }
            var url = global.groupContentsUrl + '?userID=' + encodeURIComponent(userId());
            if (groupId) {
                url += '&groupId=' + encodeURIComponent(groupId);
            }
            return fetch(url)
                .then(function (r) { return r.json(); })
                .then(function (data) {
                    if (data.status !== 'success') {
                        self.renderEmpty(true);
                        return data;
                    }
                    self.activeGroupId = data.group.id;
                    self.render(data.reactions, data.group);
                    return data;
                })
                .catch(function (err) {
                    console.error('Failed to load group', err);
                });
        },

        /** Render reaction cards into the rail. */
        render: function (reactions, group) {
            var rail = this.rail();
            if (!rail) return;
            rail.innerHTML = '';
            if (!reactions || reactions.length === 0) {
                this.renderEmpty(false);
                return;
            }
            this.hideEmpty();
            var activeId = currentReactionId();
            var self = this;
            reactions.forEach(function (rxn) {
                rail.appendChild(self.buildCard(rxn, String(rxn.id) === String(activeId)));
            });
        },

        buildCard: function (rxn, isActive) {
            var card = document.createElement('div');
            card.className = 'reaction-card' + (isActive ? ' active' : '');
            card.setAttribute('data-reaction-id', rxn.id);
            card.setAttribute('role', 'button');
            card.setAttribute('tabindex', '0');
            card.setAttribute('aria-pressed', isActive ? 'true' : 'false');

            var massOk = rxn.mass_balanced;
            var chargeOk = rxn.charge_balanced;
            var balanceClass = (massOk && chargeOk) ? 'balanced' : 'unbalanced';
            var balanceLabel = (massOk && chargeOk) ? 'Balanced'
                : ((massOk ? '' : 'Mass ') + (chargeOk ? '' : 'Charge ')).trim() + ' unbalanced';

            var title = rxn.short_name || ('Reaction ' + rxn.id);
            var formula = rxn.rxn_formula || rxn.molc_formula || '';

            card.innerHTML =
                '<button type="button" class="reaction-card-remove" title="Remove from group" aria-label="Remove from group">' +
                '<i class="fas fa-times" aria-hidden="true"></i></button>' +
                '<div class="reaction-card-title"></div>' +
                '<div class="reaction-card-formula"></div>' +
                '<div class="reaction-card-meta">' +
                '<span class="reaction-card-balance ' + balanceClass + '">' + balanceLabel + '</span>' +
                (rxn.subsystem ? '<span class="reaction-card-subsystem"></span>' : '') +
                '</div>';
            var titleEl = card.querySelector('.reaction-card-title');
            titleEl.textContent = title;
            titleEl.setAttribute('title', title);
            var formulaEl = card.querySelector('.reaction-card-formula');
            formulaEl.textContent = formula;
            // Full formula on hover; the cell itself truncates with an ellipsis.
            if (formula) formulaEl.setAttribute('title', formula);
            if (rxn.subsystem) {
                var subEl = card.querySelector('.reaction-card-subsystem');
                subEl.textContent = rxn.subsystem;
                subEl.setAttribute('title', rxn.subsystem);
            }

            var self = this;
            var select = function () { self.selectCard(rxn.id, card); };
            card.addEventListener('click', function (e) {
                if (e.target.closest('.reaction-card-remove')) return;
                select();
            });
            card.addEventListener('keydown', function (e) {
                if (e.key === 'Enter' || e.key === ' ') {
                    e.preventDefault();
                    select();
                }
            });
            card.querySelector('.reaction-card-remove').addEventListener('click', function (e) {
                e.stopPropagation();
                self.removeCard(rxn.id, rxn.short_name);
            });
            return card;
        },

        /** Select a reaction card → load it into the detail panels (FR-002). */
        selectCard: function (reactionId, card) {
            if (this.pendingSelection || this.loadingReactionId !== null) {
                return;
            }
            this.pendingSelection = true;
            var self = this;
            this.guardUnsaved().then(function (proceed) {
                if (!proceed) {
                    self.pendingSelection = false;
                    return;
                }
                self.setOpeningCard(card, reactionId);
                if (typeof window.loadReactionById !== 'function') {
                    // Fallback: navigate (still reuses the proven load path).
                    window.location.href = '/?reaction_id=' + reactionId;
                    return;
                }
                window.loadReactionById(reactionId)
                    .then(function () {
                        self.markActive(card);
                        if (window.WorkspacePanels) window.WorkspacePanels.refresh();
                        // Leaving a clone for another reaction ends clone mode;
                        // the clone persists as an ordinary saved reaction.
                        if (global.CloneReaction) CloneReaction.exitCloneUI();
                    })
                    .catch(function (err) {
                        console.error('Failed to open reaction; navigating instead', err);
                        window.location.href = '/?reaction_id=' + reactionId;
                    })
                    .finally(function () {
                        self.clearOpeningCard();
                    });
            }).catch(function (err) {
                console.error('Could not switch reactions', err);
                self.clearOpeningCard();
            });
        },

        setOpeningCard: function (card, reactionId) {
            this.loadingReactionId = reactionId;
            var rail = this.rail();
            if (rail) {
                rail.classList.add('is-loading-reaction');
                rail.querySelectorAll('.reaction-card.is-opening').forEach(function (c) {
                    c.classList.remove('is-opening');
                    c.removeAttribute('aria-busy');
                });
            }
            if (card) {
                card.classList.add('is-opening');
                card.setAttribute('aria-busy', 'true');
            }
        },

        clearOpeningCard: function () {
            this.pendingSelection = false;
            this.loadingReactionId = null;
            var rail = this.rail();
            if (!rail) return;
            rail.classList.remove('is-loading-reaction');
            rail.querySelectorAll('.reaction-card.is-opening').forEach(function (card) {
                card.classList.remove('is-opening');
                card.removeAttribute('aria-busy');
            });
        },

        markActive: function (card) {
            var rail = this.rail();
            if (rail) {
                rail.querySelectorAll('.reaction-card.active').forEach(function (c) {
                    c.classList.remove('active');
                    c.setAttribute('aria-pressed', 'false');
                });
            }
            if (card) {
                card.classList.add('active');
                card.setAttribute('aria-pressed', 'true');
            }
        },

        /**
         * Unsaved-edit guard (FR-017): if the active reaction has unsaved edits,
         * confirm before switching away. Returns Promise<boolean> (proceed?).
         */
        guardUnsaved: function () {
            var dirty = window.OpenReactionState && OpenReactionState.isDirty();
            if (!dirty) return Promise.resolve(true);
            return Notify.confirm({
                title: 'Unsaved changes',
                message: 'You have unsaved changes in the open reaction. Switch anyway and discard them?',
                confirmText: 'Discard & switch',
                cancelText: 'Stay',
                danger: true,
            });
        },

        /** Remove a reaction from the current group (FR-011, no data loss). */
        removeCard: function (reactionId, name) {
            if (!userId()) return;
            var self = this;
            Notify.confirm({
                title: 'Remove from group',
                message: 'Remove "' + (name || 'this reaction') + '" from the group?\n' +
                    'It will stay in your Saved Reactions.',
                confirmText: 'Remove',
                cancelText: 'Cancel',
            }).then(function (ok) {
                if (!ok) return;
                postJSON(global.groupRemoveUrl, {
                    userID: userId(),
                    groupId: self.activeGroupId,
                    reactionId: reactionId,
                }).then(function () {
                    self.loadGroup(self.activeGroupId);
                    if (global.TabsController) global.TabsController.refreshCounts();
                });
            });
        },

        /** Add saved reactions to the current group (FR-010, FR-012). */
        addReactions: function (reactionIds) {
            if (!userId() || !reactionIds || !reactionIds.length) return Promise.resolve();
            var self = this;
            return postJSON(global.groupAddUrl, {
                userID: userId(),
                groupId: this.activeGroupId,
                reactionIds: reactionIds,
            }).then(function (res) {
                self.loadGroup(self.activeGroupId);
                if (global.TabsController) global.TabsController.refreshCounts();
                return res;
            });
        },

        renderEmpty: function (anonymous) {
            var rail = this.rail();
            if (rail) rail.innerHTML = '';
            var empty = document.getElementById(EMPTY_ID);
            if (empty) {
                empty.style.display = 'flex';
                var msg = empty.querySelector('.group-empty-message');
                if (msg) {
                    msg.textContent = anonymous
                        ? 'Log in to organise reactions into groups. You can still create a reaction below.'
                        : 'This group is empty. Create a reaction or add one from your saved reactions.';
                }
            }
        },

        hideEmpty: function () {
            var empty = document.getElementById(EMPTY_ID);
            if (empty) empty.style.display = 'none';
        },

        // ---- Add-from-saved picker (FR-010, FR-012) ----------------------
        openPicker: function () {
            if (!userId()) {
                Notify.error('Log in to add reactions to a group.');
                return;
            }
            var list = document.getElementById('addFromSavedList');
            if (list) list.innerHTML = '<div class="picker-loading">Loading saved reactions…</div>';
            this.resetPickerFilters();
            var self = this;
            var url = global.groupSavedUrl + '?userID=' + encodeURIComponent(userId());
            if (this.activeGroupId) url += '&groupId=' + encodeURIComponent(this.activeGroupId);
            fetch(url)
                .then(function (r) { return r.json(); })
                .then(function (data) {
                    self.renderPicker(data.reactions || []);
                    if (global.jQuery) jQuery('#addFromSavedModal').modal('show');
                });
        },

        renderPicker: function (reactions) {
            this.pickerReactions = reactions || [];
            this.populatePickerFilters(this.pickerReactions);
            var list = document.getElementById('addFromSavedList');
            if (!list) return;
            list.innerHTML = '';
            if (!this.pickerReactions.length) {
                list.innerHTML = '<div class="picker-empty">You have no saved reactions yet.</div>';
                this.updatePickerSummary(0, 0);
                return;
            }
            this.pickerReactions.forEach(function (rxn) {
                var row = document.createElement('label');
                row.className = 'picker-row' + (rxn.in_group ? ' in-group' : '');
                row.setAttribute('data-reaction-id', rxn.id);
                row.setAttribute('data-in-group', rxn.in_group ? 'true' : 'false');
                row.setAttribute('data-search-all', GroupView.pickerSearchText(rxn, 'all'));
                row.setAttribute('data-search-name', GroupView.pickerSearchText(rxn, 'name'));
                row.setAttribute('data-search-subsystem', GroupView.pickerSearchText(rxn, 'subsystem'));
                row.setAttribute('data-search-flags', GroupView.pickerSearchText(rxn, 'flags'));
                row.setAttribute('data-search-substrates', GroupView.pickerSearchText(rxn, 'substrates'));
                row.setAttribute('data-search-products', GroupView.pickerSearchText(rxn, 'products'));
                row.setAttribute('data-search-formula', GroupView.pickerSearchText(rxn, 'formula'));
                row.setAttribute('data-search-direction', GroupView.pickerSearchText(rxn, 'direction'));
                row.setAttribute('data-subsystem', (rxn.subsystem || '').toLowerCase());
                row.setAttribute('data-flag-ids', (rxn.flags || []).map(function (flag) {
                    return String(flag.id);
                }).join(','));
                var disabled = rxn.in_group ? 'disabled checked' : '';
                row.innerHTML =
                    '<input type="checkbox" class="picker-check" value="' + rxn.id + '" ' + disabled + '>' +
                    '<span class="picker-main">' +
                    '<span class="picker-name"></span>' +
                    '<span class="picker-formula"></span>' +
                    '<span class="picker-meta"></span>' +
                    '</span>' +
                    (rxn.in_group ? '<span class="picker-tag">Already in group</span>' : '');
                row.querySelector('.picker-name').textContent = rxn.short_name || ('Reaction ' + rxn.id);
                row.querySelector('.picker-formula').textContent = rxn.rxn_formula || rxn.molc_formula || '';
                row.querySelector('.picker-formula').setAttribute('title', rxn.rxn_formula || rxn.molc_formula || '');
                GroupView.renderPickerMeta(row.querySelector('.picker-meta'), rxn);
                list.appendChild(row);
            });
            this.applyPickerFilters();
        },

        renderPickerMeta: function (meta, rxn) {
            if (!meta) return;
            meta.innerHTML = '';
            if (rxn.subsystem) {
                var subsystem = document.createElement('span');
                subsystem.className = 'picker-meta-chip';
                subsystem.textContent = rxn.subsystem;
                subsystem.setAttribute('title', rxn.subsystem);
                meta.appendChild(subsystem);
            }
            (rxn.flags || []).forEach(function (flag) {
                var flagChip = document.createElement('span');
                flagChip.className = 'picker-flag-chip';
                flagChip.setAttribute('title', flag.name || 'Flag');
                if (flag.color) flagChip.style.borderColor = flag.color;
                var icon = document.createElement('i');
                icon.className = 'fas fa-flag';
                if (flag.color) icon.style.color = flag.color;
                flagChip.appendChild(icon);
                flagChip.appendChild(document.createTextNode(flag.name || 'Flag'));
                meta.appendChild(flagChip);
            });
        },

        pickerSearchText: function (rxn, field) {
            var flags = (rxn.flags || []).map(function (flag) {
                return flag.name || '';
            }).join(' ');
            var formula = [rxn.rxn_formula || '', rxn.molc_formula || ''].join(' ');
            var map = {
                name: rxn.short_name || '',
                subsystem: rxn.subsystem || '',
                flags: flags,
                substrates: rxn.substrates || '',
                products: rxn.products || '',
                formula: formula,
                direction: rxn.direction || '',
            };
            if (field !== 'all') return (map[field] || '').toLowerCase();
            return [
                map.name, map.subsystem, map.flags, map.substrates,
                map.products, map.formula, map.direction,
            ].join(' ').toLowerCase();
        },

        populatePickerFilters: function (reactions) {
            var subsystemSelect = document.getElementById('addFromSavedSubsystemFilter');
            var flagMenu = document.getElementById('addFromSavedFlagMenu');
            if (subsystemSelect) {
                subsystemSelect.innerHTML = '<option value="">All subsystems</option>';
                this.uniqueSorted(reactions.map(function (rxn) {
                    return rxn.subsystem || '';
                })).forEach(function (subsystem) {
                    var option = document.createElement('option');
                    option.value = subsystem.toLowerCase();
                    option.textContent = subsystem;
                    subsystemSelect.appendChild(option);
                });
            }
            if (flagMenu) {
                flagMenu.innerHTML = '';
                var seen = {};
                reactions.forEach(function (rxn) {
                    (rxn.flags || []).forEach(function (flag) {
                        var key = String(flag.id);
                        if (!key || seen[key]) return;
                        seen[key] = flag;
                    });
                });
                flagMenu.appendChild(this.buildFlagFilterOption({
                    id: '',
                    name: 'All flags',
                    color: '#94a3b8',
                }));
                Object.keys(seen)
                    .sort(function (a, b) {
                        return (seen[a].name || '').localeCompare(seen[b].name || '');
                    })
                    .forEach(function (key) {
                        flagMenu.appendChild(GroupView.buildFlagFilterOption(seen[key]));
                    });
                var label = document.getElementById('addFromSavedFlagLabel');
                var icon = document.getElementById('addFromSavedFlagIcon');
                this.setPickerFlagFilter({
                    id: this.selectedPickerFlag(),
                    name: (label && label.textContent) || 'All flags',
                    color: (icon && icon.style.color) || '#94a3b8',
                });
            }
        },

        buildFlagFilterOption: function (flag) {
            var item = document.createElement('button');
            item.type = 'button';
            item.className = 'add-from-saved-flag-option';
            item.setAttribute('role', 'option');
            item.setAttribute('data-flag-id', flag.id || '');
            item.setAttribute('data-flag-name', flag.name || 'Flag');
            item.setAttribute('data-flag-color', flag.color || '#94a3b8');
            item.innerHTML =
                '<i class="fas fa-flag" aria-hidden="true"></i>' +
                '<span></span>';
            item.querySelector('i').style.color = flag.color || '#94a3b8';
            item.querySelector('span').textContent = flag.name || 'Flag';
            item.addEventListener('click', function () {
                GroupView.setPickerFlagFilter({
                    id: item.getAttribute('data-flag-id') || '',
                    name: item.getAttribute('data-flag-name') || 'All flags',
                    color: item.getAttribute('data-flag-color') || '#94a3b8',
                });
                GroupView.closePickerFlagMenu();
                GroupView.applyPickerFilters();
            });
            return item;
        },

        setPickerFlagFilter: function (flag) {
            var button = document.getElementById('addFromSavedFlagButton');
            var label = document.getElementById('addFromSavedFlagLabel');
            var icon = document.getElementById('addFromSavedFlagIcon');
            if (button) button.setAttribute('data-flag-id', flag.id || '');
            if (label) label.textContent = flag.name || 'All flags';
            if (icon) icon.style.color = flag.color || '#94a3b8';
            document.querySelectorAll('#addFromSavedFlagMenu .add-from-saved-flag-option').forEach(function (item) {
                var selected = (item.getAttribute('data-flag-id') || '') === (flag.id || '');
                item.setAttribute('aria-selected', selected ? 'true' : 'false');
            });
        },

        togglePickerFlagMenu: function () {
            var filter = document.getElementById('addFromSavedFlagFilter');
            var button = document.getElementById('addFromSavedFlagButton');
            if (!filter || !button) return;
            var open = filter.classList.toggle('is-open');
            button.setAttribute('aria-expanded', open ? 'true' : 'false');
        },

        closePickerFlagMenu: function () {
            var filter = document.getElementById('addFromSavedFlagFilter');
            var button = document.getElementById('addFromSavedFlagButton');
            if (filter) filter.classList.remove('is-open');
            if (button) button.setAttribute('aria-expanded', 'false');
        },

        uniqueSorted: function (values) {
            var seen = {};
            values.forEach(function (value) {
                value = (value || '').trim();
                if (value) seen[value.toLowerCase()] = value;
            });
            return Object.keys(seen).sort().map(function (key) { return seen[key]; });
        },

        selectedPickerFlag: function () {
            var button = document.getElementById('addFromSavedFlagButton');
            return button ? (button.getAttribute('data-flag-id') || '') : '';
        },

        resetPickerFilters: function () {
            var search = document.getElementById('addFromSavedSearch');
            var field = document.getElementById('addFromSavedSearchField');
            var subsystem = document.getElementById('addFromSavedSubsystemFilter');
            var onlyAddable = document.getElementById('addFromSavedOnlyAddable');
            if (search) search.value = '';
            if (field) field.value = 'all';
            if (subsystem) subsystem.value = '';
            this.setPickerFlagFilter({
                id: '',
                name: 'All flags',
                color: '#94a3b8',
            });
            if (onlyAddable) onlyAddable.checked = false;
            this.closePickerFlagMenu();
        },

        applyPickerFilters: function () {
            var termEl = document.getElementById('addFromSavedSearch');
            var fieldEl = document.getElementById('addFromSavedSearchField');
            var subsystemEl = document.getElementById('addFromSavedSubsystemFilter');
            var onlyAddableEl = document.getElementById('addFromSavedOnlyAddable');
            var term = ((termEl && termEl.value) || '').toLowerCase().trim();
            var field = (fieldEl && fieldEl.value) || 'all';
            var subsystem = ((subsystemEl && subsystemEl.value) || '').toLowerCase();
            var selectedFlag = this.selectedPickerFlag();
            var onlyAddable = !!(onlyAddableEl && onlyAddableEl.checked);
            var total = 0;
            var visible = 0;

            document.querySelectorAll('#addFromSavedList .picker-row').forEach(function (row) {
                total += 1;
                var hay = row.getAttribute('data-search-' + field) || '';
                var rowSubsystem = row.getAttribute('data-subsystem') || '';
                var rowFlags = (row.getAttribute('data-flag-ids') || '').split(',').filter(Boolean);
                var inGroup = row.getAttribute('data-in-group') === 'true';
                var matchesTerm = !term || hay.indexOf(term) !== -1;
                var matchesSubsystem = !subsystem || rowSubsystem === subsystem;
                var matchesFlags = !selectedFlag || rowFlags.indexOf(selectedFlag) !== -1;
                var matchesAddable = !onlyAddable || !inGroup;
                var show = matchesTerm && matchesSubsystem && matchesFlags && matchesAddable;
                row.style.display = show ? '' : 'none';
                if (show) visible += 1;
            });
            this.updatePickerSummary(visible, total);
        },

        filterPicker: function () {
            this.applyPickerFilters();
        },

        updatePickerSummary: function (visible, total) {
            var summary = document.getElementById('addFromSavedFilterSummary');
            if (!summary) return;
            summary.textContent = total
                ? visible + ' of ' + total + ' saved reactions shown'
                : '';
        },

        confirmPicker: function () {
            var ids = [];
            document.querySelectorAll('#addFromSavedList .picker-check').forEach(function (cb) {
                if (cb.checked && !cb.disabled) ids.push(cb.value);
            });
            if (ids.length) this.addReactions(ids);
            if (global.jQuery) jQuery('#addFromSavedModal').modal('hide');
        },
    };

    global.GroupView = GroupView;

    document.addEventListener('DOMContentLoaded', function () {
        var addBtn = document.getElementById('addFromSavedButton');
        if (addBtn) addBtn.addEventListener('click', function () { GroupView.openPicker(); });
        var confirmBtn = document.getElementById('confirmAddFromSaved');
        if (confirmBtn) confirmBtn.addEventListener('click', function () { GroupView.confirmPicker(); });
        var search = document.getElementById('addFromSavedSearch');
        if (search) search.addEventListener('input', function () { GroupView.filterPicker(this.value); });
        var flagButton = document.getElementById('addFromSavedFlagButton');
        if (flagButton) flagButton.addEventListener('click', function () { GroupView.togglePickerFlagMenu(); });
        document.addEventListener('click', function (e) {
            var flagFilter = document.getElementById('addFromSavedFlagFilter');
            if (flagFilter && !flagFilter.contains(e.target)) GroupView.closePickerFlagMenu();
        });
        document.addEventListener('keydown', function (e) {
            if (e.key === 'Escape') GroupView.closePickerFlagMenu();
        });
        [
            'addFromSavedSearchField',
            'addFromSavedSubsystemFilter',
            'addFromSavedOnlyAddable',
        ].forEach(function (id) {
            var control = document.getElementById(id);
            if (control) control.addEventListener('change', function () { GroupView.applyPickerFilters(); });
        });

        // Dirty tracking for the unsaved-edit guard (FR-017): a genuine user edit
        // anywhere in the detail panels marks the active open reaction dirty.
        // Programmatic fills during load are cleared by loadReactionById afterward.
        var workspace = document.querySelector('.main-content');
        if (workspace && window.OpenReactionState) {
            var markDirty = function (e) {
                if (e && e.isTrusted) OpenReactionState.setDirty(undefined, true);
            };
            workspace.addEventListener('input', markDirty, true);
            workspace.addEventListener('change', markDirty, true);
        }
    });
})(window);
