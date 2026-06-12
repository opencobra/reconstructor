/**
 * reactantsCompactView.js — compact (equation) presentation + quick-edit layer
 * for the Reactants panel.
 *
 * The compact view is an alternative rendering of the *existing* Reactants form,
 * which remains the single source of truth. It reads the live `.inputs-group`
 * rows to draw two equation lines (substrates --> products), and its edit modal
 * writes values back into those same row controls and re-runs the existing verify
 * flow — so the standard Create/Update Reaction path handles persistence unchanged.
 *
 * Toggle visibility / chosen view are gated on `?reaction_id=...` and persisted
 * globally in localStorage, mirroring the panel-order/size preference pattern in
 * Creatediv.js.
 */
(function (global) {
    'use strict';

    var STORAGE_KEY = 'reactantsCompactView';
    var FALLBACK_ORGANS = [
        'Adipocytes', 'Agland', 'Brain', 'Breast', 'Cervix', 'Colon', 'sIEC',
        'Uterus', 'Testis', 'Lung', 'Gall', 'Heart', 'Kidney', 'Liver',
        'Ovary', 'Pancreas', 'Pthyroidgland', 'Prostate', 'Retina',
        'Muscle', 'Skin', 'Scord', 'Spleen', 'Stomach', 'Thyroidgland',
        'Urinarybladder'
    ];

    // Side descriptors keep the substrate/product DOM names in one place.
    var SIDES = {
        substrates: { containerId: 'substratesDiv', fullname: 'substrates', sch: 'subs_sch', comps: 'subs_comps', lineId: 'compactSubstrates', nameClass: 'substrates-name' },
        products: { containerId: 'productsDiv', fullname: 'products', sch: 'prod_sch', comps: 'prod_comps', lineId: 'compactProducts', nameClass: 'products-name' },
    };

    // Row currently being edited by the modal, plus a snapshot of its values.
    var editContext = null;

    function hasReactionId() {
        return Boolean(new URLSearchParams(window.location.search).get('reaction_id'));
    }

    function isCompactPreferred() {
        return localStorage.getItem(STORAGE_KEY) === 'true';
    }

    function ensureEditOverlayInBody() {
        var overlay = document.getElementById('compactEditOverlay');
        if (overlay && overlay.parentElement !== document.body) {
            document.body.appendChild(overlay);
        }
        return overlay;
    }

    // ---- Reading a row (form is the source of truth) ---------------------

    function getIdentifierValue(group, fullname) {
        var auto = group.querySelector('.autocomplete-container');
        if (auto) {
            var vis = auto.querySelector('.autocomplete-input');
            var hidden = auto.querySelector('input[type="hidden"]');
            if (hidden && hidden.value) return hidden.value;
            return vis ? vis.value : '';
        }
        var input = group.querySelector('.cell-identifier input[name="' + fullname + '"]');
        return input ? input.value : '';
    }

    function readRow(group, side) {
        var stoichInput = group.querySelector('input[name="' + side.sch + '"]');
        var compSelect = group.querySelector('select[name="' + side.comps + '"]');
        var typeSelect = group.querySelector('select[name="' + side.fullname + '_type"]');
        var nameInput = group.querySelector('.' + side.nameClass);
        var statusInput = group.querySelector('.valid-status');
        var statusDot = group.querySelector('.status-dot');

        var found = null;
        if (statusInput) {
            found = statusInput.value === 'true';
        } else if (statusDot && statusDot.classList.contains('found')) {
            found = true;
        } else if (statusDot && statusDot.classList.contains('not-found')) {
            found = false;
        }

        return {
            group: group,
            side: side,
            stoich: stoichInput ? stoichInput.value : '1',
            identifier: getIdentifierValue(group, side.fullname),
            compartment: compSelect ? compSelect.value : '-',
            type: typeSelect ? typeSelect.value : 'VMH',
            name: nameInput ? nameInput.value : '',
            found: found,
        };
    }

    function activeRows(side) {
        var container = document.getElementById(side.containerId);
        if (!container) return [];
        var groups = Array.from(container.querySelectorAll('.inputs-group'));
        if (typeof reactionRowHasMetabolite === 'function') {
            groups = groups.filter(reactionRowHasMetabolite);
        }
        return groups;
    }

    // ---- Rendering -------------------------------------------------------

    function buildToken(rowData) {
        var token = document.createElement('span');
        token.className = 'metab-token';
        if (rowData.found === false) {
            token.classList.add('metab-token-not-found');
        }
        token._rowData = rowData; // per-token data → no cross-metabolite bleed

        var stoichInput = document.createElement('input');
        stoichInput.type = 'number';
        stoichInput.min = '1';
        stoichInput.step = 'any';
        stoichInput.className = 'compact-stoich-input';
        stoichInput.value = rowData.stoich || '1';
        stoichInput.setAttribute('aria-label', 'Stoichiometry for ' + (rowData.identifier || 'metabolite'));
        stoichInput._rowData = rowData;
        token.appendChild(stoichInput);

        var editButton = document.createElement('button');
        editButton.type = 'button';
        editButton.className = 'compact-token-edit';
        editButton.textContent = rowData.identifier + '[' + rowData.compartment + ']';
        editButton.title = 'Edit metabolite';
        editButton._rowData = rowData;
        token.appendChild(editButton);

        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'compact-token-remove';
        removeButton.setAttribute('aria-label', 'Remove ' + (rowData.identifier || 'metabolite'));
        removeButton.title = 'Remove metabolite';
        removeButton.innerHTML = '<i class="fas fa-times" aria-hidden="true"></i>';
        removeButton._rowData = rowData;
        token.appendChild(removeButton);

        return token;
    }

    function renderLine(lineEl, side) {
        if (!lineEl) return;
        lineEl.innerHTML = '';
        var rows = activeRows(side);
        if (!rows.length) {
            var empty = document.createElement('span');
            empty.className = 'compact-empty';
            empty.textContent = '—';
            lineEl.appendChild(empty);
            return;
        }
        rows.forEach(function (group, index) {
            if (index > 0) {
                var plus = document.createElement('span');
                plus.className = 'compact-plus';
                plus.textContent = '+';
                lineEl.appendChild(plus);
            }
            lineEl.appendChild(buildToken(readRow(group, side)));
        });
    }

    function readDirection() {
        var directionEl = document.getElementById('reactionDirection');
        var raw = directionEl ? directionEl.value : 'forward';
        if (typeof normalizeDirectionForSelect === 'function') {
            return normalizeDirectionForSelect(raw);
        }
        return raw === 'bidirectional' ? 'bidirectional' : 'forward';
    }

    function readOrgans() {
        return readOrganList().join(', ');
    }

    function readOrganList() {
        var container = document.getElementById('organTags');
        if (!container) return [];
        return Array.from(container.querySelectorAll('.tag'))
            .map(function (tag) { return tag.firstChild ? tag.firstChild.textContent.trim() : ''; })
            .filter(Boolean)
            .filter(function (value, index, self) { return self.indexOf(value) === index; });
    }

    function splitOrganInput(value) {
        return String(value || '')
            .split(/[,\n;]/)
            .map(function (item) { return item.trim(); })
            .filter(Boolean)
            .filter(function (value, index, self) { return self.indexOf(value) === index; });
    }

    function dispatchFieldChange(el) {
        if (!el) return;
        el.dispatchEvent(new Event('input', { bubbles: true }));
        el.dispatchEvent(new Event('change', { bubbles: true }));
    }

    function createOrganTag(tagText, container) {
        var tag = document.createElement('span');
        tag.className = 'tag';
        tag.contentEditable = false;
        tag.textContent = tagText;

        var closeBtn = document.createElement('span');
        closeBtn.className = 'close-btn';
        closeBtn.textContent = '×';
        closeBtn.addEventListener('click', function () {
            if (tag.parentNode) tag.parentNode.removeChild(tag);
            dispatchFieldChange(container);
            syncCompactMetaControls();
        });

        tag.appendChild(closeBtn);
        container.appendChild(tag);
    }

    function setOriginalOrgans(organs) {
        var container = document.getElementById('organTags');
        if (!container) return;
        container.innerHTML = '';
        organs.forEach(function (organ) { createOrganTag(organ, container); });
        dispatchFieldChange(container);
    }

    function setOriginalSubsystem(value) {
        var subsystem = document.getElementById('subsystemField');
        if (!subsystem) return;
        subsystem.value = value;
        dispatchFieldChange(subsystem);
    }

    function setOriginalDirection(value) {
        var directionEl = document.getElementById('reactionDirection');
        if (!directionEl) return;
        directionEl.value = value === 'bidirectional' ? 'bidirectional' : 'forward';
        directionEl.dispatchEvent(new Event('change', { bubbles: true }));
    }

    function getSubsystemOptions() {
        return Array.isArray(global.subsystemList) ? global.subsystemList : [];
    }

    function getOrganOptions() {
        return Array.isArray(global.ReactantsOrganList) ? global.ReactantsOrganList : FALLBACK_ORGANS;
    }

    function filterOptions(options, query) {
        var normalized = String(query || '').trim().toLowerCase();
        var source = normalized
            ? options.filter(function (item) { return String(item).toLowerCase().includes(normalized); })
            : options;
        return source.slice(0, 40);
    }

    function showCompactDropdown(dropdown, options, onSelect) {
        if (!dropdown) return;
        dropdown.innerHTML = '';
        if (!options.length) {
            dropdown.style.display = 'none';
            return;
        }
        options.forEach(function (item) {
            var option = document.createElement('button');
            option.type = 'button';
            option.className = 'compact-meta-option';
            option.textContent = item;
            option.addEventListener('mousedown', function (e) {
                e.preventDefault();
                onSelect(item);
                hideCompactDropdown(dropdown);
            });
            dropdown.appendChild(option);
        });
        dropdown.style.display = 'block';
    }

    function hideCompactDropdown(dropdown) {
        if (!dropdown) return;
        dropdown.style.display = 'none';
    }

    function renderSubsystemDropdown() {
        var input = document.getElementById('compactSubsystemInput');
        var dropdown = document.getElementById('compactSubsystemDropdown');
        if (!input || !dropdown) return;
        showCompactDropdown(dropdown, filterOptions(getSubsystemOptions(), input.value), function (item) {
            input.value = item;
            setOriginalSubsystem(item);
        });
    }

    function currentOrganQuery(value) {
        var parts = String(value || '').split(/[,\n;]/);
        return parts[parts.length - 1] || '';
    }

    function selectCompactOrgan(item) {
        var input = document.getElementById('compactOrganInput');
        if (!input) return;
        var parts = String(input.value || '').split(/[,\n;]/);
        parts.pop();
        var organs = parts.map(function (part) { return part.trim(); }).filter(Boolean);
        organs.push(item);
        organs = organs.filter(function (value, index, self) { return self.indexOf(value) === index; });
        setOriginalOrgans(organs);
        input.value = readOrgans();
    }

    function renderOrganDropdown() {
        var input = document.getElementById('compactOrganInput');
        var dropdown = document.getElementById('compactOrganDropdown');
        if (!input || !dropdown) return;
        showCompactDropdown(dropdown, filterOptions(getOrganOptions(), currentOrganQuery(input.value)), selectCompactOrgan);
    }

    function syncCompactMetaControls() {
        var subsystem = document.getElementById('subsystemField');
        var compactSubsystem = document.getElementById('compactSubsystemInput');
        if (compactSubsystem && document.activeElement !== compactSubsystem) {
            compactSubsystem.value = (subsystem && subsystem.value.trim()) || '';
        }

        var compactOrgan = document.getElementById('compactOrganInput');
        if (compactOrgan && document.activeElement !== compactOrgan) {
            compactOrgan.value = readOrgans();
        }

        var compactDirection = document.getElementById('compactDirectionSelect');
        if (compactDirection && document.activeElement !== compactDirection) {
            compactDirection.value = readDirection();
        }
    }

    function render() {
        var view = document.getElementById('reactantsCompactView');
        if (!view) return;

        renderLine(document.getElementById('compactSubstrates'), SIDES.substrates);
        renderLine(document.getElementById('compactProducts'), SIDES.products);

        var direction = readDirection();
        var arrow = document.getElementById('compactArrow');
        if (arrow) {
            arrow.textContent = direction === 'bidirectional' ? '<-->' : '-->';
        }

        syncCompactMetaControls();
    }

    function setRowStoich(rowData, value, shouldRender) {
        if (!rowData || !rowData.group) return;
        var stoichInput = rowData.group.querySelector('input[name="' + rowData.side.sch + '"]');
        if (!stoichInput) return;
        stoichInput.value = value;
        stoichInput.dispatchEvent(new Event('input', { bubbles: true }));
        if (shouldRender) render();
    }

    function removeRow(group) {
        if (!group) return;
        group.remove();
        if (typeof updateAtomChargeCounters === 'function') {
            updateAtomChargeCounters();
        }
        if (global.ReactantsFormDirty) {
            ReactantsFormDirty.scheduleEvaluate();
        }
        render();
    }

    async function confirmAndRemoveRow(rowData) {
        if (!rowData || !rowData.group) return;
        var label = rowData.identifier ? rowData.identifier + '[' + rowData.compartment + ']' : 'this metabolite';
        var confirmed = true;
        if (global.Notify && typeof Notify.confirm === 'function') {
            confirmed = await Notify.confirm({
                title: 'Remove metabolite',
                message: 'Remove ' + label + ' from the reaction?',
                confirmText: 'Remove',
                cancelText: 'Cancel',
                danger: true,
            });
        }
        if (confirmed) {
            removeRow(rowData.group);
        }
    }

    function addMetabolite(side) {
        var container = document.getElementById(side.containerId);
        if (!container) return;

        if (typeof addField === 'function') {
            addField(side.containerId, side.fullname, side.sch);
        } else {
            var fallbackButton = document.getElementById(side.containerId === 'substratesDiv' ? 'addSubstrate' : 'addProduct');
            if (fallbackButton) fallbackButton.click();
        }

        var groups = container.querySelectorAll('.inputs-group');
        var group = groups[groups.length - 1];
        if (!group) return;
        openEditModal(readRow(group, side), { isNew: true });
        render();
    }

    // ---- Hover tooltip ---------------------------------------------------

    var activeTip = null;

    function removeTip() {
        if (activeTip && activeTip.parentNode) activeTip.parentNode.removeChild(activeTip);
        activeTip = null;
    }

    function showTip(token) {
        removeTip();
        var data = token._rowData;
        if (!data) return;
        var tip = document.createElement('div');
        tip.className = 'compact-token-tip';
        var foundLabel = data.found === true ? 'Found in VMH'
            : data.found === false ? 'Not found in VMH' : 'Not verified';
        tip.innerHTML =
            '<div class="ctt-row"><span class="ctt-key">Type</span><span class="ctt-val">' + escapeHtml(data.type || '—') + '</span></div>' +
            '<div class="ctt-row"><span class="ctt-key">Name</span><span class="ctt-val">' + escapeHtml(data.name || '—') + '</span></div>' +
            '<div class="ctt-row"><span class="ctt-key">Status</span><span class="ctt-val">' + foundLabel + '</span></div>';
        document.body.appendChild(tip);
        activeTip = tip;

        var rect = token.getBoundingClientRect();
        var tipRect = tip.getBoundingClientRect();
        var top = rect.bottom + 8;
        var left = rect.left;
        if (left + tipRect.width > window.innerWidth - 12) {
            left = Math.max(8, window.innerWidth - tipRect.width - 12);
        }
        tip.style.top = top + 'px';
        tip.style.left = left + 'px';
    }

    function escapeHtml(str) {
        var div = document.createElement('div');
        div.textContent = str == null ? '' : String(str);
        return div.innerHTML;
    }

    // ---- Edit modal ------------------------------------------------------

    function openEditModal(rowData, options) {
        options = options || {};
        var overlay = ensureEditOverlayInBody();
        if (!overlay) return;
        editContext = {
            group: rowData.group,
            side: rowData.side,
            original: { identifier: rowData.identifier, type: rowData.type, compartment: rowData.compartment, stoich: rowData.stoich },
            isNew: !!options.isNew,
        };
        document.getElementById('compactEditIdentifier').value = rowData.identifier;
        document.getElementById('compactEditType').value = rowData.type;
        document.getElementById('compactEditStoich').value = rowData.stoich;
        document.getElementById('compactEditCompartment').value = rowData.compartment;
        var title = document.getElementById('compactEditTitle');
        if (title) title.textContent = options.isNew ? 'Add metabolite' : 'Edit metabolite';
        var saveButton = document.getElementById('compactEditSave');
        if (saveButton) saveButton.textContent = options.isNew ? 'Add' : 'Save';
        overlay.classList.add('open');
        var modal = overlay.querySelector('.confirm-modal');
        requestAnimationFrame(function () {
            if (modal) modal.classList.add('confirm-modal-in');
            document.getElementById('compactEditIdentifier').focus();
        });
    }

    function closeEditModal(options) {
        options = options || {};
        var context = editContext;
        var overlay = ensureEditOverlayInBody();
        if (overlay) {
            overlay.classList.remove('open');
            var modal = overlay.querySelector('.confirm-modal');
            if (modal) modal.classList.remove('confirm-modal-in');
        }
        editContext = null;
        if (options.discardNew && context && context.isNew) {
            removeRow(context.group);
        }
    }

    function setIdentifierValue(group, fullname, value) {
        var auto = group.querySelector('.autocomplete-container');
        if (auto) {
            var vis = auto.querySelector('.autocomplete-input');
            var hidden = auto.querySelector('input[type="hidden"]');
            if (vis) { vis.disabled = false; vis.value = value; }
            if (hidden) hidden.value = value;
            dispatchFieldChange(vis);
            dispatchFieldChange(hidden);
            return;
        }
        var input = group.querySelector('.cell-identifier input[name="' + fullname + '"]');
        if (input) {
            input.disabled = false;
            input.value = value;
            dispatchFieldChange(input);
        }
    }

    // Re-render whenever the verify flow mutates the edited row (async fetch).
    function watchRowAndRerender(group) {
        if (!window.MutationObserver) {
            setTimeout(render, 1500);
            return;
        }
        var observer = new MutationObserver(function () { render(); });
        observer.observe(group, { childList: true, subtree: true, attributes: true });
        setTimeout(function () { observer.disconnect(); render(); }, 8000);
    }

    function saveEdit() {
        if (!editContext) return;
        var group = editContext.group;
        var side = editContext.side;
        var orig = editContext.original;

        var newIdentifier = document.getElementById('compactEditIdentifier').value.trim();
        var newType = document.getElementById('compactEditType').value;
        var newStoich = document.getElementById('compactEditStoich').value;
        var newComp = document.getElementById('compactEditCompartment').value;

        if (!newIdentifier) {
            Notify.error('Enter a metabolite identifier.');
            return;
        }
        if (!newStoich || Number(newStoich) <= 0) {
            Notify.error('Enter a positive stoichiometry.');
            return;
        }
        if (newComp === '-') {
            Notify.error('Select a compartment.');
            return;
        }

        var stoichInput = group.querySelector('input[name="' + side.sch + '"]');
        var compSelect = group.querySelector('select[name="' + side.comps + '"]');
        var typeSelect = group.querySelector('select[name="' + side.fullname + '_type"]');

        var identityChanged = newIdentifier !== orig.identifier || newType !== orig.type || newComp !== orig.compartment;

        // Apply input type first so the identifier cell matches (text vs file vs Saved).
        if (typeSelect) {
            typeSelect.disabled = false;
            if (typeSelect.value !== newType) {
                typeSelect.value = newType;
                if (typeof toggleFileInput === 'function') toggleFileInput(group, newType);
                if (typeof handleMetaboliteTypeChange === 'function') handleMetaboliteTypeChange(typeSelect);
            }
        }
        setIdentifierValue(group, side.fullname, newIdentifier);
        if (compSelect) compSelect.value = newComp;
        if (stoichInput) stoichInput.value = newStoich;
        dispatchFieldChange(typeSelect);
        dispatchFieldChange(compSelect);
        dispatchFieldChange(stoichInput);

        closeEditModal({ discardNew: false });

        if (identityChanged) {
            // Clear stale verification state, then re-run the row's existing verify
            // so name / found-status / atom counts refresh (auto re-verify on save).
            group.querySelectorAll('.valid-status').forEach(function (el) { el.remove(); });
            delete group.dataset.atomCounts;
            delete group.dataset.charge;
            var statusDot = group.querySelector('.status-dot');
            if (statusDot) statusDot.style.display = 'none';
            var nameInput = group.querySelector('.' + side.nameClass);
            if (nameInput) nameInput.disabled = false;

            var doneBtn = group.querySelector('.done-field-btn');
            if (doneBtn && typeof handleDoneButtonClick === 'function') {
                doneBtn.style.display = '';
                handleDoneButtonClick({ target: doneBtn });
                watchRowAndRerender(group);
            }
        } else if (stoichInput) {
            // Stoichiometry-only change: let the delegated listener recompute balance.
            stoichInput.dispatchEvent(new Event('input', { bubbles: true }));
        }

        render();
    }

    // ---- View toggling ---------------------------------------------------

    function applyView(isCompact, persistPreference) {
        var container = document.getElementById('reactionFormContainer');
        if (!container) return;
        var toolbar = container.querySelector('.reactants-toolbar');
        var form = document.getElementById('reactionForm');
        var compact = document.getElementById('reactantsCompactView');

        if (isCompact) {
            if (toolbar) toolbar.style.display = 'none';
            if (form) form.style.display = 'none';
            if (compact) compact.style.display = 'block';
            render();
        } else {
            if (toolbar) toolbar.style.display = '';
            if (form) form.style.display = '';
            if (compact) compact.style.display = 'none';
            removeTip();
        }
        var toggle = document.getElementById('reactantsCompactToggle');
        if (toggle) toggle.checked = isCompact;
        if (persistPreference !== false) {
            localStorage.setItem(STORAGE_KEY, isCompact ? 'true' : 'false');
        }
    }

    function syncToggleVisibility() {
        var wrap = document.getElementById('reactantsViewToggle');
        if (!wrap) return false;
        var show = hasReactionId();
        wrap.style.display = show ? '' : 'none';
        if (!show) {
            // New/uncreated reaction: never present the compact view, but keep the
            // saved preference intact for the next created reaction.
            applyView(false, false);
        }
        return show;
    }

    // ---- Wiring ----------------------------------------------------------

    function bindEvents() {
        var view = document.getElementById('reactantsCompactView');
        if (view) {
            view.addEventListener('mouseover', function (e) {
                var token = e.target.closest('.metab-token');
                if (token && !e.target.closest('.compact-stoich-input, .compact-token-remove')) showTip(token);
            });
            view.addEventListener('mouseout', function (e) {
                if (e.target.closest('.metab-token')) removeTip();
            });
            view.addEventListener('input', function (e) {
                var stoichInput = e.target.closest('.compact-stoich-input');
                if (stoichInput && stoichInput._rowData) {
                    setRowStoich(stoichInput._rowData, stoichInput.value, false);
                }
            });
            view.addEventListener('change', function (e) {
                var stoichInput = e.target.closest('.compact-stoich-input');
                if (!stoichInput || !stoichInput._rowData) return;
                if (!stoichInput.value || Number(stoichInput.value) <= 0) {
                    stoichInput.value = '1';
                    Notify.error('Enter a positive stoichiometry.');
                }
                setRowStoich(stoichInput._rowData, stoichInput.value, true);
            });
            view.addEventListener('click', async function (e) {
                var removeButton = e.target.closest('.compact-token-remove');
                if (removeButton && removeButton._rowData) {
                    e.preventDefault();
                    e.stopPropagation();
                    removeTip();
                    await confirmAndRemoveRow(removeButton._rowData);
                    return;
                }

                var editButton = e.target.closest('.compact-token-edit');
                if (editButton && editButton._rowData) {
                    e.preventDefault();
                    e.stopPropagation();
                    removeTip();
                    openEditModal(editButton._rowData);
                    return;
                }

                var token = e.target.closest('.metab-token');
                if (token && token._rowData && !e.target.closest('.compact-stoich-input, .compact-token-remove')) {
                    e.preventDefault();
                    removeTip();
                    openEditModal(token._rowData);
                }
            });
            view.addEventListener('keydown', function (e) {
                if ((e.key === 'Enter' || e.key === ' ') && e.target.closest('.compact-token-edit')) {
                    var editButton = e.target.closest('.compact-token-edit');
                    if (editButton && editButton._rowData) {
                        e.preventDefault();
                        openEditModal(editButton._rowData);
                    }
                }
            });
        }

        var toggle = document.getElementById('reactantsCompactToggle');
        if (toggle) {
            toggle.addEventListener('change', function () { applyView(toggle.checked); });
        }

        var compactAddSubstrate = document.getElementById('compactAddSubstrate');
        if (compactAddSubstrate) {
            compactAddSubstrate.addEventListener('click', function () { addMetabolite(SIDES.substrates); });
        }

        var compactAddProduct = document.getElementById('compactAddProduct');
        if (compactAddProduct) {
            compactAddProduct.addEventListener('click', function () { addMetabolite(SIDES.products); });
        }

        var compactSubsystem = document.getElementById('compactSubsystemInput');
        if (compactSubsystem) {
            compactSubsystem.addEventListener('input', function () {
                setOriginalSubsystem(compactSubsystem.value);
                renderSubsystemDropdown();
            });
            compactSubsystem.addEventListener('focus', renderSubsystemDropdown);
            compactSubsystem.addEventListener('change', function () {
                compactSubsystem.value = compactSubsystem.value.trim();
                setOriginalSubsystem(compactSubsystem.value);
            });
        }

        var compactOrgan = document.getElementById('compactOrganInput');
        if (compactOrgan) {
            compactOrgan.addEventListener('input', function () {
                setOriginalOrgans(splitOrganInput(compactOrgan.value));
                renderOrganDropdown();
            });
            compactOrgan.addEventListener('focus', renderOrganDropdown);
            compactOrgan.addEventListener('change', function () {
                setOriginalOrgans(splitOrganInput(compactOrgan.value));
                compactOrgan.value = readOrgans();
            });
            compactOrgan.addEventListener('keydown', function (e) {
                if (e.key === 'Enter') {
                    e.preventDefault();
                    setOriginalOrgans(splitOrganInput(compactOrgan.value));
                    compactOrgan.value = readOrgans();
                    compactOrgan.blur();
                }
            });
        }

        document.addEventListener('click', function (e) {
            var subsystemWrap = document.getElementById('compactSubsystemInput')?.closest('.compact-meta-input');
            var organWrap = document.getElementById('compactOrganInput')?.closest('.compact-meta-input');
            if (!subsystemWrap || !subsystemWrap.contains(e.target)) {
                hideCompactDropdown(document.getElementById('compactSubsystemDropdown'));
            }
            if (!organWrap || !organWrap.contains(e.target)) {
                hideCompactDropdown(document.getElementById('compactOrganDropdown'));
            }
        });

        var compactDirection = document.getElementById('compactDirectionSelect');
        if (compactDirection) {
            compactDirection.addEventListener('change', function () {
                setOriginalDirection(compactDirection.value);
                render();
            });
        }

        var cancelBtn = document.getElementById('compactEditCancel');
        var saveBtn = document.getElementById('compactEditSave');
        var overlay = ensureEditOverlayInBody();
        if (cancelBtn) cancelBtn.addEventListener('click', function () { closeEditModal({ discardNew: true }); });
        if (saveBtn) saveBtn.addEventListener('click', saveEdit);
        if (overlay) {
            overlay.addEventListener('click', function (e) { if (e.target === overlay) closeEditModal({ discardNew: true }); });
        }
        document.addEventListener('keydown', function (e) {
            if (e.key === 'Escape' && editContext) closeEditModal({ discardNew: true });
        });
    }

    function init() {
        ensureEditOverlayInBody();
        bindEvents();
        var showToggle = syncToggleVisibility();
        if (showToggle) {
            applyView(isCompactPreferred());
        }
    }

    // Called by displayDivs / loadReactionById once reaction data is loaded, so
    // the compact lines reflect freshly populated rows and tab switches.
    function refresh() {
        var showToggle = syncToggleVisibility();
        if (showToggle && isCompactPreferred()) {
            applyView(true);
        }
    }

    global.ReactantsCompactView = { init: init, refresh: refresh, render: render };

    document.addEventListener('DOMContentLoaded', init);
})(window);
