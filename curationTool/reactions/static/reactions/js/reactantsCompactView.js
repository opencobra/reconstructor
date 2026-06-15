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
 * The toggle is always available (both new and existing reactions), since the
 * compact view operates entirely on the live form rows and needs no saved
 * reaction. The chosen view is persisted globally in localStorage, mirroring the
 * panel-order/size preference pattern in Creatediv.js.
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
    var editSaving = false;
    var compactSavedMetabolitesPromise = null;

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

    function savedMetaboliteLabel(metabolite) {
        if (!metabolite) return '';
        var name = metabolite.name || '';
        var abbr = metabolite.vmh_abbr || '';
        if (abbr && name) return abbr + ' - ' + name;
        return name || abbr || String(metabolite.id || '');
    }

    function findSavedMetaboliteById(id) {
        if (!id || !Array.isArray(global.savedMetabolitesCache)) return null;
        return global.savedMetabolitesCache.find(function (metabolite) {
            return String(metabolite.id) === String(id);
        }) || null;
    }

    function getIdentifierDisplayValue(group, fullname, type, fallbackName) {
        if (type === 'Saved') {
            var auto = group.querySelector('.autocomplete-container');
            if (auto) {
                var visible = auto.querySelector('.autocomplete-input');
                if (visible && visible.value) return visible.value;
            }
            var savedId = getIdentifierValue(group, fullname);
            return fallbackName || savedMetaboliteLabel(findSavedMetaboliteById(savedId)) || savedId;
        }

        if (type === 'MDL Mol file') {
            var fileInput = group.querySelector('.cell-identifier input[type="file"]');
            if (fileInput && fileInput.files && fileInput.files.length > 0) {
                return fileInput.files[0].name;
            }
            var raw = getIdentifierValue(group, fullname);
            return raw ? raw.split('/').pop() : '';
        }

        return getIdentifierValue(group, fullname);
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

        var type = typeSelect ? typeSelect.value : 'VMH';
        var identifier = getIdentifierValue(group, side.fullname);

        return {
            group: group,
            side: side,
            stoich: stoichInput ? stoichInput.value : '1',
            identifier: identifier,
            displayIdentifier: getIdentifierDisplayValue(group, side.fullname, type, nameInput ? nameInput.value : ''),
            compartment: compSelect ? compSelect.value : '-',
            type: type,
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
        stoichInput.setAttribute('aria-label', 'Stoichiometry for ' + (rowData.displayIdentifier || rowData.identifier || 'metabolite'));
        stoichInput._rowData = rowData;
        token.appendChild(stoichInput);

        var editButton = document.createElement('button');
        editButton.type = 'button';
        editButton.className = 'compact-token-edit';
        editButton.textContent = (rowData.displayIdentifier || rowData.identifier) + '[' + rowData.compartment + ']';
        editButton.title = 'Edit metabolite';
        editButton._rowData = rowData;
        token.appendChild(editButton);

        var removeButton = document.createElement('button');
        removeButton.type = 'button';
        removeButton.className = 'compact-token-remove';
        removeButton.setAttribute('aria-label', 'Remove ' + (rowData.displayIdentifier || rowData.identifier || 'metabolite'));
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

    function setCompactEditError(message) {
        var error = document.getElementById('compactEditError');
        if (!error) return;
        error.textContent = message || '';
        error.style.display = message ? '' : 'none';
    }

    function setCompactSaveBusy(isBusy) {
        editSaving = !!isBusy;
        var saveButton = document.getElementById('compactEditSave');
        var cancelButton = document.getElementById('compactEditCancel');
        if (saveButton) {
            saveButton.disabled = !!isBusy;
            saveButton.textContent = isBusy ? 'Verifying...' : (editContext && editContext.isNew ? 'Add' : 'Save');
        }
        if (cancelButton) {
            cancelButton.disabled = !!isBusy;
        }
    }

    function compactSavedSearchValue(metaboliteId, fallback) {
        return savedMetaboliteLabel(findSavedMetaboliteById(metaboliteId)) || fallback || '';
    }

    function ensureSavedMetabolites() {
        if (Array.isArray(global.savedMetabolitesCache)) {
            return Promise.resolve(global.savedMetabolitesCache);
        }
        if (compactSavedMetabolitesPromise) {
            return compactSavedMetabolitesPromise;
        }
        var userId = sessionStorage.getItem('userID');
        compactSavedMetabolitesPromise = fetch('/get_saved_metabolites/?user_id=' + encodeURIComponent(userId || ''))
            .then(function (response) { return response.json(); })
            .then(function (data) {
                global.savedMetabolitesCache = data.metabolites || [];
                return global.savedMetabolitesCache;
            })
            .catch(function () {
                global.savedMetabolitesCache = [];
                return global.savedMetabolitesCache;
            });
        return compactSavedMetabolitesPromise;
    }

    function hideCompactSavedDropdown() {
        var dropdown = document.getElementById('compactEditSavedDropdown');
        if (dropdown) dropdown.style.display = 'none';
    }

    function chooseCompactSavedMetabolite(metabolite) {
        var search = document.getElementById('compactEditSavedSearch');
        var hidden = document.getElementById('compactEditSavedId');
        if (search) search.value = savedMetaboliteLabel(metabolite);
        if (hidden) hidden.value = metabolite && metabolite.id != null ? metabolite.id : '';
        hideCompactSavedDropdown();
        setCompactEditError('');
    }

    function renderCompactSavedDropdown() {
        var search = document.getElementById('compactEditSavedSearch');
        var dropdown = document.getElementById('compactEditSavedDropdown');
        if (!search || !dropdown) return;
        ensureSavedMetabolites().then(function (metabolites) {
            var query = String(search.value || '').trim().toLowerCase();
            var matches = (metabolites || []).filter(function (metabolite) {
                return !query ||
                    String(metabolite.name || '').toLowerCase().includes(query) ||
                    String(metabolite.vmh_abbr || '').toLowerCase().includes(query) ||
                    String(metabolite.id || '').toLowerCase() === query;
            }).slice(0, 40);

            dropdown.innerHTML = '';
            matches.forEach(function (metabolite) {
                var option = document.createElement('button');
                option.type = 'button';
                option.className = 'compact-edit-saved-option';
                option.textContent = savedMetaboliteLabel(metabolite);
                option.addEventListener('mousedown', function (event) {
                    event.preventDefault();
                    chooseCompactSavedMetabolite(metabolite);
                });
                dropdown.appendChild(option);
            });
            dropdown.style.display = matches.length ? 'block' : 'none';
        });
    }

    function inferSavedSelectionFromSearch() {
        var search = document.getElementById('compactEditSavedSearch');
        var hidden = document.getElementById('compactEditSavedId');
        if (!search || !hidden || hidden.value) return;
        var query = String(search.value || '').trim().toLowerCase();
        if (!query || !Array.isArray(global.savedMetabolitesCache)) return;
        var match = global.savedMetabolitesCache.find(function (metabolite) {
            return String(metabolite.name || '').toLowerCase() === query ||
                String(metabolite.vmh_abbr || '').toLowerCase() === query ||
                savedMetaboliteLabel(metabolite).toLowerCase() === query;
        });
        if (match) {
            chooseCompactSavedMetabolite(match);
        }
    }

    function updateCompactFileName() {
        var fileInput = document.getElementById('compactEditFile');
        var fileName = document.getElementById('compactEditFileName');
        if (!fileName) return;
        var existingName = editContext && editContext.original ? editContext.original.fileName : '';
        if (fileInput && fileInput.files && fileInput.files.length > 0) {
            fileName.textContent = fileInput.files[0].name;
        } else {
            fileName.textContent = existingName ? 'Current file: ' + existingName : '';
        }
    }

    function setCompactIdentifierMode(type, rowData) {
        var label = document.getElementById('compactEditIdentifierLabel');
        var textInput = document.getElementById('compactEditIdentifier');
        var savedContainer = document.getElementById('compactEditSavedContainer');
        var savedSearch = document.getElementById('compactEditSavedSearch');
        var savedId = document.getElementById('compactEditSavedId');
        var fileContainer = document.getElementById('compactEditFileContainer');
        var fileInput = document.getElementById('compactEditFile');

        if (label) label.textContent = type === 'Saved' ? 'My metabolite' : (type === 'MDL Mol file' ? 'MOL file' : 'Metabolite identifier');
        if (textInput) textInput.style.display = type === 'Saved' || type === 'MDL Mol file' ? 'none' : '';
        if (savedContainer) savedContainer.style.display = type === 'Saved' ? '' : 'none';
        if (fileContainer) fileContainer.style.display = type === 'MDL Mol file' ? '' : 'none';
        hideCompactSavedDropdown();

        if (type === 'Saved') {
            var id = rowData ? rowData.identifier : '';
            var display = rowData ? (rowData.displayIdentifier || rowData.name || '') : '';
            if (savedId) savedId.value = id || '';
            if (savedSearch) savedSearch.value = compactSavedSearchValue(id, display);
            ensureSavedMetabolites().then(function () {
                if (savedSearch && savedId) {
                    savedSearch.value = compactSavedSearchValue(savedId.value, savedSearch.value);
                }
            });
        } else if (type === 'MDL Mol file') {
            if (fileInput) fileInput.value = '';
            updateCompactFileName();
        }
    }

    function selectedCompactIdentifier(type) {
        if (type === 'Saved') {
            inferSavedSelectionFromSearch();
            var savedId = document.getElementById('compactEditSavedId');
            return savedId ? savedId.value.trim() : '';
        }
        if (type === 'MDL Mol file') {
            var fileInput = document.getElementById('compactEditFile');
            if (fileInput && fileInput.files && fileInput.files.length > 0) {
                return fileInput.files[0].name;
            }
            return editContext && editContext.original ? editContext.original.identifier : '';
        }
        var textInput = document.getElementById('compactEditIdentifier');
        return textInput ? textInput.value.trim() : '';
    }

    function selectedCompactFile() {
        var fileInput = document.getElementById('compactEditFile');
        return fileInput && fileInput.files && fileInput.files.length > 0 ? fileInput.files[0] : null;
    }

    // ---- Edit modal ------------------------------------------------------

    function openEditModal(rowData, options) {
        options = options || {};
        var overlay = ensureEditOverlayInBody();
        if (!overlay) return;
        var modalType = rowData.type === 'Draw' ? 'VMH' : rowData.type;
        editContext = {
            group: rowData.group,
            side: rowData.side,
            original: {
                identifier: rowData.identifier,
                displayIdentifier: rowData.displayIdentifier,
                type: rowData.type,
                compartment: rowData.compartment,
                stoich: rowData.stoich,
                fileName: rowData.type === 'MDL Mol file' ? (rowData.displayIdentifier || rowData.identifier || '') : '',
            },
            isNew: !!options.isNew,
        };
        setCompactEditError('');
        document.getElementById('compactEditIdentifier').value = modalType === 'Saved' ? '' : (rowData.displayIdentifier || rowData.identifier || '');
        document.getElementById('compactEditType').value = modalType;
        document.getElementById('compactEditStoich').value = rowData.stoich;
        document.getElementById('compactEditCompartment').value = rowData.compartment;
        setCompactIdentifierMode(modalType, rowData);
        var title = document.getElementById('compactEditTitle');
        if (title) title.textContent = options.isNew ? 'Add metabolite' : 'Edit metabolite';
        var saveButton = document.getElementById('compactEditSave');
        if (saveButton) saveButton.textContent = options.isNew ? 'Add' : 'Save';
        setCompactSaveBusy(false);
        overlay.classList.add('open');
        var modal = overlay.querySelector('.confirm-modal');
        requestAnimationFrame(function () {
            if (modal) modal.classList.add('confirm-modal-in');
            var focusTarget = modalType === 'Saved'
                ? document.getElementById('compactEditSavedSearch')
                : modalType === 'MDL Mol file'
                    ? document.getElementById('compactEditFile')
                    : document.getElementById('compactEditIdentifier');
            if (focusTarget) focusTarget.focus();
        });
    }

    function closeEditModal(options) {
        options = options || {};
        if (editSaving) return;
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

    function sourceFileInput(group, side) {
        var fileInput = group.querySelector('.cell-identifier input[type="file"]');
        if (!fileInput) {
            var identifierCell = group.querySelector('.cell-identifier');
            var textInput = group.querySelector('.cell-identifier input[name="' + side.fullname + '"]');
            fileInput = document.createElement('input');
            fileInput.type = 'file';
            fileInput.name = textInput ? textInput.name : side.fullname;
            fileInput.id = textInput ? textInput.id : side.fullname;
            if (identifierCell) identifierCell.appendChild(fileInput);
        }
        return fileInput;
    }

    function copyFileToInput(file, fileInput) {
        if (!file || !fileInput || typeof DataTransfer === 'undefined') return;
        var transfer = new DataTransfer();
        transfer.items.add(file);
        fileInput.files = transfer.files;
    }

    function clearRowVerification(group) {
        group.querySelectorAll('.valid-status').forEach(function (el) { el.remove(); });
        delete group.dataset.atomCounts;
        delete group.dataset.charge;
        var statusDot = group.querySelector('.status-dot');
        if (statusDot) {
            statusDot.style.display = 'none';
            statusDot.className = 'status-dot';
            statusDot.onclick = null;
            statusDot.style.cursor = 'default';
        }
    }

    function buildCompactVerifyData(group, side, values) {
        var data = new FormData();
        data.append('metabolite', values.identifier);
        data.append('type', values.type);
        data.append('compartment', values.compartment);
        data.append('stoichiometry', values.stoich);
        data.append('userID', sessionStorage.getItem('userID'));

        if (values.type === 'MDL Mol file') {
            var existingFileInput = group.querySelector('.cell-identifier input[type="file"]');
            var file = values.file || (existingFileInput && existingFileInput.files && existingFileInput.files[0]);
            if (file) {
                data.append('file', file);
            }
        }

        return data;
    }

    function verifyCompactValues(group, side, values) {
        return fetch(verifyMetabolite, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
            },
            body: buildCompactVerifyData(group, side, values),
        })
            .then(function (response) { return response.json(); })
            .then(function (data) {
                if (data.error) {
                    throw new Error(data.message || 'Verification failed.');
                }
                return data;
            });
    }

    function ensureSourceInputMode(group, side, type) {
        var typeSelect = group.querySelector('select[name="' + side.fullname + '_type"]');
        if (typeSelect) {
            typeSelect.disabled = false;
            typeSelect.value = type;
            if (typeof toggleFileInput === 'function') toggleFileInput(group, type);
            if (typeof handleMetaboliteTypeChange === 'function') handleMetaboliteTypeChange(typeSelect);
        }
        return typeSelect;
    }

    function applyVerifiedCompactEdit(group, side, values, verifyData) {
        var stoichInput = group.querySelector('input[name="' + side.sch + '"]');
        var compSelect = group.querySelector('select[name="' + side.comps + '"]');
        var typeSelect = ensureSourceInputMode(group, side, values.type);
        var mainInput = null;

        if (values.type === 'Saved') {
            var auto = group.querySelector('.autocomplete-container');
            var search = auto ? auto.querySelector('.autocomplete-input') : null;
            var hidden = auto ? auto.querySelector('input[type="hidden"]') : null;
            if (search) {
                search.disabled = false;
                search.value = values.savedLabel || values.identifier;
            }
            if (hidden) {
                hidden.value = values.identifier;
            }
            mainInput = auto;
            dispatchFieldChange(search);
            dispatchFieldChange(hidden);
        } else if (values.type === 'MDL Mol file') {
            mainInput = group.querySelector('.cell-identifier input[type="text"]:not(.autocomplete-input)');
            if (mainInput) {
                mainInput.disabled = false;
                mainInput.value = values.identifier;
            }
            var fileInput = sourceFileInput(group, side);
            if (values.file) {
                copyFileToInput(values.file, fileInput);
            }
            dispatchFieldChange(mainInput);
            dispatchFieldChange(fileInput);
        } else {
            mainInput = group.querySelector('.cell-identifier input[name="' + side.fullname + '"]');
            if (mainInput) {
                mainInput.disabled = false;
                mainInput.value = values.identifier;
            }
            dispatchFieldChange(mainInput);
        }

        if (compSelect) compSelect.value = values.compartment;
        if (stoichInput) stoichInput.value = values.stoich;
        dispatchFieldChange(typeSelect);
        dispatchFieldChange(compSelect);
        dispatchFieldChange(stoichInput);

        clearRowVerification(group);
        var doneBtn = group.querySelector('.done-field-btn');
        if (doneBtn && typeof updateNameFields === 'function') {
            updateNameFields(verifyData, mainInput, typeSelect, doneBtn);
        }

        group.dataset.atomCounts = JSON.stringify(verifyData.atom_counts || {});
        group.dataset.charge = verifyData.charge;

        var hiddenStatus = document.createElement('input');
        hiddenStatus.style.display = 'none';
        hiddenStatus.value = verifyData.found;
        hiddenStatus.className = 'valid-status';
        group.appendChild(hiddenStatus);

        if (typeof updateAtomChargeCounters === 'function') {
            updateAtomChargeCounters();
        }
        if (global.ReactantsFormDirty) {
            ReactantsFormDirty.scheduleEvaluate();
        }
        watchRowAndRerender(group);
    }

    async function saveEdit() {
        if (!editContext) return;
        var group = editContext.group;
        var side = editContext.side;
        var newType = document.getElementById('compactEditType').value;
        var newIdentifier = selectedCompactIdentifier(newType);
        var newStoich = document.getElementById('compactEditStoich').value;
        var newComp = document.getElementById('compactEditCompartment').value;
        var selectedFile = selectedCompactFile();

        if (!newIdentifier) {
            var identifierMessage = newType === 'Saved'
                ? 'Select a saved metabolite.'
                : newType === 'MDL Mol file'
                    ? 'Choose a MOL file before saving.'
                    : 'Enter a metabolite identifier.';
            setCompactEditError(identifierMessage);
            Notify.error(identifierMessage);
            return;
        }
        if (!newStoich || Number(newStoich) <= 0) {
            setCompactEditError('Enter a positive stoichiometry.');
            Notify.error('Enter a positive stoichiometry.');
            return;
        }
        if (newComp === '-') {
            setCompactEditError('Select a compartment.');
            Notify.error('Select a compartment.');
            return;
        }

        var values = {
            identifier: newIdentifier,
            savedLabel: document.getElementById('compactEditSavedSearch')?.value || '',
            type: newType,
            stoich: newStoich,
            compartment: newComp,
            file: selectedFile,
        };

        setCompactEditError('');
        setCompactSaveBusy(true);
        try {
            var verifyData = await verifyCompactValues(group, side, values);
            applyVerifiedCompactEdit(group, side, values, verifyData);
            setCompactSaveBusy(false);
            closeEditModal({ discardNew: false });
            render();
        } catch (error) {
            var message = error && error.message ? error.message : 'Verification failed.';
            setCompactEditError('Verification failed: ' + message);
            Notify.error('Verification failed: ' + message);
            setCompactSaveBusy(false);
        }
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
        // The compact view works off the live form rows, so it is available for
        // both new and existing reactions. Show the toggle whenever it exists.
        var wrap = document.getElementById('reactantsViewToggle');
        if (!wrap) return false;
        wrap.style.display = '';
        return true;
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

        var compactEditType = document.getElementById('compactEditType');
        if (compactEditType) {
            compactEditType.addEventListener('change', function () {
                setCompactEditError('');
                setCompactIdentifierMode(compactEditType.value, editContext ? {
                    identifier: editContext.original.identifier,
                    displayIdentifier: editContext.original.displayIdentifier,
                    name: editContext.original.displayIdentifier,
                } : null);
            });
        }

        var compactSavedSearch = document.getElementById('compactEditSavedSearch');
        if (compactSavedSearch) {
            compactSavedSearch.addEventListener('input', function () {
                var hidden = document.getElementById('compactEditSavedId');
                if (hidden) hidden.value = '';
                setCompactEditError('');
                renderCompactSavedDropdown();
            });
            compactSavedSearch.addEventListener('focus', renderCompactSavedDropdown);
            compactSavedSearch.addEventListener('change', inferSavedSelectionFromSearch);
        }

        var compactFile = document.getElementById('compactEditFile');
        if (compactFile) {
            compactFile.addEventListener('change', function () {
                setCompactEditError('');
                updateCompactFileName();
            });
        }

        document.addEventListener('click', function (e) {
            var savedWrap = document.getElementById('compactEditSavedContainer');
            if (savedWrap && !savedWrap.contains(e.target)) {
                hideCompactSavedDropdown();
            }
        });

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
