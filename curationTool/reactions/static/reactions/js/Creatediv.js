var button_to_div = {
    'reactants-button': 'reactants',
    'atommapping-button': 'atommapping',
    'cheminfo-button': 'cheminfo',
    'metinfo-button': 'metaboliteinfo',
    'references-button': 'references',
    'extlinks-button': 'extlinks',
    'reactinfo-button': 'reactioninfo',
    'geneinfo-button': 'gene_info',
    'comments-button': 'comments',
    'reactiontemps-button': 'reaction_temps'
    // 'viewsavedreactions-button' is excluded
};

(function () {
    const PANEL_ORDER_KEY = 'panelOrder';
    const PANEL_SIZE_KEY = 'panelSizes';

    const WorkspacePanels = {
        init() {
            this.container = document.getElementById('workspacePanels') || document.querySelector('.divs-container');
            this.buttons = document.querySelectorAll('.dynamic-button-side-button');
            if (!this.container) {
                return;
            }

            this.enhancePanels();
            this.applyStoredOrder();
            this.applyStoredSizes();

            restoreVisibleDivs();
            this.syncInitialVisibility();
            refreshSideButtons();

            this.bindButtons();
            this.initializeDragAndDrop();
            this.setupAboutLink();
        },

        panelList() {
            return Array.from(this.container.querySelectorAll('.content-div'));
        },

        bindButtons() {
            this.buttons.forEach((button) => {
                if (!(button.id in button_to_div)) return;
                button.addEventListener('click', () => this.togglePanel(button));
            });
        },

        togglePanel(button) {
            const targetName = button_to_div[button.id];
            const panel = document.getElementsByName(targetName + '-div')[0];
            if (!panel) return;

            const isVisible = panel.style.display === 'block';
            if (isVisible) {
                this.hidePanel(panel, button);
            } else {
                this.showPanel(panel, button);
            }
            updateLocalStorageState();
        },

        showPanel(panel, button) {
            panel.style.display = 'block';
            panel.classList.add('is-visible');
            if (button) {
                button.classList.add('active');
            }
            this.activatePanel(panel);
        },

        hidePanel(panel, button) {
            panel.style.display = 'none';
            panel.classList.remove('is-visible');
            if (button) {
                button.classList.remove('active');
            }
        },

        activatePanel(panel) {
            this.container.appendChild(panel);
            panel.classList.add('just-activated');
            setTimeout(() => panel.classList.remove('just-activated'), 450);
            this.persistPanelOrder();
        },

        enhancePanels() {
            this.panelList().forEach((panel) => {
                const header = panel.querySelector('.div-header');
                if (header) {
                    header.classList.add('panel-drag-handle');
                    const actions = this.ensureHeaderActions(panel, header);
                    if (actions && !actions.querySelector('.panel-handle-icon')) {
                        const icon = document.createElement('span');
                        icon.className = 'panel-control panel-control--drag panel-handle-icon';
                        icon.setAttribute('title', 'Drag to reorder');
                        icon.setAttribute('aria-label', 'Drag to reorder');
                        icon.innerHTML = '<i class="fas fa-grip-vertical" aria-hidden="true"></i>';
                        actions.appendChild(icon);
                    }
                }
                this.addResetButton(panel, header);
                this.addResizeHandle(panel);
            });
        },

        ensureHeaderActions(panel, header) {
            if (!header) return null;
            let actions = header.querySelector('.panel-header-actions');
            if (!actions) {
                actions = document.createElement('span');
                actions.className = 'panel-header-actions';
                actions.setAttribute('role', 'group');
                actions.setAttribute('aria-label', 'Panel controls');
                header.appendChild(actions);
            } else if (!actions.getAttribute('role')) {
                actions.setAttribute('role', 'group');
                actions.setAttribute('aria-label', 'Panel controls');
            }
            return actions;
        },

        addResetButton(panel, header) {
            if (!header) return;
            const actions = this.ensureHeaderActions(panel, header);
            if (!actions || actions.querySelector('.panel-reset-button')) {
                return;
            }
            const resetBtn = document.createElement('button');
            resetBtn.type = 'button';
            resetBtn.className = 'panel-control panel-reset-button';
            resetBtn.title = 'Reset panel size';
            resetBtn.setAttribute('aria-label', 'Reset panel size');
            resetBtn.innerHTML = '<i class="fas fa-redo-alt" aria-hidden="true"></i>';
            resetBtn.addEventListener('click', (event) => {
                event.preventDefault();
                event.stopPropagation();
                this.resetPanelSize(panel);
            });
            actions.appendChild(resetBtn);
        },

        addResizeHandle(panel) {
            if (!panel || panel.querySelector('.panel-resize-handle')) {
                return;
            }
            const handle = document.createElement('button');
            handle.type = 'button';
            handle.className = 'panel-resize-handle';
            handle.title = 'Resize panel';
            handle.setAttribute('aria-label', 'Resize panel');
            handle.innerHTML = '<i class="fas fa-expand-alt" aria-hidden="true"></i>';
            panel.appendChild(handle);

            let startWidth = 0;
            let startHeight = 0;
            let startX = 0;
            let startY = 0;

            const onPointerMove = (event) => {
                const deltaX = event.clientX - startX;
                const deltaY = event.clientY - startY;
                const newWidth = Math.max(320, startWidth + deltaX);
                const newHeight = Math.max(280, startHeight + deltaY);
                panel.style.flex = '0 0 auto';
                panel.style.width = newWidth + 'px';
                panel.style.height = newHeight + 'px';
            };

            const onPointerUp = () => {
                panel.classList.remove('is-resizing');
                document.removeEventListener('pointermove', onPointerMove);
                document.removeEventListener('pointerup', onPointerUp);
                this.persistPanelSizes();
            };

            handle.addEventListener('pointerdown', (event) => {
                event.preventDefault();
                event.stopPropagation();
                startWidth = panel.offsetWidth;
                startHeight = panel.offsetHeight;
                startX = event.clientX;
                startY = event.clientY;
                panel.classList.add('is-resizing');

                document.addEventListener('pointermove', onPointerMove);
                document.addEventListener('pointerup', onPointerUp);
            });
        },

        initializeDragAndDrop() {
            if (!window.Sortable || !this.container) {
                console.warn('Sortable.js is unavailable. Drag-and-drop ordering disabled.');
                return;
            }
            Sortable.create(this.container, {
                animation: 180,
                handle: '.panel-drag-handle',
                draggable: '.content-div',
                ghostClass: 'workspace-panel-ghost',
                filter: '.panel-resize-handle, .panel-reset-button, .panel-reset-button *',
                preventOnFilter: false,
                onStart: (evt) => evt.item.classList.add('is-dragging'),
                onEnd: (evt) => {
                    evt.item.classList.remove('is-dragging');
                    this.persistPanelOrder();
                },
            });
        },

        persistPanelOrder() {
            if (!this.container) return;
            const order = this.panelList().map((panel) => panel.getAttribute('name'));
            localStorage.setItem(PANEL_ORDER_KEY, JSON.stringify(order));
        },

        applyStoredOrder() {
            const stored = localStorage.getItem(PANEL_ORDER_KEY);
            if (!stored) return;
            try {
                const order = JSON.parse(stored);
                order.forEach((name) => {
                    const panel = this.container.querySelector(`.content-div[name="${name}"]`);
                    if (panel) {
                        this.container.appendChild(panel);
                    }
                });
            } catch (error) {
                console.warn('Unable to restore panel order', error);
            }
        },

        persistPanelSizes() {
            const sizes = {};
            this.panelList().forEach((panel) => {
                if (panel.style.width || panel.style.height) {
                    sizes[panel.getAttribute('name')] = {
                        width: panel.style.width || panel.offsetWidth + 'px',
                        height: panel.style.height || panel.offsetHeight + 'px',
                    };
                }
            });
            localStorage.setItem(PANEL_SIZE_KEY, JSON.stringify(sizes));
        },

        applyStoredSizes() {
            const saved = localStorage.getItem(PANEL_SIZE_KEY);
            if (!saved) return;
            try {
                const sizes = JSON.parse(saved);
                this.panelList().forEach((panel) => {
                    const key = panel.getAttribute('name');
                    if (sizes[key]) {
                        const { width, height } = sizes[key];
                        if (width) {
                            panel.style.flex = '0 0 auto';
                            panel.style.width = width;
                        }
                        if (height) {
                            panel.style.height = height;
                        }
                    }
                });
            } catch (error) {
                console.warn('Unable to restore panel sizes', error);
            }
        },

        resetPanelSize(panel) {
            if (!panel) return;
            panel.classList.remove('is-resizing');
            panel.style.removeProperty('width');
            panel.style.removeProperty('height');
            panel.style.removeProperty('max-width');
            panel.style.removeProperty('flex');
            panel.style.removeProperty('flex-basis');
            panel.style.removeProperty('flex-grow');
            panel.style.removeProperty('flex-shrink');
            this.persistPanelSizes();
        },

        syncInitialVisibility() {
            this.panelList().forEach((panel) => {
                if (panel.style.display === 'block') {
                    panel.classList.add('is-visible');
                }
            });
        },

        setupAboutLink() {
            const aboutButton = document.getElementById('item about-item');
            if (!aboutButton) return;
            aboutButton.addEventListener('click', function () {
                if (typeof aboutUrl === 'string') {
                    window.open(aboutUrl, '_blank', 'noopener');
                }
            });
        },
    };

    document.addEventListener('DOMContentLoaded', function () {
        WorkspacePanels.init();
    });
})();

function refreshSideButtons() {
    for (const [button, div] of Object.entries(button_to_div)) {
        const buttonElement = document.getElementById(button);
        const divElement = document.getElementsByName(div + '-div')[0];
        if (buttonElement && divElement) {
            const isVisible = divElement.style.display === 'block' || divElement.classList.contains('is-visible');
            if (isVisible) {
                buttonElement.classList.add('active');
            } else {
                buttonElement.classList.remove('active');
            }
        }
    }
}