/**
 * Graph Visualization Module
 * 
 * Professional hypergraph visualization for metabolic reaction networks.
 * Uses Cytoscape.js for rendering with custom styling and interactions.
 * 
 * Features:
 * - Metabolites as circular nodes (VMH: blue, Saved: purple)
 * - Reactions as diamond nodes connecting substrates to products
 * - Interactive tooltips with detailed information
 * - Click-to-edit reactions, click-to-view VMH metabolites
 * - Search/filter functionality
 * - Pan, zoom, and fit controls
 * - Responsive design
 */

(function() {
    'use strict';
    const vmhBaseUrl = (window.vmhBaseUrl || '').replace(/\/$/, '');

    // Configuration
    const CONFIG = {
        api: {
            graphInfo: '/get_graph_info/'
        },
        urls: {
            editReaction: (id) => `/?reaction_id=${id}&action=edit`,
            vmhMetabolite: (abbr) => vmhBaseUrl ? `${vmhBaseUrl}/metabolite/${encodeURIComponent(abbr)}` : null
        },
        colors: {
            vmhNode: '#4fc3f7',
            savedNode: '#ba68c8',
            reactionBalanced: '#81c784',
            reactionUnbalanced: '#ff8a65',
            edge: '#4a5568',
            edgeHighlight: '#667eea',
            background: '#1a1a2e'
        },
        // Default layout - uses cose-bilkent for better metabolic network visualization
        // This is overridden by the layout control panel
        defaultLayout: {
            name: 'cose-bilkent',
            quality: 'proof',
            animate: true,
            animationDuration: 1000,
            animationEasing: 'ease-out-cubic',
            fit: true,
            padding: 60,
            nodeDimensionsIncludeLabels: true,
            // Tuned for metabolic networks - spread out, not jumbled
            nodeRepulsion: 6000,
            idealEdgeLength: 150,
            edgeElasticity: 0.45,
            // Gentle gravity - keeps nodes from flying too far
            gravity: 0.15,
            gravityRange: 3.8,
            gravityCompound: 1.0,
            // Spacing
            tile: true,
            tilingPaddingVertical: 40,
            tilingPaddingHorizontal: 40,
            // Quality settings
            numIter: 2500,
            nestingFactor: 0.1,
            randomize: false
        }
    };

    // State
    let cy = null;
    let graphData = null;
    let tooltip = null;
    let contextMenu = null;
    let selectedElement = null;
    let isInfoPanelVisible = false;

    /**
     * Initialize the graph visualization module
     */
    function init() {
        createModalHTML();
        bindGlobalEvents();
    }

    /**
     * Create the modal HTML structure
     */
    function createModalHTML() {
        const modalHTML = `
            <div class="graph-modal-overlay" id="graphModalOverlay">
                <div class="graph-modal">
                    <!-- Header -->
                    <div class="graph-header">
                        <div class="graph-header-left">
                            <h2 class="graph-title">
                                <i class="fa fa-project-diagram"></i>
                                Reaction Network
                                <span class="graph-filter-badge" id="graphFilterBadge" style="display: none;">
                                    <i class="fa fa-filter"></i> Selected Only
                                </span>
                            </h2>
                            <div class="graph-stats" id="graphStats">
                                <div class="graph-stat">
                                    <span class="graph-stat-icon vmh"></span>
                                    <span class="graph-stat-value" id="statVmh">0</span>
                                    <span>VMH</span>
                                </div>
                                <div class="graph-stat">
                                    <span class="graph-stat-icon saved"></span>
                                    <span class="graph-stat-value" id="statSaved">0</span>
                                    <span>New</span>
                                </div>
                                <div class="graph-stat">
                                    <span class="graph-stat-icon reaction"></span>
                                    <span class="graph-stat-value" id="statReactions">0</span>
                                    <span>Reactions</span>
                                </div>
                            </div>
                        </div>
                        <div class="graph-header-controls">
                            <div class="graph-search-container">
                                <span class="graph-search-icon-wrapper">
                                    <i class="fa fa-search"></i>
                                </span>
                                <input type="text" class="graph-search" id="graphSearch" 
                                       placeholder="Search metabolites... (/)">
                            </div>
                            <button class="graph-btn" id="btnFit" title="Fit to screen (F)">
                                <i class="fa fa-expand"></i>
                            </button>
                            <button class="graph-btn" id="btnRelayout" title="Re-layout (R)">
                                <i class="fa fa-sync-alt"></i>
                            </button>
                            <button class="graph-btn" id="btnToggleLabels" title="Toggle labels (L)">
                                <i class="fa fa-font"></i>
                            </button>
                            <button class="graph-btn" id="btnHelp" title="Keyboard shortcuts">
                                <i class="fa fa-keyboard"></i>
                            </button>
                            <button class="graph-btn-close" id="btnCloseGraph" title="Close (Esc)">
                                <i class="fa fa-times"></i>
                            </button>
                        </div>
                    </div>

                    <!-- Main Graph Container -->
                    <div class="graph-container">
                        <div id="cy"></div>

                        <!-- Sidebar Controls -->
                        <div class="graph-sidebar">
                            <button class="graph-sidebar-btn" id="btnZoomIn" title="Zoom in">
                                <i class="fa fa-plus"></i>
                            </button>
                            <button class="graph-sidebar-btn" id="btnZoomOut" title="Zoom out">
                                <i class="fa fa-minus"></i>
                            </button>
                            <button class="graph-sidebar-btn" id="btnCenter" title="Center graph">
                                <i class="fa fa-crosshairs"></i>
                            </button>
                            <button class="graph-sidebar-btn" id="btnScreenshot" title="Download image">
                                <i class="fa fa-camera"></i>
                            </button>
                            <div class="sidebar-divider"></div>
                            <button class="graph-sidebar-btn" id="btnLayoutPanel" title="Layout controls">
                                <i class="fa fa-sliders-h"></i>
                            </button>
                        </div>

                        <!-- Layout Control Panel -->
                        <div class="layout-control-panel" id="layoutControlPanel">
                            <div class="layout-panel-header">
                                <span class="layout-panel-title">
                                    <i class="fa fa-project-diagram"></i> Layout
                                </span>
                                <button class="layout-panel-close" id="btnCloseLayoutPanel">
                                    <i class="fa fa-times"></i>
                                </button>
                            </div>
                            
                            <div class="layout-panel-section">
                                <label class="layout-label">Algorithm</label>
                                <select class="layout-select" id="layoutAlgorithm">
                                    <option value="cose-bilkent">Force-Directed (COSE)</option>
                                    <option value="cola">Cola (Constraint-Based)</option>
                                    <option value="dagre">Dagre (Hierarchical)</option>
                                    <option value="concentric">Concentric (by Degree)</option>
                                    <option value="circle">Circle (Reactions Center)</option>
                                    <option value="grid">Grid</option>
                                </select>
                            </div>
                            
                            <div class="layout-panel-section" id="spacingSection">
                                <label class="layout-label">
                                    <span>Node Spacing</span>
                                    <span class="layout-value" id="nodeSpacingValue">5</span>
                                </label>
                                <input type="range" class="layout-slider" id="nodeSpacing" 
                                       min="1" max="300" value="150" step="5">
                            </div>
                            
                            <div class="layout-panel-section" id="edgeLengthSection">
                                <label class="layout-label">
                                    <span>Edge Length</span>
                                    <span class="layout-value" id="edgeLengthValue">100</span>
                                </label>
                                <input type="range" class="layout-slider" id="edgeLength" 
                                       min="1" max="300" value="100" step="10">
                            </div>
                            
                            <div class="layout-panel-section physics-controls" id="physicsSection">
                                <label class="layout-label">
                                    <span>Gravity</span>
                                    <span class="layout-value" id="gravityValue">0.25</span>
                                </label>
                                <input type="range" class="layout-slider" id="gravitySlider" 
                                       min="0" max="2" value="0.25" step="0.05">
                            </div>
                            
                            <div class="layout-panel-section physics-controls" id="repulsionSection">
                                <label class="layout-label">
                                    <span>Node Repulsion</span>
                                    <span class="layout-value" id="repulsionValue">4500</span>
                                </label>
                                <input type="range" class="layout-slider" id="repulsionSlider" 
                                       min="1000" max="20000" value="4500" step="500">
                            </div>
                            
                            <div class="layout-panel-divider"></div>
                            
                            <div class="layout-panel-section">
                                <label class="layout-label">Quick Presets</label>
                                <div class="layout-presets">
                                    <button class="layout-preset-btn" data-preset="compact" title="Compact layout">
                                        <i class="fa fa-compress-arrows-alt"></i>
                                        <span>Compact</span>
                                    </button>
                                    <button class="layout-preset-btn" data-preset="spread" title="Spread out layout">
                                        <i class="fa fa-expand-arrows-alt"></i>
                                        <span>Spread</span>
                                    </button>
                                    <button class="layout-preset-btn" data-preset="hierarchical" title="Hierarchical layout">
                                        <i class="fa fa-sitemap"></i>
                                        <span>Hierarchy</span>
                                    </button>
                                    <button class="layout-preset-btn" data-preset="circular" title="Circular layout">
                                        <i class="fa fa-circle-notch"></i>
                                        <span>Circular</span>
                                    </button>
                                </div>
                            </div>
                            
                            <div class="layout-panel-footer">
                                <button class="layout-apply-btn" id="btnApplyLayout">
                                    <i class="fa fa-play"></i> Apply Layout
                                </button>
                            </div>
                        </div>

                        <!-- Legend -->
                        <div class="graph-legend">
                            <div class="graph-legend-title">Legend</div>
                            <div class="graph-legend-items">
                                <div class="graph-legend-item">
                                    <div class="legend-node vmh"></div>
                                    <span>VMH Metabolite</span>
                                </div>
                                <div class="graph-legend-item">
                                    <div class="legend-node saved"></div>
                                    <span>New Metabolite</span>
                                </div>
                                <div class="graph-legend-item">
                                    <div class="legend-node reaction"></div>
                                    <span>Balanced Reaction</span>
                                </div>
                                <div class="graph-legend-item">
                                    <div class="legend-node reaction reaction-unbalanced"></div>
                                    <span>Unbalanced Reaction</span>
                                </div>
                            </div>
                        </div>

                        <!-- Info Panel -->
                        <div class="graph-info-panel" id="infoPanel">
                            <div class="info-panel-header">
                                <span class="info-panel-title" id="infoPanelTitle">Details</span>
                                <button class="info-panel-close" id="btnCloseInfo">
                                    <i class="fa fa-times"></i>
                                </button>
                            </div>
                            <div class="info-panel-content" id="infoPanelContent">
                            </div>
                        </div>

                        <!-- Loading Overlay -->
                        <div class="graph-loading" id="graphLoading">
                            <div class="graph-loading-spinner"></div>
                            <div class="graph-loading-text">Loading reaction network...</div>
                        </div>

                        <!-- Tooltip -->
                        <div class="graph-tooltip" id="graphTooltip"></div>

                        <!-- Context Menu -->
                        <div class="graph-context-menu" id="graphContextMenu"></div>
                        
                        <!-- Help Modal -->
                        <div class="graph-help-overlay" id="graphHelpOverlay">
                            <div class="graph-help-modal">
                                <div class="graph-help-header">
                                    <h3>Keyboard Shortcuts & Tips</h3>
                                    <button class="info-panel-close" id="btnCloseHelp">
                                        <i class="fa fa-times"></i>
                                    </button>
                                </div>
                                <div class="graph-help-content">
                                    <div class="help-section">
                                        <div class="help-section-title">Navigation</div>
                                        <div class="help-row"><kbd>F</kbd> Fit graph to screen</div>
                                        <div class="help-row"><kbd>R</kbd> Re-layout graph</div>
                                        <div class="help-row"><kbd>L</kbd> Toggle labels</div>
                                        <div class="help-row"><kbd>P</kbd> Toggle layout panel</div>
                                        <div class="help-row"><kbd>+</kbd> / <kbd>-</kbd> Zoom in/out</div>
                                        <div class="help-row"><kbd>/</kbd> Focus search</div>
                                        <div class="help-row"><kbd>Esc</kbd> Close panel/modal</div>
                                    </div>
                                    <div class="help-section">
                                        <div class="help-section-title">Mouse Interactions</div>
                                        <div class="help-row"><strong>Hover</strong> on node to see details</div>
                                        <div class="help-row"><strong>Click</strong> on node to open info panel</div>
                                        <div class="help-row"><strong>Double-click</strong> reaction to edit</div>
                                        <div class="help-row"><strong>Double-click</strong> VMH metabolite to view on VMH</div>
                                        <div class="help-row"><strong>Right-click</strong> for context menu</div>
                                        <div class="help-row"><strong>Drag</strong> to pan, <strong>Scroll</strong> to zoom</div>
                                    </div>
                                    <div class="help-section">
                                        <div class="help-section-title">Layout Controls</div>
                                        <div class="help-row">
                                            <i class="fa fa-sliders-h" style="color: #667eea;"></i>
                                            Click slider icon in sidebar to open layout panel
                                        </div>
                                        <div class="help-row">
                                            <strong>Algorithms:</strong> COSE (force), Cola, Dagre (hierarchy)
                                        </div>
                                        <div class="help-row">
                                            <strong>Presets:</strong> Compact, Spread, Hierarchy, Circular
                                        </div>
                                        <div class="help-row">
                                            <strong>Sliders:</strong> Adjust spacing, gravity, repulsion
                                        </div>
                                    </div>
                                    <div class="help-section">
                                        <div class="help-section-title">Visual Guide</div>
                                        <div class="help-row">
                                            <span class="legend-dot" style="background: #4fc3f7;"></span>
                                            VMH Metabolites (existing in database)
                                        </div>
                                        <div class="help-row">
                                            <span class="legend-dot" style="background: #ba68c8;"></span>
                                            New Metabolites (user-created)
                                        </div>
                                        <div class="help-row">
                                            <span class="legend-diamond" style="background: #81c784;"></span>
                                            Balanced Reactions
                                        </div>
                                        <div class="help-row">
                                            <span class="legend-diamond" style="background: #ff8a65;"></span>
                                            Unbalanced Reactions
                                        </div>
                                        <div class="help-row">
                                            <i class="fa fa-circle" style="font-size: 8px; color: #667eea;"></i>
                                            Node size = number of connections
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        document.body.insertAdjacentHTML('beforeend', modalHTML);
        
        tooltip = document.getElementById('graphTooltip');
        contextMenu = document.getElementById('graphContextMenu');
    }

    /**
     * Bind global event handlers
     */
    function bindGlobalEvents() {
        // Close button
        document.getElementById('btnCloseGraph').addEventListener('click', closeGraph);

        // Keyboard shortcuts
        document.addEventListener('keydown', (e) => {
            // Only handle if graph modal is active
            if (!document.getElementById('graphModalOverlay').classList.contains('active')) {
                return;
            }
            
            // Escape key to close context menu or modal
            if (e.key === 'Escape') {
                const helpOverlay = document.getElementById('graphHelpOverlay');
                if (helpOverlay.classList.contains('visible')) {
                    hideHelp();
                } else if (contextMenu.classList.contains('visible')) {
                    hideContextMenu();
                } else if (isInfoPanelVisible) {
                    hideInfoPanel();
                } else {
                    closeGraph();
                }
                return;
            }
            
            // Don't trigger shortcuts when typing in search
            if (document.activeElement === document.getElementById('graphSearch')) {
                return;
            }
            
            // Keyboard shortcuts for graph controls
            switch (e.key.toLowerCase()) {
                case 'f':
                    // Fit to screen
                    if (cy) cy.fit(50);
                    break;
                case 'r':
                    // Re-layout
                    relayout();
                    break;
                case 'l':
                    // Toggle labels
                    toggleLabels();
                    break;
                case 'p':
                    // Toggle layout panel
                    toggleLayoutPanel();
                    break;
                case '+':
                case '=':
                    // Zoom in
                    if (cy) cy.zoom(cy.zoom() * 1.3);
                    break;
                case '-':
                    // Zoom out
                    if (cy) cy.zoom(cy.zoom() / 1.3);
                    break;
                case '/':
                    // Focus search
                    e.preventDefault();
                    document.getElementById('graphSearch').focus();
                    break;
                case '?':
                    // Show help
                    toggleHelp();
                    break;
                case 'h':
                    // Show help (alternative)
                    toggleHelp();
                    break;
            }
        });

        // Zoom controls
        document.getElementById('btnZoomIn').addEventListener('click', () => {
            if (cy) cy.zoom(cy.zoom() * 1.3);
        });
        document.getElementById('btnZoomOut').addEventListener('click', () => {
            if (cy) cy.zoom(cy.zoom() / 1.3);
        });
        document.getElementById('btnCenter').addEventListener('click', () => {
            if (cy) cy.center();
        });
        document.getElementById('btnFit').addEventListener('click', () => {
            if (cy) cy.fit(50);
        });
        document.getElementById('btnRelayout').addEventListener('click', relayout);
        document.getElementById('btnToggleLabels').addEventListener('click', toggleLabels);
        document.getElementById('btnScreenshot').addEventListener('click', downloadScreenshot);
        document.getElementById('btnHelp').addEventListener('click', toggleHelp);
        document.getElementById('btnCloseHelp').addEventListener('click', hideHelp);
        
        // Close help on overlay click
        document.getElementById('graphHelpOverlay').addEventListener('click', (e) => {
            if (e.target.id === 'graphHelpOverlay') {
                hideHelp();
            }
        });

        // Search
        const searchInput = document.getElementById('graphSearch');
        searchInput.addEventListener('input', debounce(handleSearch, 200));

        // Close info panel
        document.getElementById('btnCloseInfo').addEventListener('click', hideInfoPanel);

        // Click outside to close context menu
        document.addEventListener('click', (e) => {
            if (!contextMenu.contains(e.target)) {
                hideContextMenu();
            }
        });
        
        // Layout Control Panel Events
        bindLayoutControlEvents();
    }
    
    // =========================================
    // LAYOUT CONTROL SYSTEM
    // =========================================
    
    // Current layout state
    let currentLayoutConfig = {
        algorithm: 'cose-bilkent',
        nodeSpacing: 150,
        edgeLength: 100,
        gravity: 0.25,
        repulsion: 4500
    };
    
    /**
     * Bind layout control panel events
     */
    function bindLayoutControlEvents() {
        // Toggle layout panel
        document.getElementById('btnLayoutPanel').addEventListener('click', toggleLayoutPanel);
        document.getElementById('btnCloseLayoutPanel').addEventListener('click', hideLayoutPanel);
        
        // Algorithm selector
        document.getElementById('layoutAlgorithm').addEventListener('change', (e) => {
            currentLayoutConfig.algorithm = e.target.value;
            updatePhysicsControlsVisibility(e.target.value);
        });
        
        // Sliders with live value updates
        const sliders = [
            { id: 'nodeSpacing', prop: 'nodeSpacing', display: 'nodeSpacingValue' },
            { id: 'edgeLength', prop: 'edgeLength', display: 'edgeLengthValue' },
            { id: 'gravitySlider', prop: 'gravity', display: 'gravityValue' },
            { id: 'repulsionSlider', prop: 'repulsion', display: 'repulsionValue' }
        ];
        
        sliders.forEach(slider => {
            const element = document.getElementById(slider.id);
            element.addEventListener('input', (e) => {
                const value = parseFloat(e.target.value);
                currentLayoutConfig[slider.prop] = value;
                document.getElementById(slider.display).textContent = 
                    slider.prop === 'gravity' ? value.toFixed(2) : value;
            });
        });
        
        // Preset buttons
        document.querySelectorAll('.layout-preset-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                applyLayoutPreset(btn.dataset.preset);
                // Update active state
                document.querySelectorAll('.layout-preset-btn').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
            });
        });
        
        // Apply layout button
        document.getElementById('btnApplyLayout').addEventListener('click', applyCurrentLayout);
    }
    
    /**
     * Toggle layout panel visibility
     */
    function toggleLayoutPanel() {
        const panel = document.getElementById('layoutControlPanel');
        const btn = document.getElementById('btnLayoutPanel');
        
        if (panel.classList.contains('visible')) {
            hideLayoutPanel();
        } else {
            panel.classList.add('visible');
            btn.classList.add('active');
        }
    }
    
    /**
     * Hide layout panel
     */
    function hideLayoutPanel() {
        document.getElementById('layoutControlPanel').classList.remove('visible');
        document.getElementById('btnLayoutPanel').classList.remove('active');
    }
    
    /**
     * Update physics controls visibility based on algorithm
     */
    function updatePhysicsControlsVisibility(algorithm) {
        const physicsControls = document.querySelectorAll('.physics-controls');
        const spacingSection = document.getElementById('spacingSection');
        const edgeLengthSection = document.getElementById('edgeLengthSection');
        
        // Force-directed algorithms show physics controls
        const isPhysicsBased = ['cose-bilkent', 'cola'].includes(algorithm);
        physicsControls.forEach(ctrl => {
            ctrl.classList.toggle('hidden', !isPhysicsBased);
        });
        
        // Adjust spacing label based on algorithm
        if (['grid', 'circle', 'concentric'].includes(algorithm)) {
            spacingSection.querySelector('span:first-child').textContent = 'Spacing';
        } else {
            spacingSection.querySelector('span:first-child').textContent = 'Node Spacing';
        }
    }
    
    /**
     * Apply a layout preset
     */
    function applyLayoutPreset(preset) {
        const presets = {
            compact: {
                algorithm: 'cose-bilkent',
                nodeSpacing: 25,
                edgeLength: 60,
                gravity: 0.8,
                repulsion: 2000
            },
            spread: {
                algorithm: 'cose-bilkent',
                nodeSpacing: 100,
                edgeLength: 200,
                gravity: 0.1,
                repulsion: 10000
            },
            hierarchical: {
                algorithm: 'dagre',
                nodeSpacing: 60,
                edgeLength: 100,
                gravity: 0.25,
                repulsion: 4500
            },
            circular: {
                algorithm: 'concentric',
                nodeSpacing: 80,
                edgeLength: 120,
                gravity: 0.25,
                repulsion: 4500
            }
        };
        
        const config = presets[preset];
        if (!config) return;
        
        // Update state
        currentLayoutConfig = { ...config };
        
        // Update UI controls to match preset
        document.getElementById('layoutAlgorithm').value = config.algorithm;
        document.getElementById('nodeSpacing').value = config.nodeSpacing;
        document.getElementById('nodeSpacingValue').textContent = config.nodeSpacing;
        document.getElementById('edgeLength').value = config.edgeLength;
        document.getElementById('edgeLengthValue').textContent = config.edgeLength;
        document.getElementById('gravitySlider').value = config.gravity;
        document.getElementById('gravityValue').textContent = config.gravity.toFixed(2);
        document.getElementById('repulsionSlider').value = config.repulsion;
        document.getElementById('repulsionValue').textContent = config.repulsion;
        
        // Update physics controls visibility
        updatePhysicsControlsVisibility(config.algorithm);
        
        // Auto-apply the layout
        applyCurrentLayout();
    }
    
    /**
     * Build layout options for current configuration
     */
    function buildLayoutOptions() {
        const { algorithm, nodeSpacing, edgeLength, gravity, repulsion } = currentLayoutConfig;
        
        const baseOptions = {
            fit: true,
            padding: 50,
            animate: true,
            animationDuration: 800,
            animationEasing: 'ease-out'
        };
        
        switch (algorithm) {
            case 'cose-bilkent':
                return {
                    name: 'cose-bilkent',
                    ...baseOptions,
                    quality: 'proof',
                    nodeDimensionsIncludeLabels: true,
                    // Node spacing
                    nodeRepulsion: repulsion,
                    idealEdgeLength: edgeLength,
                    edgeElasticity: 0.45,
                    // Physics
                    gravity: gravity,
                    gravityRange: 3.8,
                    gravityCompound: 1.0,
                    gravityRangeCompound: 1.5,
                    // Layout quality
                    numIter: 2500,
                    tile: true,
                    tilingPaddingVertical: nodeSpacing,
                    tilingPaddingHorizontal: nodeSpacing,
                    // Nesting factor for compound nodes
                    nestingFactor: 0.1,
                    // Randomize initial positions for variety
                    randomize: false,
                    // Animation
                    animationEasing: 'ease-out-cubic'
                };
                
            case 'cola':
                return {
                    name: 'cola',
                    ...baseOptions,
                    maxSimulationTime: 4000,
                    ungrabifyWhileSimulating: false,
                    // Node spacing
                    nodeSpacing: function() { return nodeSpacing; },
                    edgeLength: function() { return edgeLength; },
                    // Forces
                    edgeSymDiffLength: edgeLength * 0.3,
                    avoidOverlap: true,
                    convergenceThreshold: 0.01,
                    // Flow (optional - for hierarchical tendencies)
                    flow: gravity > 0.5 ? { axis: 'y', minSeparation: 30 } : undefined
                };
                
            case 'dagre':
                return {
                    name: 'dagre',
                    ...baseOptions,
                    // Direction
                    rankDir: 'TB',
                    // Spacing
                    nodeSep: nodeSpacing,
                    edgeSep: Math.max(10, edgeLength * 0.3),
                    rankSep: edgeLength,
                    // Alignment
                    align: 'UL',
                    // Acyclic algorithm
                    acyclicer: 'greedy',
                    // Rank assignment
                    ranker: 'network-simplex'
                };
                
            case 'concentric':
                // Radial layout - high-degree nodes in center, but also group by type
                return {
                    name: 'concentric',
                    ...baseOptions,
                    // Combine degree with type preference (reactions slightly more central)
                    concentric: function(node) {
                        const degree = node.degree();
                        const type = node.data('type');
                        // Boost reactions slightly toward center
                        const typeBonus = type === 'reaction' ? 5 : 0;
                        return degree + typeBonus;
                    },
                    levelWidth: function(nodes) {
                        return Math.max(2, Math.ceil(nodes.length / 8));
                    },
                    minNodeSpacing: nodeSpacing,
                    spacingFactor: 1 + (edgeLength / 100),
                    equidistant: false,
                    startAngle: 3 / 2 * Math.PI,
                    sweep: 2 * Math.PI,
                    clockwise: true
                };
                
            case 'circle':
                // Use concentric layout with reactions in center, metabolites on perimeter
                return {
                    name: 'concentric',
                    ...baseOptions,
                    // Place reactions in center ring, metabolites on outer rings
                    concentric: function(node) {
                        const type = node.data('type');
                        if (type === 'reaction') {
                            return 100; // Highest value = center
                        } else if (type === 'saved') {
                            return 50;  // Middle ring for new metabolites
                        } else {
                            return 10;  // Outer ring for VMH metabolites
                        }
                    },
                    levelWidth: function(nodes) {
                        // Create distinct rings for each type
                        return 1;
                    },
                    minNodeSpacing: nodeSpacing,
                    spacingFactor: 1 + (edgeLength / 150),
                    equidistant: true,
                    startAngle: 3 / 2 * Math.PI,
                    sweep: 2 * Math.PI,
                    clockwise: true
                };
                
            case 'grid':
                return {
                    name: 'grid',
                    ...baseOptions,
                    spacingFactor: 1 + (nodeSpacing / 50),
                    avoidOverlap: true,
                    avoidOverlapPadding: nodeSpacing / 2,
                    condense: nodeSpacing < 50,
                    rows: undefined,
                    cols: undefined,
                    sort: function(a, b) {
                        // Group by type
                        if (a.data('type') !== b.data('type')) {
                            const order = ['reaction', 'vmh', 'saved'];
                            return order.indexOf(a.data('type')) - order.indexOf(b.data('type'));
                        }
                        return 0;
                    }
                };
                
            default:
                // Fallback to basic COSE
                return {
                    name: 'cose',
                    ...baseOptions,
                    nodeRepulsion: repulsion,
                    idealEdgeLength: edgeLength,
                    gravity: gravity,
                    numIter: 1500
                };
        }
    }
    
    /**
     * Apply the current layout configuration
     */
    function applyCurrentLayout() {
        if (!cy) return;
        
        const btn = document.getElementById('btnApplyLayout');
        const originalContent = btn.innerHTML;
        
        // Show running state
        btn.innerHTML = '<i class="fa fa-spinner"></i> Running...';
        btn.classList.add('running');
        
        // Build and run layout
        const layoutOptions = buildLayoutOptions();
        
        // Add completion callback
        layoutOptions.stop = function() {
            btn.innerHTML = originalContent;
            btn.classList.remove('running');
        };
        
        try {
            const layout = cy.layout(layoutOptions);
            layout.run();
        } catch (error) {
            console.error('Layout error:', error);
            btn.innerHTML = originalContent;
            btn.classList.remove('running');
            
            // Fallback to basic cose if layout fails
            const fallbackLayout = cy.layout({
                name: 'cose',
                animate: true,
                fit: true,
                padding: 50
            });
            fallbackLayout.run();
        }
    }
    
    /**
     * Toggle help modal
     */
    function toggleHelp() {
        const overlay = document.getElementById('graphHelpOverlay');
        overlay.classList.toggle('visible');
    }
    
    /**
     * Hide help modal
     */
    function hideHelp() {
        document.getElementById('graphHelpOverlay').classList.remove('visible');
    }

    /**
     * Open the graph visualization
     */
    async function openGraph() {
        const overlay = document.getElementById('graphModalOverlay');
        overlay.classList.add('active');
        document.body.style.overflow = 'hidden';

        showLoading(true);

        try {
            // Get CSRF token and user ID from global scope
            const csrfToken = window.csrfToken || document.querySelector('[name=csrfmiddlewaretoken]')?.value || '';
            const userID = window.userID || '';
            
            // Get checked/selected reaction IDs (if any are selected)
            const checkedReactions = window.checkedReactions || [];
            
            // Build request body
            let bodyParams = 'userID=' + encodeURIComponent(userID);
            
            // If reactions are checked, only show those in the graph
            if (checkedReactions.length > 0) {
                bodyParams += '&reactionIDs=' + encodeURIComponent(JSON.stringify(checkedReactions));
            }
            
            const response = await fetch(CONFIG.api.graphInfo, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/x-www-form-urlencoded',
                    'X-CSRFToken': csrfToken
                },
                body: bodyParams
            });
            
            if (!response.ok) throw new Error('Failed to fetch graph data');
            
            const result = await response.json();
            
            if (result.status !== 'success') {
                throw new Error(result.message || 'Unknown error');
            }
            
            graphData = result.data;
            
            // Check for empty graph
            if (graphData.stats.edge_count === 0) {
                showEmptyState();
                return;
            }
            
            // Show filter badge if viewing selected reactions only
            const filterBadge = document.getElementById('graphFilterBadge');
            if (checkedReactions.length > 0) {
                filterBadge.style.display = 'inline-flex';
            } else {
                filterBadge.style.display = 'none';
            }
            
            updateStats(graphData.stats);
            initCytoscape(graphData);
            
        } catch (error) {
            console.error('Error loading graph:', error);
            showError('Failed to load reaction network. Please try again.');
        } finally {
            showLoading(false);
        }
    }

    /**
     * Close the graph visualization
     */
    function closeGraph() {
        const overlay = document.getElementById('graphModalOverlay');
        overlay.classList.remove('active');
        document.body.style.overflow = '';
        
        hideTooltip();
        hideContextMenu();
        hideInfoPanel();
        
        // Hide filter badge
        document.getElementById('graphFilterBadge').style.display = 'none';
        
        if (cy) {
            cy.destroy();
            cy = null;
        }
    }

    /**
     * Initialize Cytoscape.js
     */
    function initCytoscape(data) {
        const elements = buildElements(data);

        cy = cytoscape({
            container: document.getElementById('cy'),
            elements: elements,
            style: getCytoscapeStyle(),
            layout: { name: 'preset' }, // We'll run layout after
            minZoom: 0.001,  // Allow zooming out much further for large networks
            maxZoom: 5,
            wheelSensitivity: 0.3
        });

        bindCytoscapeEvents();
        runLayout();
    }

    /**
     * Build Cytoscape elements from graph data
     */
    function buildElements(data) {
        const elements = [];
        const nodes = data.nodes || {};
        const edges = data.edges || [];

        // Add metabolite nodes
        Object.entries(nodes).forEach(([id, node]) => {
            // Determine label: prefer abbr, then name
            const label = node.abbr || node.name || id;
            
            elements.push({
                data: {
                    id: id,
                    label: label,
                    type: 'metabolite',
                    subtype: node.type, // 'vmh' or 'saved'
                    formula: node.formula,
                    name: node.name,
                    compartment: node.compartment,
                    source: node.source_type,
                    originalId: node.saved_id,
                    inVmh: node.in_vmh,
                    vmhUrl: node.vmh_url,
                    inchiKey: node.inchi_key,
                    degree: 0 // Will be calculated
                }
            });
        });

        // Add reaction nodes and edges
        edges.forEach((edge, index) => {
            const reactionId = `reaction_${edge.reaction_id}`;
            
            // Determine if balanced (can be true, false, or null)
            const isBalanced = edge.balanced === true;
            
            // Add reaction node
            elements.push({
                data: {
                    id: reactionId,
                    label: edge.name || `R${edge.reaction_id}`,
                    type: 'reaction',
                    reactionDbId: edge.reaction_id,
                    abbreviation: edge.name,
                    description: edge.description,
                    balanced: isBalanced,
                    balancedRaw: edge.balanced, // Keep original for display (could be null)
                    confidenceScore: edge.confidence_score,
                    substrateCount: edge.substrates.length,
                    productCount: edge.products.length,
                    subsystem: edge.subsystem,
                    direction: edge.direction,
                    reversible: edge.reversible,
                    vmhFound: edge.vmh_found,
                    flags: edge.flags || []
                }
            });

            // Add edges from substrates to reaction
            edge.substrates.forEach((subId) => {
                elements.push({
                    data: {
                        id: `${subId}_to_${reactionId}`,
                        source: subId,
                        target: reactionId,
                        type: 'substrate'
                    }
                });
                
                // Increment degree
                const nodeData = elements.find(e => e.data.id === subId);
                if (nodeData) nodeData.data.degree = (nodeData.data.degree || 0) + 1;
            });

            // Add edges from reaction to products
            edge.products.forEach((prodId) => {
                elements.push({
                    data: {
                        id: `${reactionId}_to_${prodId}`,
                        source: reactionId,
                        target: prodId,
                        type: 'product'
                    }
                });
                
                // Increment degree
                const nodeData = elements.find(e => e.data.id === prodId);
                if (nodeData) nodeData.data.degree = (nodeData.data.degree || 0) + 1;
            });
        });

        return elements;
    }

    /**
     * Get Cytoscape style configuration
     */
    function getCytoscapeStyle() {
        return [
            // Metabolite nodes - base style
            {
                selector: 'node[type="metabolite"]',
                style: {
                    'shape': 'ellipse',
                    'width': 'mapData(degree, 1, 10, 30, 70)',
                    'height': 'mapData(degree, 1, 10, 30, 70)',
                    'label': 'data(label)',
                    'text-valign': 'bottom',
                    'text-halign': 'center',
                    'text-margin-y': 5,
                    'font-size': '11px',
                    'color': '#e8e8e8',
                    'text-outline-width': 2,
                    'text-outline-color': '#1a1a2e',
                    'border-width': 2,
                    'border-color': '#ffffff',
                    'border-opacity': 0.3,
                    'transition-property': 'background-color, border-color, width, height',
                    'transition-duration': '0.2s'
                }
            },
            // VMH metabolites
            {
                selector: 'node[type="metabolite"][subtype="vmh"]',
                style: {
                    'background-color': CONFIG.colors.vmhNode,
                    'background-opacity': 0.9
                }
            },
            // Saved metabolites
            {
                selector: 'node[type="metabolite"][subtype="saved"]',
                style: {
                    'background-color': CONFIG.colors.savedNode,
                    'background-opacity': 0.9
                }
            },
            // Reaction nodes
            {
                selector: 'node[type="reaction"]',
                style: {
                    'shape': 'diamond',
                    'width': 24,
                    'height': 24,
                    'label': 'data(label)',
                    'text-valign': 'bottom',
                    'text-halign': 'center',
                    'text-margin-y': 5,
                    'font-size': '10px',
                    'color': '#a0a0a0',
                    'text-outline-width': 2,
                    'text-outline-color': '#1a1a2e',
                    'border-width': 2,
                    'border-color': '#ffffff',
                    'border-opacity': 0.4,
                    'transition-property': 'background-color, border-color, width, height',
                    'transition-duration': '0.2s'
                }
            },
            // Balanced reactions
            {
                selector: 'node[type="reaction"][?balanced]',
                style: {
                    'background-color': CONFIG.colors.reactionBalanced
                }
            },
            // Unbalanced reactions
            {
                selector: 'node[type="reaction"][!balanced]',
                style: {
                    'background-color': CONFIG.colors.reactionUnbalanced
                }
            },
            // Edges - base style
            {
                selector: 'edge',
                style: {
                    'width': 2,
                    'line-color': CONFIG.colors.edge,
                    'curve-style': 'bezier',
                    'opacity': 0.6,
                    'transition-property': 'line-color, opacity, width',
                    'transition-duration': '0.2s'
                }
            },
            // Substrate edges (arrows pointing to reaction)
            {
                selector: 'edge[type="substrate"]',
                style: {
                    'target-arrow-shape': 'triangle',
                    'target-arrow-color': CONFIG.colors.edge,
                    'arrow-scale': 0.8
                }
            },
            // Product edges (arrows pointing to product)
            {
                selector: 'edge[type="product"]',
                style: {
                    'target-arrow-shape': 'triangle',
                    'target-arrow-color': CONFIG.colors.edge,
                    'arrow-scale': 0.8
                }
            },
            // Highlighted state
            {
                selector: 'node.highlighted',
                style: {
                    'border-color': CONFIG.colors.edgeHighlight,
                    'border-width': 4,
                    'border-opacity': 1
                }
            },
            {
                selector: 'edge.highlighted',
                style: {
                    'line-color': CONFIG.colors.edgeHighlight,
                    'target-arrow-color': CONFIG.colors.edgeHighlight,
                    'opacity': 1,
                    'width': 3
                }
            },
            // Faded state (when something else is highlighted)
            {
                selector: '.faded',
                style: {
                    'opacity': 0.15
                }
            },
            // Selected state
            {
                selector: 'node:selected',
                style: {
                    'border-color': '#ffd700',
                    'border-width': 4,
                    'border-opacity': 1
                }
            },
            // Search match
            {
                selector: '.search-match',
                style: {
                    'border-color': '#ffd700',
                    'border-width': 4,
                    'border-opacity': 1,
                    'z-index': 999
                }
            },
            // Labels hidden state
            {
                selector: '.labels-hidden',
                style: {
                    'label': ''
                }
            }
        ];
    }

    /**
     * Bind Cytoscape event handlers
     */
    function bindCytoscapeEvents() {
        // Hover effects for metabolite nodes
        cy.on('mouseover', 'node[type="metabolite"]', (e) => {
            const node = e.target;
            highlightConnected(node);
            showMetaboliteTooltip(node, e);
        });

        cy.on('mouseout', 'node[type="metabolite"]', () => {
            clearHighlight();
            hideTooltip();
        });

        // Hover effects for reaction nodes
        cy.on('mouseover', 'node[type="reaction"]', (e) => {
            const node = e.target;
            highlightConnected(node);
            showReactionTooltip(node, e);
        });

        cy.on('mouseout', 'node[type="reaction"]', () => {
            clearHighlight();
            hideTooltip();
        });

        // Click on metabolite - show info panel
        cy.on('tap', 'node[type="metabolite"]', (e) => {
            const node = e.target;
            showMetaboliteInfoPanel(node);
        });

        // Click on reaction - show info panel
        cy.on('tap', 'node[type="reaction"]', (e) => {
            const node = e.target;
            showReactionInfoPanel(node);
        });
        
        // Double-click on metabolite - open VMH if available
        cy.on('dbltap', 'node[type="metabolite"]', (e) => {
            const node = e.target;
            const data = node.data();
            if (data.subtype === 'vmh' && data.inVmh) {
                const abbr = data.label.replace(/\[.*\]$/, '');
                window.open(CONFIG.urls.vmhMetabolite(abbr), '_blank');
            }
        });
        
        // Double-click on reaction - open edit page
        cy.on('dbltap', 'node[type="reaction"]', (e) => {
            const node = e.target;
            const data = node.data();
            window.open(CONFIG.urls.editReaction(data.reactionDbId), '_blank');
        });

        // Right-click context menu
        cy.on('cxttap', 'node[type="metabolite"]', (e) => {
            showMetaboliteContextMenu(e);
        });

        cy.on('cxttap', 'node[type="reaction"]', (e) => {
            showReactionContextMenu(e);
        });

        // Click on background to deselect
        cy.on('tap', (e) => {
            if (e.target === cy) {
                hideInfoPanel();
            }
        });

        // Track mouse position for tooltip positioning
        cy.on('mousemove', (e) => {
            if (tooltip.classList.contains('visible')) {
                positionTooltip(e.originalEvent);
            }
        });
    }

    /**
     * Highlight a node and its connected elements
     */
    function highlightConnected(node) {
        // Get neighborhood (connected nodes and edges)
        const neighborhood = node.closedNeighborhood();
        
        // Fade all elements
        cy.elements().addClass('faded');
        
        // Highlight neighborhood
        neighborhood.removeClass('faded').addClass('highlighted');
    }

    /**
     * Clear all highlighting
     */
    function clearHighlight() {
        cy.elements().removeClass('faded highlighted');
    }

    /**
     * Show tooltip for metabolite node
     */
    function showMetaboliteTooltip(node, e) {
        const data = node.data();
        const isVmh = data.subtype === 'vmh';
        
        let html = `
            <div class="tooltip-header">
                <div class="tooltip-icon ${data.subtype}"></div>
                <span class="tooltip-title">${escapeHtml(data.label)}</span>
                <span class="tooltip-subtitle">${isVmh ? 'VMH' : 'New'}</span>
            </div>
        `;
        
        if (data.name && data.name !== data.label) {
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-label">Name</span>
                    <span class="tooltip-value">${escapeHtml(truncate(data.name, 50))}</span>
                </div>
            `;
        }

        if (data.formula) {
            html += `<div class="tooltip-formula">${formatFormula(data.formula)}</div>`;
        }
        
        if (data.compartment) {
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-label">Compartment</span>
                    <span class="tooltip-value">${escapeHtml(data.compartment)}</span>
                </div>
            `;
        }

        if (!isVmh && data.source) {
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-label">Source</span>
                    <span class="tooltip-value">${escapeHtml(data.source)}</span>
                </div>
            `;
        }

        const connections = node.connectedEdges().length;
        html += `
            <div class="tooltip-row">
                <span class="tooltip-label">Connections</span>
                <span class="tooltip-value">${connections} reaction${connections !== 1 ? 's' : ''}</span>
            </div>
        `;

        if (isVmh && data.inVmh) {
            html += `
                <div class="tooltip-action">
                    <i class="fa fa-external-link-alt"></i>
                    Click to view on VMH
                </div>
            `;
        }

        tooltip.innerHTML = html;
        tooltip.classList.add('visible');
        positionTooltip(e.originalEvent);
    }

    /**
     * Show tooltip for reaction node
     */
    function showReactionTooltip(node, e) {
        const data = node.data();
        
        // Determine balanced status text
        let balancedText = 'Unknown';
        let balancedClass = '';
        if (data.balancedRaw === true) {
            balancedText = 'Balanced';
            balancedClass = 'balanced';
        } else if (data.balancedRaw === false) {
            balancedText = 'Unbalanced';
            balancedClass = 'unbalanced';
        }
        
        let html = `
            <div class="tooltip-header">
                <div class="tooltip-icon reaction" style="background: ${data.balanced ? CONFIG.colors.reactionBalanced : CONFIG.colors.reactionUnbalanced}"></div>
                <span class="tooltip-title">${escapeHtml(data.abbreviation || data.label)}</span>
                <span class="tooltip-subtitle">${balancedText}</span>
            </div>
        `;

        if (data.description) {
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-value" style="text-align: left; max-width: 100%;">${escapeHtml(truncate(data.description, 100))}</span>
                </div>
            `;
        }

        html += `
            <div class="tooltip-row">
                <span class="tooltip-label">Substrates</span>
                <span class="tooltip-value">${data.substrateCount}</span>
            </div>
            <div class="tooltip-row">
                <span class="tooltip-label">Products</span>
                <span class="tooltip-value">${data.productCount}</span>
            </div>
        `;
        
        if (data.subsystem) {
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-label">Subsystem</span>
                    <span class="tooltip-value">${escapeHtml(truncate(data.subsystem, 30))}</span>
                </div>
            `;
        }

        if (data.confidenceScore !== undefined && data.confidenceScore !== null) {
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-label">Confidence</span>
                    <span class="tooltip-value">${Math.round(data.confidenceScore * 100)}%</span>
                </div>
            `;
        }
        
        if (data.flags && data.flags.length > 0) {
            const flagNames = data.flags.map(f => f.name).join(', ');
            html += `
                <div class="tooltip-row">
                    <span class="tooltip-label">Flags</span>
                    <span class="tooltip-value">${escapeHtml(truncate(flagNames, 40))}</span>
                </div>
            `;
        }

        html += `
            <div class="tooltip-action">
                <i class="fa fa-edit"></i>
                Click to view details • Right-click for options
            </div>
        `;

        tooltip.innerHTML = html;
        tooltip.classList.add('visible');
        positionTooltip(e.originalEvent);
    }

    /**
     * Position tooltip near cursor
     */
    function positionTooltip(e) {
        const rect = document.getElementById('graphModalOverlay').getBoundingClientRect();
        const tooltipRect = tooltip.getBoundingClientRect();
        
        let x = e.clientX - rect.left + 15;
        let y = e.clientY - rect.top + 15;
        
        // Keep tooltip within bounds
        if (x + tooltipRect.width > rect.width - 20) {
            x = e.clientX - rect.left - tooltipRect.width - 15;
        }
        if (y + tooltipRect.height > rect.height - 20) {
            y = e.clientY - rect.top - tooltipRect.height - 15;
        }
        
        tooltip.style.left = x + 'px';
        tooltip.style.top = y + 'px';
    }

    /**
     * Hide tooltip
     */
    function hideTooltip() {
        tooltip.classList.remove('visible');
    }

    /**
     * Show context menu for metabolite
     */
    function showMetaboliteContextMenu(e) {
        const node = e.target;
        const data = node.data();
        const isVmh = data.subtype === 'vmh';
        
        let html = '';
        
        if (isVmh && data.inVmh) {
            // Extract abbr without compartment for VMH URL
            const abbr = data.label.replace(/\[.*\]$/, '');
            html += `
                <div class="context-menu-item" data-action="vmh" data-id="${escapeHtml(abbr)}">
                    <i class="fa fa-external-link-alt"></i>
                    View on VMH
                </div>
            `;
        }
        
        html += `
            <div class="context-menu-item" data-action="focus" data-id="${data.id}">
                <i class="fa fa-crosshairs"></i>
                Focus on this metabolite
            </div>
            <div class="context-menu-item" data-action="neighbors" data-id="${data.id}">
                <i class="fa fa-project-diagram"></i>
                Show connected reactions
            </div>
        `;

        showContextMenuAt(e, html, node);
    }

    /**
     * Show context menu for reaction
     */
    function showReactionContextMenu(e) {
        const node = e.target;
        const data = node.data();
        
        const html = `
            <div class="context-menu-item" data-action="edit" data-id="${data.reactionDbId}">
                <i class="fa fa-edit"></i>
                Edit reaction
            </div>
            <div class="context-menu-item" data-action="focus" data-id="${data.id}">
                <i class="fa fa-crosshairs"></i>
                Focus on this reaction
            </div>
            <div class="context-menu-divider"></div>
            <div class="context-menu-item" data-action="neighbors" data-id="${data.id}">
                <i class="fa fa-project-diagram"></i>
                Show connected metabolites
            </div>
        `;

        showContextMenuAt(e, html, node);
    }

    /**
     * Show context menu at position
     */
    function showContextMenuAt(e, html, node) {
        selectedElement = node;
        contextMenu.innerHTML = html;
        contextMenu.classList.add('visible');
        
        const rect = document.getElementById('graphModalOverlay').getBoundingClientRect();
        const pos = e.position || e.cyPosition;
        const rendered = e.renderedPosition;
        
        let x = rendered.x + rect.left;
        let y = rendered.y + rect.top;
        
        // Adjust for header
        y += 60; // Approximate header height
        
        // Keep within bounds
        contextMenu.style.left = x + 'px';
        contextMenu.style.top = y + 'px';
        
        // Bind click handlers
        contextMenu.querySelectorAll('.context-menu-item').forEach(item => {
            item.addEventListener('click', handleContextMenuClick);
        });
    }

    /**
     * Handle context menu item click
     */
    function handleContextMenuClick(e) {
        const action = e.currentTarget.dataset.action;
        const id = e.currentTarget.dataset.id;
        
        hideContextMenu();
        
        switch (action) {
            case 'edit':
                window.open(CONFIG.urls.editReaction(id), '_blank');
                break;
            case 'vmh':
                window.open(CONFIG.urls.vmhMetabolite(id), '_blank');
                break;
            case 'focus':
                focusOnNode(id);
                break;
            case 'neighbors':
                showNeighbors(id);
                break;
        }
    }

    /**
     * Hide context menu
     */
    function hideContextMenu() {
        contextMenu.classList.remove('visible');
        selectedElement = null;
    }

    /**
     * Focus on a specific node
     */
    function focusOnNode(nodeId) {
        const node = cy.getElementById(nodeId);
        if (node.length) {
            cy.animate({
                center: { eles: node },
                zoom: 2
            }, {
                duration: 500
            });
        }
    }

    /**
     * Show only neighbors of a node
     */
    function showNeighbors(nodeId) {
        const node = cy.getElementById(nodeId);
        if (node.length) {
            const neighborhood = node.closedNeighborhood();
            cy.elements().addClass('faded');
            neighborhood.removeClass('faded');
            
            cy.animate({
                fit: { eles: neighborhood, padding: 50 }
            }, {
                duration: 500
            });
            
            // Add reset button if not exists
            setTimeout(() => {
                clearHighlight();
            }, 3000);
        }
    }

    /**
     * Show metabolite info in side panel
     */
    function showMetaboliteInfoPanel(node) {
        const data = node.data();
        const isVmh = data.subtype === 'vmh';
        
        const panelTitle = document.getElementById('infoPanelTitle');
        const panelContent = document.getElementById('infoPanelContent');
        
        panelTitle.textContent = data.label;
        
        let html = `
            <div class="info-section">
                <div class="info-row">
                    <span class="info-label">Type</span>
                    <span class="info-badge ${data.subtype}">${isVmh ? 'VMH' : 'New'}</span>
                </div>
        `;
        
        if (data.name && data.name !== data.label) {
            html += `
                <div class="info-row">
                    <span class="info-label">Full Name</span>
                    <span class="info-value">${escapeHtml(data.name)}</span>
                </div>
            `;
        }
        
        if (data.formula) {
            html += `
                <div class="info-row">
                    <span class="info-label">Formula</span>
                    <span class="info-value">${formatFormula(data.formula)}</span>
                </div>
            `;
        }
        
        if (data.compartment) {
            html += `
                <div class="info-row">
                    <span class="info-label">Compartment</span>
                    <span class="info-value">${escapeHtml(data.compartment)}</span>
                </div>
            `;
        }
        
        if (!isVmh && data.source) {
            html += `
                <div class="info-row">
                    <span class="info-label">Source</span>
                    <span class="info-value">${escapeHtml(data.source)}</span>
                </div>
            `;
        }
        
        if (data.inchiKey) {
            html += `
                <div class="info-row">
                    <span class="info-label">InChI Key</span>
                    <span class="info-value" style="font-size: 10px;">${escapeHtml(data.inchiKey)}</span>
                </div>
            `;
        }
        
        html += `</div>`;
        
        // Connected reactions section
        const connectedReactions = node.connectedEdges()
            .connectedNodes()
            .filter('[type="reaction"]');
        
        if (connectedReactions.length > 0) {
            html += `
                <div class="info-section">
                    <div class="info-section-title">Connected Reactions (${connectedReactions.length})</div>
            `;
            
            connectedReactions.forEach(rxn => {
                const rxnData = rxn.data();
                const balanceIcon = rxnData.balancedRaw === true ? '✓' : (rxnData.balancedRaw === false ? '!' : '?');
                const balanceClass = rxnData.balancedRaw === true ? 'balanced' : (rxnData.balancedRaw === false ? 'unbalanced' : '');
                html += `
                    <div class="info-row">
                        <a class="info-link" href="${CONFIG.urls.editReaction(rxnData.reactionDbId)}" target="_blank">
                            ${escapeHtml(rxnData.abbreviation || rxnData.label)}
                        </a>
                        <span class="info-badge ${balanceClass}">
                            ${balanceIcon}
                        </span>
                    </div>
                `;
            });
            
            html += `</div>`;
        }
        
        // VMH link
        if (isVmh && data.inVmh) {
            const abbr = data.label.replace(/\[.*\]$/, ''); // Remove compartment suffix if present
            html += `
                <div class="info-section">
                    <a class="info-link" href="${CONFIG.urls.vmhMetabolite(abbr)}" target="_blank">
                        <i class="fa fa-external-link-alt"></i> View on VMH
                    </a>
                </div>
            `;
        }
        
        panelContent.innerHTML = html;
        showInfoPanel();
    }

    /**
     * Show reaction info in side panel
     */
    function showReactionInfoPanel(node) {
        const data = node.data();
        
        const panelTitle = document.getElementById('infoPanelTitle');
        const panelContent = document.getElementById('infoPanelContent');
        
        panelTitle.textContent = data.abbreviation || data.label;
        
        // Determine balanced status
        let balancedText = 'Unknown';
        let balancedClass = '';
        if (data.balancedRaw === true) {
            balancedText = 'Balanced';
            balancedClass = 'balanced';
        } else if (data.balancedRaw === false) {
            balancedText = 'Unbalanced';
            balancedClass = 'unbalanced';
        }
        
        let html = `
            <div class="info-section">
                <div class="info-row">
                    <span class="info-label">Status</span>
                    <span class="info-badge ${balancedClass}">
                        ${balancedText}
                    </span>
                </div>
        `;
        
        if (data.confidenceScore !== undefined && data.confidenceScore !== null) {
            html += `
                <div class="info-row">
                    <span class="info-label">Confidence</span>
                    <span class="info-value">${Math.round(data.confidenceScore * 100)}%</span>
                </div>
            `;
        }
        
        html += `
                <div class="info-row">
                    <span class="info-label">Substrates</span>
                    <span class="info-value">${data.substrateCount}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Products</span>
                    <span class="info-value">${data.productCount}</span>
                </div>
        `;
        
        if (data.direction) {
            html += `
                <div class="info-row">
                    <span class="info-label">Direction</span>
                    <span class="info-value">${escapeHtml(data.direction)}${data.reversible ? ' (reversible)' : ''}</span>
                </div>
            `;
        }
        
        if (data.subsystem) {
            html += `
                <div class="info-row">
                    <span class="info-label">Subsystem</span>
                    <span class="info-value">${escapeHtml(data.subsystem)}</span>
                </div>
            `;
        }
        
        html += `</div>`;
        
        if (data.description) {
            html += `
                <div class="info-section">
                    <div class="info-section-title">Description</div>
                    <p style="color: #e8e8e8; font-size: 13px; margin: 0; line-height: 1.5;">
                        ${escapeHtml(data.description)}
                    </p>
                </div>
            `;
        }
        
        // Flags
        if (data.flags && data.flags.length > 0) {
            html += `
                <div class="info-section">
                    <div class="info-section-title">Flags</div>
            `;
            data.flags.forEach(flag => {
                html += `
                    <div class="info-row">
                        <span style="display: flex; align-items: center; gap: 6px;">
                            <span style="width: 10px; height: 10px; border-radius: 50%; background: ${flag.color || '#888'};"></span>
                            <span class="info-value">${escapeHtml(flag.name)}</span>
                        </span>
                    </div>
                `;
            });
            html += `</div>`;
        }
        
        // Substrates
        const substrates = node.incomers('node[type="metabolite"]');
        if (substrates.length > 0) {
            html += `
                <div class="info-section">
                    <div class="info-section-title">Substrates</div>
            `;
            substrates.forEach(met => {
                const metData = met.data();
                html += `
                    <div class="info-row">
                        <span class="info-value">${escapeHtml(metData.label)}</span>
                        <span class="info-badge ${metData.subtype}">${metData.subtype === 'vmh' ? 'VMH' : 'New'}</span>
                    </div>
                `;
            });
            html += `</div>`;
        }
        
        // Products
        const products = node.outgoers('node[type="metabolite"]');
        if (products.length > 0) {
            html += `
                <div class="info-section">
                    <div class="info-section-title">Products</div>
            `;
            products.forEach(met => {
                const metData = met.data();
                html += `
                    <div class="info-row">
                        <span class="info-value">${escapeHtml(metData.label)}</span>
                        <span class="info-badge ${metData.subtype}">${metData.subtype === 'vmh' ? 'VMH' : 'New'}</span>
                    </div>
                `;
            });
            html += `</div>`;
        }
        
        // Edit link
        html += `
            <div class="info-section">
                <a class="info-link" href="${CONFIG.urls.editReaction(data.reactionDbId)}" target="_blank">
                    <i class="fa fa-edit"></i> Edit this reaction
                </a>
            </div>
        `;
        
        panelContent.innerHTML = html;
        showInfoPanel();
    }

    /**
     * Show info panel
     */
    function showInfoPanel() {
        document.getElementById('infoPanel').classList.add('visible');
        isInfoPanelVisible = true;
    }

    /**
     * Hide info panel
     */
    function hideInfoPanel() {
        document.getElementById('infoPanel').classList.remove('visible');
        isInfoPanelVisible = false;
    }

    /**
     * Run graph layout - uses the default layout or current layout panel settings
     */
    function runLayout() {
        if (!cy) return;
        
        // If layout panel has been used, apply those settings
        // Otherwise use the default optimized layout
        const layoutOptions = currentLayoutConfig.algorithm !== 'cose-bilkent' 
            ? buildLayoutOptions() 
            : CONFIG.defaultLayout;
            
        const layout = cy.layout(layoutOptions);
        layout.run();
    }

    /**
     * Re-layout the graph
     */
    function relayout() {
        const btn = document.getElementById('btnRelayout');
        btn.classList.add('active');
        
        // Use current layout settings from panel
        applyCurrentLayout();
        
        setTimeout(() => {
            btn.classList.remove('active');
        }, 1200);
    }

    /**
     * Toggle node labels
     */
    let labelsVisible = true;
    function toggleLabels() {
        const btn = document.getElementById('btnToggleLabels');
        labelsVisible = !labelsVisible;
        
        if (labelsVisible) {
            cy.elements().removeClass('labels-hidden');
            btn.classList.remove('active');
        } else {
            cy.elements().addClass('labels-hidden');
            btn.classList.add('active');
        }
    }

    /**
     * Handle search input
     */
    function handleSearch(e) {
        const query = e.target.value.trim().toLowerCase();
        
        // Clear previous search highlights
        cy.elements().removeClass('search-match faded');
        
        if (!query) return;
        
        // Find matching nodes (search in label, name, formula, description for reactions)
        const matches = cy.nodes().filter(node => {
            const label = (node.data('label') || '').toLowerCase();
            const name = (node.data('name') || '').toLowerCase();
            const formula = (node.data('formula') || '').toLowerCase();
            const description = (node.data('description') || '').toLowerCase();
            const abbreviation = (node.data('abbreviation') || '').toLowerCase();
            
            return label.includes(query) || 
                   name.includes(query) || 
                   formula.includes(query) ||
                   description.includes(query) ||
                   abbreviation.includes(query);
        });
        
        if (matches.length > 0) {
            cy.elements().addClass('faded');
            matches.removeClass('faded').addClass('search-match');
            
            // Also show connected elements
            matches.closedNeighborhood().removeClass('faded');
            
            // Fit to matches if there are few results
            if (matches.length <= 10) {
                cy.animate({
                    fit: { eles: matches.closedNeighborhood(), padding: 50 }
                }, {
                    duration: 500
                });
            }
        }
    }

    /**
     * Download graph as PNG
     */
    function downloadScreenshot() {
        if (!cy) return;
        
        const png = cy.png({
            output: 'blob',
            bg: CONFIG.colors.background,
            scale: 2,
            full: true
        });
        
        const link = document.createElement('a');
        link.href = URL.createObjectURL(png);
        link.download = 'reaction-network.png';
        link.click();
        URL.revokeObjectURL(link.href);
    }

    /**
     * Show/hide loading overlay
     */
    function showLoading(show) {
        const loading = document.getElementById('graphLoading');
        loading.style.display = show ? 'flex' : 'none';
    }

    /**
     * Show error message
     */
    function showError(message) {
        const loading = document.getElementById('graphLoading');
        loading.innerHTML = `
            <i class="fa fa-exclamation-triangle" style="font-size: 48px; color: #ff6b6b;"></i>
            <div class="graph-loading-text">${escapeHtml(message)}</div>
            <button class="graph-btn" onclick="window.GraphVisualization.closeGraph()" style="margin-top: 16px;">
                Close
            </button>
        `;
    }
    
    /**
     * Show empty state when no reactions exist
     */
    function showEmptyState() {
        const loading = document.getElementById('graphLoading');
        loading.innerHTML = `
            <i class="fa fa-flask" style="font-size: 64px; color: #667eea; opacity: 0.5;"></i>
            <div class="graph-loading-text" style="margin-top: 20px; font-size: 18px; color: #e8e8e8;">
                No reactions to display
            </div>
            <div class="graph-loading-text" style="margin-top: 8px; color: #8a8aaa; font-size: 14px;">
                Save some reactions first to see them visualized here.
            </div>
            <button class="graph-btn" onclick="window.GraphVisualization.closeGraph()" style="margin-top: 24px;">
                <i class="fa fa-arrow-left"></i> Go Back
            </button>
        `;
        loading.style.display = 'flex';
    }

    /**
     * Update stats display
     */
    function updateStats(stats) {
        document.getElementById('statVmh').textContent = stats.vmh_metabolites || 0;
        document.getElementById('statSaved').textContent = stats.saved_metabolites || 0;
        document.getElementById('statReactions').textContent = stats.edge_count || 0;
    }

    // Utility functions
    function escapeHtml(text) {
        if (!text) return '';
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }

    function formatFormula(formula) {
        if (!formula) return '';
        // Convert numbers to subscript
        return formula.replace(/(\d+)/g, '<sub>$1</sub>');
    }

    function truncate(str, maxLength) {
        if (!str || str.length <= maxLength) return str;
        return str.substring(0, maxLength) + '...';
    }

    function debounce(func, wait) {
        let timeout;
        return function(...args) {
            clearTimeout(timeout);
            timeout = setTimeout(() => func.apply(this, args), wait);
        };
    }

    // Initialize on DOM ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }

    // Expose public API
    window.GraphVisualization = {
        open: openGraph,
        close: closeGraph,
        closeGraph: closeGraph
    };

})();
