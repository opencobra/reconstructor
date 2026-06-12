/**
 * tabsController.js — header tab strip where each tab is a ReactionGroup
 * (001-multi-reaction-tabbed, US2; FR-005..FR-009).
 *
 * Renders the user's groups as tabs, marks the active one, and switches between
 * them by fetching the target group's contents and re-rendering the group view
 * via GroupView — no full page reload (FR-006 / SC-002). Create / rename / delete
 * call the group endpoints; active group + membership persist per user
 * server-side, so reloading restores the last state (FR-007). At least one
 * default tab is always present (FR-008). Anonymous users get one ephemeral tab
 * and the existing single-reaction behaviour (FR-018).
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

    var TabsController = {
        groups: [],
        activeGroupId: null,

        strip: function () { return document.getElementById('groupTabs'); },

        init: function () {
            if (!this.strip()) return;
            if (!userId()) {
                this.renderEphemeral();
                if (global.GroupView) GroupView.renderEmpty(true);
                return;
            }
            this.reload(true);
        },

        /** Fetch groups and render tabs. If `loadActive`, also load active group. */
        reload: function (loadActive) {
            var self = this;
            var url = global.groupsListUrl + '?userID=' + encodeURIComponent(userId());
            return fetch(url)
                .then(function (r) { return r.json(); })
                .then(function (data) {
                    if (data.status === 'anonymous') {
                        self.renderEphemeral();
                        return;
                    }
                    self.groups = data.groups || [];
                    self.activeGroupId = data.active_group_id;
                    self.render();
                    if (loadActive && global.GroupView) {
                        // If the URL already names a reaction, the page load path
                        // opens it; we just render the group rail around it.
                        GroupView.loadGroup(self.activeGroupId);
                    }
                });
        },

        render: function () {
            var strip = this.strip();
            if (!strip) return;
            strip.innerHTML = '';
            var self = this;
            this.groups.forEach(function (g) {
                strip.appendChild(self.buildTab(g));
            });
            strip.appendChild(this.buildAddButton());
        },

        buildTab: function (group) {
            var isActive = String(group.id) === String(this.activeGroupId);
            var tab = document.createElement('div');
            tab.className = 'group-tab' + (isActive ? ' active' : '');
            tab.setAttribute('data-group-id', group.id);
            tab.setAttribute('role', 'tab');
            tab.setAttribute('tabindex', '0');
            tab.setAttribute('aria-selected', isActive ? 'true' : 'false');

            tab.innerHTML =
                '<span class="group-tab-name"></span>' +
                '<span class="group-tab-count">' + (group.count || 0) + '</span>' +
                (isActive
                    ? '<button type="button" class="group-tab-action group-tab-rename" title="Rename group" aria-label="Rename group"><i class="fas fa-pen" aria-hidden="true"></i></button>' +
                      '<button type="button" class="group-tab-action group-tab-delete" title="Delete group" aria-label="Delete group"><i class="fas fa-trash" aria-hidden="true"></i></button>'
                    : '');
            tab.querySelector('.group-tab-name').textContent = group.name;

            var self = this;
            tab.addEventListener('click', function (e) {
                if (e.target.closest('.group-tab-action')) return;
                self.switchTo(group.id);
            });
            tab.addEventListener('keydown', function (e) {
                if (e.key === 'Enter' || e.key === ' ') {
                    e.preventDefault();
                    self.switchTo(group.id);
                }
            });
            var renameBtn = tab.querySelector('.group-tab-rename');
            if (renameBtn) {
                renameBtn.addEventListener('click', function (e) {
                    e.stopPropagation();
                    self.renameGroup(group);
                });
            }
            var deleteBtn = tab.querySelector('.group-tab-delete');
            if (deleteBtn) {
                deleteBtn.addEventListener('click', function (e) {
                    e.stopPropagation();
                    self.deleteGroup(group);
                });
            }
            return tab;
        },

        buildAddButton: function () {
            var btn = document.createElement('button');
            btn.type = 'button';
            btn.className = 'group-tab-add';
            btn.title = 'New group';
            btn.setAttribute('aria-label', 'New group');
            btn.innerHTML = '<i class="fas fa-plus" aria-hidden="true"></i>';
            var self = this;
            btn.addEventListener('click', function () { self.createGroup(); });
            return btn;
        },

        /** Switch active group without a page reload (FR-006). */
        switchTo: function (groupId) {
            if (String(groupId) === String(this.activeGroupId)) return;
            var self = this;
            var proceed = global.GroupView ? GroupView.guardUnsaved() : Promise.resolve(true);
            proceed.then(function (ok) {
                if (!ok) return;
                self.activeGroupId = groupId;
                self.render();
                if (global.GroupView) GroupView.loadGroup(groupId);
            });
        },

        createGroup: function () {
            if (!userId()) return;
            var self = this;
            Notify.prompt({
                title: 'New group',
                message: 'Name for the new group:',
                defaultValue: 'New group',
                placeholder: 'Group name',
                confirmText: 'Create',
            }).then(function (name) {
                if (name === null || !name.trim()) return;
                postJSON(global.groupCreateUrl, { userID: userId(), name: name.trim() })
                    .then(function (res) {
                        if (res.status === 'success') {
                            self.activeGroupId = res.group.id;
                            self.reload(false).then(function () {
                                if (global.GroupView) GroupView.loadGroup(self.activeGroupId);
                            });
                        }
                    });
            });
        },

        renameGroup: function (group) {
            var self = this;
            Notify.prompt({
                title: 'Rename group',
                message: 'Enter a new name for this group:',
                defaultValue: group.name,
                placeholder: 'Group name',
                confirmText: 'Rename',
            }).then(function (name) {
                if (name === null || !name.trim()) return;
                postJSON(global.groupRenameUrl, {
                    userID: userId(), groupId: group.id, name: name.trim(),
                }).then(function (res) {
                    if (res.status === 'success') self.reload(false);
                });
            });
        },

        deleteGroup: function (group) {
            if (this.groups.length <= 1) {
                Notify.warning('You must keep at least one group.');  // FR-008
                return;
            }
            var self = this;
            Notify.confirm({
                title: 'Delete group',
                message: 'Delete group "' + group.name + '"?\n' +
                    'The reactions in it stay in your Saved Reactions.',  // FR-009
                confirmText: 'Delete group',
                cancelText: 'Cancel',
                danger: true,
            }).then(function (ok) {
                if (!ok) return;
                postJSON(global.groupDeleteUrl, { userID: userId(), groupId: group.id })
                    .then(function (res) {
                        if (res.status === 'success') {
                            self.activeGroupId = res.active_group_id;
                            self.reload(false).then(function () {
                                if (global.GroupView) GroupView.loadGroup(self.activeGroupId);
                            });
                        }
                    });
            });
        },

        /** Refresh tab count badges after add/remove (keeps tabs in step). */
        refreshCounts: function () {
            if (!userId()) return;
            this.reload(false);
        },

        renderEphemeral: function () {
            var strip = this.strip();
            if (!strip) return;
            strip.innerHTML =
                '<div class="group-tab active ephemeral" role="tab" aria-selected="true">' +
                '<span class="group-tab-name">Workspace</span></div>';
        },
    };

    global.TabsController = TabsController;

    document.addEventListener('DOMContentLoaded', function () {
        // Defer slightly so sessionStorage userID (set by reactantsdisplay.js) and
        // CSRF globals are in place.
        setTimeout(function () { TabsController.init(); }, 0);
    });
})(window);
