/**
 * notify.js — neat, professional toast notifications + confirm/prompt modals
 * that replace the browser-native alert()/confirm()/prompt() popups across the
 * workspace.
 *
 * API (all on window.Notify):
 *   Notify.success(message, title?)   → green toast
 *   Notify.error(message, title?)     → red toast
 *   Notify.warning(message, title?)   → amber toast
 *   Notify.info(message, title?)      → blue toast
 *   Notify.toast(message, opts)       → opts: { type, title, duration }
 *   Notify.confirm(opts) → Promise<boolean>
 *       opts: { title, message, confirmText, cancelText, danger }
 *   Notify.prompt(opts)  → Promise<string|null>  (null = cancelled)
 *       opts: { title, message, defaultValue, placeholder, confirmText, cancelText }
 *
 * The module is self-bootstrapping: it injects its own toast container and modal
 * markup on first use, so it works on any page that loads this script. Styling
 * lives in home_page.css (.notify-* / .confirm-modal-*).
 */
(function (global) {
    'use strict';

    var ICONS = {
        success: 'fas fa-check-circle',
        error: 'fas fa-exclamation-circle',
        warning: 'fas fa-exclamation-triangle',
        info: 'fas fa-info-circle',
    };

    function ensureToastContainer() {
        var container = document.getElementById('notifyContainer');
        if (!container) {
            container = document.createElement('div');
            container.id = 'notifyContainer';
            container.className = 'notify-container';
            container.setAttribute('role', 'status');
            container.setAttribute('aria-live', 'polite');
            document.body.appendChild(container);
        }
        return container;
    }

    function toast(message, opts) {
        opts = opts || {};
        var type = opts.type || 'info';
        var duration = opts.duration != null ? opts.duration
            : (type === 'error' ? 6000 : 4000);
        var container = ensureToastContainer();

        var el = document.createElement('div');
        el.className = 'notify-toast notify-' + type;
        el.setAttribute('role', type === 'error' ? 'alert' : 'status');

        var icon = document.createElement('i');
        icon.className = 'notify-icon ' + (ICONS[type] || ICONS.info);
        icon.setAttribute('aria-hidden', 'true');
        el.appendChild(icon);

        var body = document.createElement('div');
        body.className = 'notify-body';
        if (opts.title) {
            var t = document.createElement('div');
            t.className = 'notify-title';
            t.textContent = opts.title;
            body.appendChild(t);
        }
        var m = document.createElement('div');
        m.className = 'notify-message';
        m.textContent = message == null ? '' : String(message);
        body.appendChild(m);
        el.appendChild(body);

        var close = document.createElement('button');
        close.type = 'button';
        close.className = 'notify-close';
        close.setAttribute('aria-label', 'Dismiss');
        close.innerHTML = '<i class="fas fa-times" aria-hidden="true"></i>';
        el.appendChild(close);

        var timer = null;
        var dismiss = function () {
            if (timer) clearTimeout(timer);
            el.classList.add('notify-leaving');
            el.addEventListener('animationend', function () {
                if (el.parentNode) el.parentNode.removeChild(el);
            });
            // Fallback removal in case animationend doesn't fire.
            setTimeout(function () { if (el.parentNode) el.parentNode.removeChild(el); }, 400);
        };
        close.addEventListener('click', dismiss);

        container.appendChild(el);
        // Trigger enter animation on next frame.
        requestAnimationFrame(function () { el.classList.add('notify-in'); });

        if (duration > 0) {
            timer = setTimeout(dismiss, duration);
            el.addEventListener('mouseenter', function () { if (timer) clearTimeout(timer); });
            el.addEventListener('mouseleave', function () { timer = setTimeout(dismiss, 1500); });
        }
        return el;
    }

    // ---- Confirm / prompt modal ------------------------------------------
    function ensureModal() {
        var overlay = document.getElementById('notifyModalOverlay');
        if (overlay) return overlay;
        overlay = document.createElement('div');
        overlay.id = 'notifyModalOverlay';
        overlay.className = 'confirm-modal-overlay';
        overlay.innerHTML =
            '<div class="confirm-modal" role="dialog" aria-modal="true" aria-labelledby="confirmModalTitle">' +
            '<div class="confirm-modal-header">' +
            '<i class="confirm-modal-icon fas fa-question-circle" aria-hidden="true"></i>' +
            '<h3 class="confirm-modal-title" id="confirmModalTitle"></h3>' +
            '</div>' +
            '<div class="confirm-modal-message"></div>' +
            '<div class="confirm-modal-field" style="display:none">' +
            '<input type="text" class="confirm-modal-input" />' +
            '</div>' +
            '<div class="confirm-modal-actions">' +
            '<button type="button" class="confirm-modal-cancel ui button">Cancel</button>' +
            '<button type="button" class="confirm-modal-confirm ui button">Confirm</button>' +
            '</div>' +
            '</div>';
        document.body.appendChild(overlay);
        return overlay;
    }

    function openModal(opts, isPrompt) {
        opts = opts || {};
        return new Promise(function (resolve) {
            var overlay = ensureModal();
            var modal = overlay.querySelector('.confirm-modal');
            var titleEl = overlay.querySelector('.confirm-modal-title');
            var msgEl = overlay.querySelector('.confirm-modal-message');
            var iconEl = overlay.querySelector('.confirm-modal-icon');
            var fieldEl = overlay.querySelector('.confirm-modal-field');
            var inputEl = overlay.querySelector('.confirm-modal-input');
            var cancelBtn = overlay.querySelector('.confirm-modal-cancel');
            var confirmBtn = overlay.querySelector('.confirm-modal-confirm');

            titleEl.textContent = opts.title || (isPrompt ? 'Enter a value' : 'Please confirm');
            msgEl.textContent = opts.message || '';
            msgEl.style.display = opts.message ? '' : 'none';
            cancelBtn.textContent = opts.cancelText || 'Cancel';
            confirmBtn.textContent = opts.confirmText || (isPrompt ? 'OK' : 'Confirm');

            var danger = !!opts.danger;
            confirmBtn.className = 'confirm-modal-confirm ui button' + (danger ? ' danger' : '');
            iconEl.className = 'confirm-modal-icon ' + (danger
                ? 'fas fa-exclamation-triangle' : (isPrompt ? 'fas fa-pen' : 'fas fa-question-circle'));

            if (isPrompt) {
                fieldEl.style.display = '';
                inputEl.value = opts.defaultValue != null ? opts.defaultValue : '';
                inputEl.setAttribute('placeholder', opts.placeholder || '');
            } else {
                fieldEl.style.display = 'none';
            }

            var closed = false;
            function cleanup() {
                overlay.classList.remove('open');
                document.removeEventListener('keydown', onKey);
            }
            function done(value) {
                if (closed) return;
                closed = true;
                cleanup();
                resolve(value);
            }
            function onKey(e) {
                if (e.key === 'Escape') done(isPrompt ? null : false);
                else if (e.key === 'Enter' && (isPrompt || document.activeElement !== cancelBtn)) {
                    e.preventDefault();
                    done(isPrompt ? inputEl.value : true);
                }
            }

            cancelBtn.onclick = function () { done(isPrompt ? null : false); };
            confirmBtn.onclick = function () { done(isPrompt ? inputEl.value : true); };
            overlay.onclick = function (e) { if (e.target === overlay) done(isPrompt ? null : false); };
            document.addEventListener('keydown', onKey);

            overlay.classList.add('open');
            requestAnimationFrame(function () {
                modal.classList.add('confirm-modal-in');
                if (isPrompt) { inputEl.focus(); inputEl.select(); }
                else confirmBtn.focus();
            });
        });
    }

    global.Notify = {
        toast: toast,
        success: function (msg, title) { return toast(msg, { type: 'success', title: title }); },
        error: function (msg, title) { return toast(msg, { type: 'error', title: title }); },
        warning: function (msg, title) { return toast(msg, { type: 'warning', title: title }); },
        info: function (msg, title) { return toast(msg, { type: 'info', title: title }); },
        confirm: function (opts) {
            if (typeof opts === 'string') opts = { message: opts };
            return openModal(opts, false);
        },
        prompt: function (opts) {
            if (typeof opts === 'string') opts = { message: opts };
            return openModal(opts, true);
        },
    };
})(window);
