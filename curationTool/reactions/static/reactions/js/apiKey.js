// API key management page (feature 002): generate/rotate the session user's key.
$(function () {
    $('#generateKeyBtn').on('click', async function () {
        const btn = $(this);
        btn.addClass('loading disabled');
        $('#keyError').hide();
        try {
            const response = await fetch('/api-key/generate/', {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': window.csrfToken,
                },
            });
            const data = await response.json();
            if (response.ok && data.status === 'success') {
                $('#newKeyValue').val(data.api_key);
                $('#newKeySegment').show();
                document.getElementById('newKeyValue').select();
            } else {
                $('#keyError').text(data.message || 'Failed to generate key.').show();
            }
        } catch (error) {
            $('#keyError').text('Failed to generate key: ' + error).show();
        } finally {
            btn.removeClass('loading disabled');
        }
    });

    let copyFeedbackTimer = null;

    async function copyNewKey() {
        const input = document.getElementById('newKeyValue');
        if (!input.value) return;
        input.select();
        input.setSelectionRange(0, input.value.length);
        let copied = false;
        if (navigator.clipboard && window.isSecureContext) {
            try {
                await navigator.clipboard.writeText(input.value);
                copied = true;
            } catch (e) { /* fall through to execCommand */ }
        }
        if (!copied) copied = document.execCommand('copy');
        const btn = $('#copyKeyBtn');
        if (copied) {
            btn.html('<i class="check icon"></i> Copied');
            clearTimeout(copyFeedbackTimer);
            copyFeedbackTimer = setTimeout(function () {
                btn.html('<i class="copy icon"></i> Copy');
            }, 2000);
        }
    }

    $('#copyKeyBtn').on('click', copyNewKey);
    $('#newKeyValue').on('click', copyNewKey);
});
