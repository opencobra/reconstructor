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
            } else {
                $('#keyError').text(data.message || 'Failed to generate key.').show();
            }
        } catch (error) {
            $('#keyError').text('Failed to generate key: ' + error).show();
        } finally {
            btn.removeClass('loading disabled');
        }
    });

    $('#copyKeyBtn').on('click', function () {
        const input = document.getElementById('newKeyValue');
        input.select();
        navigator.clipboard.writeText(input.value);
        $(this).find('i').removeClass('copy').addClass('check');
    });
});
