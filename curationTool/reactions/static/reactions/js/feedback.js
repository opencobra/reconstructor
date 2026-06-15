/**
 * feedback.js — opens the feedback modal and submits user feedback to the
 * backend, which stores it locally and mirrors it as a GitHub issue.
 */
document.addEventListener('DOMContentLoaded', function () {
    const feedbackNav = document.getElementById('feedback-nav');
    if (feedbackNav) {
        feedbackNav.addEventListener('click', function (e) {
            e.preventDefault();
            const msg = document.getElementById('feedbackMessage');
            if (msg) {
                msg.textContent = '';
            }
            $('#feedbackModal').modal('show');
        });
    }

    const feedbackForm = document.getElementById('feedbackForm');
    if (feedbackForm) {
        feedbackForm.addEventListener('submit', async function (e) {
            e.preventDefault();

            const typeEl = document.getElementById('feedbackType');
            const textEl = document.getElementById('feedbackText');
            const submitButton = document.getElementById('feedbackSubmitButton');
            const messageEl = document.getElementById('feedbackMessage');

            const text = (textEl.value || '').trim();
            if (!text) {
                if (window.Notify) {
                    Notify.warning('Please enter your feedback before submitting.');
                }
                return;
            }

            const userID = sessionStorage.getItem('userID');

            submitButton.classList.add('loading', 'disabled');
            if (messageEl) {
                messageEl.textContent = '';
            }

            try {
                const response = await fetch(submitFeedbackUrl, {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                        'X-CSRFToken': csrfToken
                    },
                    body: JSON.stringify({
                        feedback_type: typeEl.value,
                        text: text,
                        user_id: userID
                    })
                });
                const data = await response.json();

                if (data.status === 'success') {
                    if (window.Notify) {
                        Notify.success(data.message || 'Thank you for your feedback!');
                    }
                    feedbackForm.reset();
                    $('#feedbackModal').modal('hide');
                } else {
                    const errMsg = data.message || 'Failed to submit feedback.';
                    if (window.Notify) {
                        Notify.error(errMsg);
                    } else if (messageEl) {
                        messageEl.textContent = errMsg;
                    }
                }
            } catch (err) {
                if (window.Notify) {
                    Notify.error('Could not submit feedback. Please try again.');
                } else if (messageEl) {
                    messageEl.textContent = 'Could not submit feedback. Please try again.';
                }
            } finally {
                submitButton.classList.remove('loading', 'disabled');
            }
        });
    }
});
