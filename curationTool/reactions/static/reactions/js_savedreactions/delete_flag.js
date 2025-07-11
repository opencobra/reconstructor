document.addEventListener('DOMContentLoaded', () => {
  const modal      = document.getElementById('deleteFlagModal');
  const msg        = document.getElementById('deleteFlagMessage');
  const closeBtn   = document.getElementById('closeDeleteFlagModal');
  const confirmBtn = document.getElementById('confirmDeleteFlagBtn');

  /* ---------- 1. delegated listener for every current/future flag ---------- */
  document.addEventListener('click', evt => {
    const icon = evt.target.closest('.flag-icon');
    if (!icon) return;                             // not a flag click

    modal.dataset.flagId      = icon.dataset.flagId;
    modal.dataset.flagName   = icon.dataset.flagName;
    modal.dataset.reactionId = icon.dataset.reactionId;

    msg.textContent =
      `Remove flag "${modal.dataset.flagName}" from this reaction?`;
    modal.style.display = 'block';
  });

  /* ---------- 2. close modal ---------- */
  closeBtn.addEventListener('click', () => {
    modal.style.display = 'none';
  });

  /* ---------- 3. confirm deletion ---------- */
  confirmBtn.addEventListener('click', () => {
    fetch('/saved_reactions/remove_flag/', {
      method : 'POST',
      headers: {
        'Content-Type': 'application/json',
        'X-CSRFToken' : csrfToken
      },
      body: JSON.stringify({
        user_id    : userID,
        reaction_id: modal.dataset.reactionId,
        flag_id    : modal.dataset.flagId
      })
    })
    .then(r => r.json())
    .then(data => {
      if (data.status === 'success') {
        location.reload();                 // or update the row in place
      } else {
        alert(data.message || 'Error removing flag.');
      }
    })
    .catch(err => {
      console.error(err);
      alert('Something went wrong.');
    });

    modal.style.display = 'none';
  });
});
