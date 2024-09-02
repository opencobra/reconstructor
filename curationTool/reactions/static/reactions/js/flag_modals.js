

document.addEventListener('DOMContentLoaded', function () {
    const createFlagModal = document.getElementById('createFlagModalCustom');
    const modalOverlay = document.getElementById('modalOverlayCustom');
    const closeSaveReaction = document.getElementById('closeSaveReactionModal');
    const saveReactionModal = document.getElementById('saveReactionModal');
    const createFlagButton = document.getElementById('createFlagButtonCustom');
    const closeCreateFlagModal = document.getElementById('closeCreateFlagModalCustom');

    function openCreateFlagModal() {
        createFlagModal.classList.add('active');
        modalOverlay.classList.add('active');
        saveReactionModal.style.display = 'none';
        document.getElementById('modalBackground').style.display = 'none';
    }

    function closeCreateFlagModalAndReopenSaveReaction() {
        createFlagModal.classList.remove('active');
        modalOverlay.classList.remove('active');
        saveReactionModal.style.display = 'block';
        document.getElementById('modalBackground').style.display = 'block';
    }

    createFlagButton.addEventListener('click', openCreateFlagModal);
    closeCreateFlagModal.addEventListener('click', closeCreateFlagModalAndReopenSaveReaction);
    modalOverlay.addEventListener('click', closeCreateFlagModalAndReopenSaveReaction);
    closeSaveReaction.addEventListener('click', function() {
        saveReactionModal.style.display = 'none';
        document.getElementById('modalBackground').style.display = 'none';
    });

});
