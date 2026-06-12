document.addEventListener('DOMContentLoaded', function () {
  const shareTemplatesContainer = document.getElementById('shareTemplatesContainer');
  const manageTemplateBtn = document.getElementById('manageTemplateBtn');
  const shareTemplateBtn = document.getElementById('shareTemplateBtn');

  // Open Manage Templates Modal
  manageTemplateBtn.addEventListener('click', function () {
    userID = sessionStorage.getItem('userID');
    if (!userID) {
        Notify.error('Please log in to manage templates.');
        return;
    }
    $('#manageTemplatesModal').modal('show');
  });

  // Open Share Templates Modal when Share button is clicked
  shareTemplateBtn.addEventListener('click', function () {
      $('#manageTemplatesModal').modal('hide'); // Close Manage Templates Modal
      $('#shareTemplateModal').modal('show');  // Open Share Templates Modal
  });

  const confirmShareBtn = document.getElementById('confirmShare');
  let templateList = [];
  let templatesFetched = { value: false };

  // Fetch templates and populate buttons
  shareTemplateBtn.addEventListener('click', async function () {
      $('#shareTemplateModal').modal('show');
      shareTemplatesContainer.innerHTML = ''; // Clear existing buttons

      // Fetch templates (reuse fetchTemplates if available)
      await window.ReactionUtils.fetchTemplates(templateList, templatesFetched);

      // Create buttons for each template
      templateList.forEach(templateName => {
          const button = document.createElement('div');
          button.textContent = templateName;
          button.className = 'template-button';
          button.addEventListener('click', function () {
              button.classList.toggle('selected'); // Toggle selection
          });
          shareTemplatesContainer.appendChild(button);
      });
  });

  // Handle share action
  confirmShareBtn.addEventListener('click', async function () {
      const selectedTemplates = [...shareTemplatesContainer.getElementsByClassName('selected')].map(
          button => button.textContent
      );
      const shareWithUser = document.getElementById('shareWithUser').value.trim();
      const userID = sessionStorage.getItem('userID');

      if (!userID) {
          Notify.error('Please log in before sharing templates.');
          return;
      }
      if (!selectedTemplates.length) {
          Notify.error('Please select at least one template to share.');
          return;
      }
      if (!shareWithUser) {
          Notify.error('Please enter the username of the user to share with.');
          return;
      }

      try {
          const response = await fetch('/share_template/', {
              method: 'POST',
              headers: {
                  'Content-Type': 'application/json',
                  'X-Requested-With': 'XMLHttpRequest',
                  'X-CSRFToken': csrfToken,
              },
              body: JSON.stringify({
                  userID: userID,
                  share_with_user: shareWithUser,
                  template_names: selectedTemplates,
              }),
          });

          const data = await response.json();
          if (response.ok && data.status === 'success') {
              Notify.success('Templates shared with ' + shareWithUser + ' successfully');
              $('#shareTemplateModal').modal('hide');
          } else {
              Notify.error(`Error sharing templates: ${data.message || 'Unknown error'}`);
          }
      } catch (error) {
          console.error('Error:', error);
          Notify.error('An error occurred while sharing the templates.');
      }
  });
});

document.addEventListener('DOMContentLoaded', function () {
    const editTemplateBtn = document.getElementById('editTemplateBtn');
    const editTemplatesContainer = document.getElementById('editTemplatesContainer');
    const deleteTemplateBtn = document.getElementById('deleteTemplateBtn');
    const submitTemplateChangesBtn = document.getElementById('submitTemplateChanges');

    let userID = sessionStorage.getItem('userID');
    let selectedTemplate = null;
    let templateList = [];
    let templatesFetched = { value: false };
    function resetTemplateList() {
        templateList = [];
        templatesFetched.value = false;
    }
    let nameInput = document.getElementById('editTemplateName');
    let descInput = document.getElementById('editTemplateDescription');
    // Open Edit Modal
    editTemplateBtn.addEventListener('click', async function () {
        $('#editTemplateModal').modal('show');
        await loadTemplates();
    });

    async function loadTemplates() {
        console.log('Loading templates...');
        editTemplatesContainer.innerHTML = '';
        try {
            await window.ReactionUtils.fetchTemplates(templateList, templatesFetched);
            templateList.forEach(templateName => {
                const button = document.createElement('div');
                button.textContent = templateName;
                button.className = 'template-button';
                button.addEventListener('click', function () {
                    // Highlight the selected template
                    document.querySelectorAll('.template-button').forEach(btn => btn.classList.remove('selected'));
                    button.classList.add('selected');
                    selectedTemplate = templateName; // Update selected template
                    loadTemplateDetails(userID, templateName);
                });
                editTemplatesContainer.appendChild(button);
            });
        } catch (error) {
            console.error('Error loading templates:', error);
        }
    }

    async function loadTemplateDetails(userID, templateName) {
        try {
            const response = await fetch('/get_rxn_template/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                },
                body: JSON.stringify({
                    reaction_type: templateName,
                    userID: userID
                })
            });
            
            const details = await response.json();
            nameInput.value = details.name;
            descInput.value = details.description || '';
        } catch (error) {
            console.error('Error loading template details:', error);
        }
    }
    // Save Changes
    submitTemplateChangesBtn.addEventListener('click', async function () {
        if (!selectedTemplate) {
            Notify.error('Please select a template to edit.');
            return;
        }

        const newName = nameInput.value.trim();
        const newDesc = descInput.value.trim();

        if (!newName) {
            Notify.error('Template name is required.');
            return;
        }

        try {
            const response = await fetch('/update_template/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                },
                body: JSON.stringify({
                    old_name: selectedTemplate,
                    new_name: newName,
                    description: newDesc,
                    userID: sessionStorage.getItem('userID')
                })
            });

            const data = await response.json();
            if (response.ok && data.status === 'success') {
                Notify.success('Template updated successfully!');
                resetTemplateList();
                loadTemplates();
                nameInput.value = '';
                descInput.value = '';
            } else {
                Notify.error(`Error: ${data.message || 'Failed to update template'}`);
            }
        } catch (error) {
            console.error('Error updating template:', error);
            Notify.error('An error occurred while updating the template.');
        }
    });
    // Delete Template
    deleteTemplateBtn.addEventListener('click', async function () {
        if (!selectedTemplate) {
            Notify.error('Please select a template to delete.');
            return;
        }

        const confirmed = await Notify.confirm({
            title: 'Delete template',
            message: `Are you sure you want to delete the template "${selectedTemplate}"?`,
            confirmText: 'Delete',
            cancelText: 'Cancel',
            danger: true,
        });
        if (!confirmed) return;

        try {
            const response = await fetch('/delete_template/', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                },
                body: JSON.stringify({
                    template_name: selectedTemplate,
                    userID: userID,
                }),
            });

            const data = await response.json();
            if (response.ok && data.status === 'success') {
                Notify.success('Template deleted successfully!');
                resetTemplateList();
                loadTemplates();
            } else {
                Notify.error(`Error deleting template: ${data.message || 'Unknown error'}`);
            }
        } catch (error) {
            console.error('Error:', error);
            Notify.error('An error occurred while deleting the template.');
        }
    });
});
