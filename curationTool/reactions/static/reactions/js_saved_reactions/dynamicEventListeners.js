function attachDynamicEventListeners() {
	// Attach a single event listener to the parent container that holds all your dynamic content
	const modalList = document.getElementById('wsReactionDetails');

	modalList.addEventListener('click', function (e) {
		// Handle the "Add" button clicks using event delegation
		if (e.target && e.target.matches('.add-reference, .add-ext-link, .add-comment')) {
			const type = e.target.classList.contains('add-reference') ? 'reference' : e.target.classList.contains('add-ext-link') ? 'ext-link' : 'comment';
			const parentSection = e.target.parentNode;
			const reactionId = e.target.getAttribute('data-reaction-id');
			const newItem = document.createElement('div');
			newItem.className = `${type}-item`;
			const items = parentSection.querySelectorAll(`.${type}-item`);
			const newIndex = items.length; // Calculate new index based on existing items
			newItem.setAttribute('data-reaction-id', reactionId);
			newItem.setAttribute('data-index', newIndex);
			newItem.innerHTML = type === 'ext-link' ? createExtLinkSelect({}, reactionId, newIndex) : '';
			newItem.innerHTML += type === 'reference' ? createRefSelect({}, reactionId, newIndex) : '';
			newItem.innerHTML += `
                <input type="text" class="${type}-input" placeholder="Enter ${type}" value="" data-reaction-id="${reactionId}" data-index="${newIndex}">
                <button class="remove-${type}" data-reaction-id="${reactionId}" data-index="${newIndex}">Remove</button>
            `;
			parentSection.insertBefore(newItem, e.target);
		}

		// Handle the "Remove" button clicks using event delegation
		if (e.target && e.target.matches('.remove-reference, .remove-ext-link, .remove-comment,.remove-gene-info')) {
			e.target.parentElement.remove();
		}
	});
}
$(document).ready(function () {
	$('.view-description-btn').popup();
	
	// Initialize user dropdown
	const userDropdown = document.getElementById('userDropdown');
	const userDropdownTrigger = document.getElementById('userDisplayBtn');
	
	if (userDropdownTrigger && userDropdown) {
		userDropdownTrigger.addEventListener('click', function(e) {
			e.stopPropagation();
			userDropdown.classList.toggle('active');
		});
		
		// Close dropdown when clicking outside
		document.addEventListener('click', function(e) {
			if (!userDropdown.contains(e.target)) {
				userDropdown.classList.remove('active');
			}
		});
		
		// Handle dropdown item clicks
		const dropdownHome = document.getElementById('dropdown-home');
		if (dropdownHome) {
			dropdownHome.addEventListener('click', function() {
				window.location.href = '/';
				userDropdown.classList.remove('active');
			});
		}
		
		const dropdownSavedMetabolites = document.getElementById('dropdown-saved-metabolites');
		if (dropdownSavedMetabolites) {
			dropdownSavedMetabolites.addEventListener('click', function() {
				// Redirect to home page and trigger saved metabolites modal
				// Store flag in session storage to open modal after redirect
				sessionStorage.setItem('openSavedMetabolites', 'true');
				window.location.href = '/';
				userDropdown.classList.remove('active');
			});
		}
	}
});

document.addEventListener('DOMContentLoaded', function () {
	// Function to enable name editing
	function enableNameEditing(event) {
		let span = event.target;
		let input = span.nextElementSibling;
		if (input && input.classList.contains('reaction-name-input')) {
			span.style.display = 'none';
			input.style.display = 'inline-block';
			input.focus();
		}
	}

	// Function to save edited name
	async function saveEditedName(event) {
		let input = event.target;
		let newName = input.value.trim();
		let span = input.previousElementSibling;
		let reactionId = input.getAttribute('data-reaction-id');
		let userId = sessionStorage.getItem('userID');
		if (!newName || newName === span.textContent) {
			span.style.display = 'inline-block';
			input.style.display = 'none';
			return;
		}

		try {
			let response = await fetch(editReactionURL, {
				method: 'POST',
				headers: {
					'Content-Type': 'application/json',
					'X-CSRFToken': csrfToken,
				},
				body: JSON.stringify({
					reaction_id: reactionId,
					new_name: newName,
					user_id: userId,
				}),
			});

			let data = await response.json();

			if (data.success) {
				span.textContent = newName;
			} else {
				alert('Error updating reaction name: \n' + data.error);
				if (data.original_name) {
					input.value = data.original_name;
				}
			}
		} catch (error) {
			console.error('Failed to update reaction name:', error);
			if (data.original_name) {
				input.value = data.original_name;
			}
		}

		// Hide input and show updated name
		span.style.display = 'inline-block';
		input.style.display = 'none';
	}

	// Attach event listeners
	document.querySelectorAll('.reaction-name').forEach((span) => {
		span.addEventListener('click', enableNameEditing);
	});

	document.querySelectorAll('.reaction-name-input').forEach((input) => {
		// Save on pressing Enter
		input.addEventListener('keypress', function (event) {
			if (event.key === 'Enter') {
				saveEditedName(event);
			}
		});

		// Save when clicking outside input
		input.addEventListener('blur', saveEditedName);
	});
});

document.addEventListener('DOMContentLoaded', function () {
	document.querySelectorAll('.view-description-btn, .no-description').forEach((button) => {
		button.addEventListener('click', function () {
			const reactionId = button.getAttribute('data-reaction-id');
			const description = button.getAttribute('data-description') || '';

			// Create modal dynamically
			const modal = document.createElement('div');
			modal.className = 'ui modal';
			modal.id = 'descriptionModal';
			modal.innerHTML = `
                <div class="header">Reaction Description</div>
                <div class="content">
                    <p id="descriptionText">${description || 'No Description Available'}</p>
                </div>
                <div class="actions">
                    <button class="ui button" id="editDescriptionModalBtn" data-reaction-id="${reactionId}">Edit Description</button>
                    <button class="ui button" id="closeDescriptionModal">Close</button>
                </div>
            `;
			document.body.appendChild(modal);

			// Attach event listener to close the modal
			modal.querySelector('#closeDescriptionModal').addEventListener('click', function () {
				$(modal).modal('hide');
			});

			// Attach event listener to edit description
			modal.querySelector('#editDescriptionModalBtn').addEventListener('click', function () {
				const reactionId = this.getAttribute('data-reaction-id');
				const currentDescription = modal.querySelector('#descriptionText').textContent;
				const contentDiv = modal.querySelector('.content');

				// Replace content with an input field
				contentDiv.innerHTML = `
                    <textarea id="editDescriptionInput" style="width:100%;height:100px;">${
											currentDescription === 'No Description Available' ? '' : currentDescription
										}</textarea>
                    <button class="ui blue button" id="saveDescriptionBtn">Save</button>
                `;

				// Save button logic
				modal.querySelector('#saveDescriptionBtn').addEventListener('click', async function () {
					const newDescription = modal.querySelector('#editDescriptionInput').value.trim();
					if (!newDescription) return;

					try {
						const response = await fetch(editReactionURL, {
							method: 'POST',
							headers: {
								'Content-Type': 'application/json',
								'X-CSRFToken': csrfToken,
							},
							body: JSON.stringify({
								reaction_id: reactionId,
								new_description: newDescription,
							}),
						});

						const data = await response.json();
						if (data.success) {
							// Update UI dynamically
							const descriptionSpan = document.querySelector(`.no-description[data-reaction-id="${reactionId}"]`);
							if (descriptionSpan) {
								descriptionSpan.textContent = newDescription;
								descriptionSpan.classList.remove('no-description'); // Remove class since it's no longer empty
							}

							// Hide modal
							$(modal).modal('hide');
							// Refresh the page to show updated description
							window.location.reload();
						} else {
							alert('Error updating description.');
						}
					} catch (error) {
						console.error('Failed to update description:', error);
					}
				});
			});

			// Initialize and remove modal from DOM after closing
			$(modal)
				.modal({
					onHidden: function () {
						modal.remove();
					},
				})
				.modal('show');
		});
	});
});
