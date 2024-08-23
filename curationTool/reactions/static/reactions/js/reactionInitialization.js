// // Executes on window load, fetches reaction details based on the 'reaction_id' URL parameter, and displays them.
// window.onload = function() {
//     var errorMessageContainer = document.getElementById('error-message');
//     if (sessionStorage.getItem('userID') !== null) {
//         errorMessageContainer.innerHTML = 'Logging in...' + sessionStorage.getItem('userName');
//         username = sessionStorage.getItem('userName');
//         userID = sessionStorage.getItem('userID');
//         document.getElementById('userDisplay').textContent = `User: ${username}`;
//         document.getElementById('loginButton').textContent = 'Log out';
//         fetch(setSessionUser, {
//             method: 'POST',
//             headers: {
//                 'X-Requested-With': 'XMLHttpRequest',
//                 'X-CSRFToken': csrfToken
//             },
//             body: JSON.stringify({ 'userID': userID })
//         })
//         .then(response => response.json())
//         .then(data => {
//             if (data.status === 'success') {
//                 errorMessageContainer.style.display = 'none';
//                 errorMessageContainer.innerHTML = '';
//             } else {
//                 errorMessageContainer.innerHTML = 'Error in setting session user: ' + data.message;
//                 window.scrollTo(0, 0);
//             }
//         })
//     }
//     const urlParams = new URLSearchParams(window.location.search);
//     const reactionId = urlParams.get('reaction_id');
//     console.log('reactionId:', reactionId);
//     console.log('usr',userID);
//     const action = urlParams.get('action');
//     console.log('action:', action);
//     if (reactionId) {
//         fetch(getReaction + reactionId)
//             .then(response => response.json())
//             .then(async reactionData => { 
//                 let reactionTitle = document.getElementById('reactionTitle');
//                 reactionTitle.textContent = `Reaction: ${reactionData.short_name}`;
//                 reactionTitle.style.display = 'block';
//                 await updateFormFields(reactionData); 
//                 console.log('Reaction Data:', reactionData);
//                 updateReactionData(reactionData,userID,reactionId);
//                 hideDoneFieldButtons();
//                 hideDoneFieldButtonall();
//                 if (action === 'edit'){
//                 console.log('Editing reaction');
//                 updateSubmitButtonText(action);}
//                 else{createnewreaction();}
//                 // confirmAll(reactionData); // Call confirmAll after the other functions are done

//             })
//             .catch(error => {
//                 console.error('Error fetching reaction data:', error);
//             });
//     }    
//     // Get subsystems from the VMH (Virtual Metabolic Human) database.
//     if (subsystemList.length === 0) {
//         errorMessageContainer.innerHTML = 'Fetching subsystems from VMH...';
//         errorMessageContainer.style.display = 'block';
//         submitBtn = document.getElementById('submitBtn');
//         submitBtn.disabled = true;
//         window.scrollTo(0, 0);
//         fetch(getVMHsubsystems, {
//             method: 'GET',
//             headers: {
//                 'X-Requested-With': 'XMLHttpRequest',
//                 'X-CSRFToken': csrfToken
//             }
//         })
//         .then(response => response.json())
//         .then(data => {
//             if (data.error) {
//                 errorMessageContainer.innerHTML = 'Error in retreiving subsystems from VMH: ' + data.message;
//                 window.scrollTo(0, 0);
//             }
//             else {
//                 subsystemList = data.subsystem_list; // Save the returned list
//                 errorMessageContainer.style.display = 'none';
//                 errorMessageContainer.innerHTML = '';
//             }
//             submitBtn.disabled = false;
//         })
//         .catch(error => console.error('Error:', error));
//     }
//     setupTooltips();
// };

// function setupTooltips() {
//     document.querySelectorAll('.info-symbol').forEach(item => {
//         item.addEventListener('mouseenter', function() {
//             const tooltipContent = this.getAttribute('data-tooltip-content');
//             const tooltip = document.createElement('div');
//             tooltip.className = 'tooltip';
//             tooltip.innerHTML = tooltipContent;
//             this.appendChild(tooltip);
//         });
//         item.addEventListener('mouseleave', function() {
//             this.removeChild(this.querySelector('.tooltip'));
//         });
//     });
// }

// function hideDoneFieldButtons() {
//     console.log('Hiding Done Field Buttons');
//     const buttons = document.querySelectorAll('.done-field-btn');
//     buttons.forEach(button => {
//         button.style.display = 'none';
//     });

// }

// function hideDoneFieldButtonall(){
//     const buttonsall = document.querySelectorAll('.done-field-btn-all');
//     buttonsall.forEach(buttonsall => {
//         buttonsall.style.display = 'none';  // Hide the button
//     });
// }
// // Ensure this function is called after fetching reaction data
// fetch(getReaction + reactionId)
//     .then(response => response.json())
//     .then(async reactionData => { 
//         let reactionTitle = document.getElementById('reactionTitle');
//         reactionTitle.textContent = `Reaction: ${reactionData.short_name}`;
//         reactionTitle.style.display = 'block';
//         await updateFormFields(reactionData); 
//         displayReactionDetails(reactionData);

//         // Confirm All after the other functions are done
//         // confirmAll(reactionData); 
//     })
//     .catch(error => {
//         console.error('Error fetching reaction data:', error);
//     });

//     function createnewreaction() {
//         // Get the input element by its ID
//         const submitBtn = document.getElementById('submitBtn');
        
//         // Change the value of the input element
//         submitBtn.value = 'Create New Reaction';
        
//         // Add an event listener to the input element to handle the click event
//         submitBtn.addEventListener('click', function(event) {
//             // Prevent the default form submission behavior
//             event.preventDefault();
            
//             // Redirect to the homepage
//             window.location.href = window.location.origin;
//         });
//     }
    


//     function updateSubmitButtonText(action) {
//         const submitBtn = document.getElementById('submitBtn');
//         if (action === 'edit') {
//             submitBtn.value = 'Edit Reaction';
//         } else {
//             submitBtn.value = 'Create New Reaction';
//         }
//     }
    

//     // function updateReactionData(reactionData, userID, reactionId) {
//     //     // Send reaction data to the Python view function
//     //     fetch('update_reaction_data', {
//     //         method: 'POST',
//     //         headers: {
//     //             'Content-Type': 'application/json',
//     //             'X-Requested-With': 'XMLHttpRequest',
//     //             'X-CSRFToken': csrfToken // Add CSRF token if required
//     //         },
//     //         body: JSON.stringify({
//     //             reaction_id: reactionId,
//     //             user_id: userID,
//     //             reaction_data: reactionData
//     //         })
//     //     })
//     //     .then(response => response.json())
//     //     .then(updatedData => {
//     //         // Call displayReactionDetails with the updated data
//     //         displayReactionDetails(updatedData);
//     //     })
//     //     .catch(error => {
//     //         console.error('Error updating reaction data:', error);
//     //     });
//     // }
    