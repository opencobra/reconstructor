// Variable to hold the type of info based on the clicked div
var infoType = '';
var extLinkType = '';
var refType = '';

// Event listener for divs with class 'info-div'
document.querySelectorAll('.content-div').forEach(function(div) {
    div.addEventListener('click', function() {

        const nameAttr = this.getAttribute('name');
        
        if (nameAttr === 'comments-div') {
            infoType = 'Comment';
            extLinkType = '';  // Set external link type to an empty string
            refType = '';      // Set reference type to an empty string
            setupSubmitHandler("submitAddInfocom","infoTextInputcom");

        } else if (nameAttr === 'references-div') {
            infoType = 'Reference';
            const identifierType = document.getElementById('identifierType').value;

            extLinkType = '';  // Set external link type to an empty string
            refType = identifierType;      // Set reference type to an empty string
            setupSubmitHandler("submitAddInforef","infoTextInputref");
            // update placeholder text based on the selected reference type
            if (identifierType.trim() !== '') {
                document.getElementById('infoTextInputref').placeholder = `Enter ${identifierType}`;
            }
            else {
                document.getElementById('infoTextInputref').placeholder = 'Select PMID or DOI from dropdown';
            }

        } else if (nameAttr === 'extlinks-div') {
            const extType = document.getElementById('extType').value;

            infoType = 'External Link';
            extLinkType = extType;  // Set external link type to the selected value
            refType = '';      // Set reference type to an empty string
            setupSubmitHandler("submitAddInfoext","infoTextInputext");
            // update placeholder text based on the selected external link type (only if not empty)
            if (extType.trim() !== '') {
                document.getElementById('infoTextInputext').placeholder = `Enter ${extType} ID`;
            }
            else {
                document.getElementById('infoTextInputext').placeholder = 'Select External link type from dropdown';
            }

        } else {
            infoType = 'Gene Info';
            extLinkType = '';  // Set external link type to an empty string
            refType = '';      // Set reference type to an empty string
    
            setupSubmitHandler("submitAddInfogene","geneInfoInput");

        }
    });
});



function setupSubmitHandler(submitButtonId, Infotextid) {

    const submitButton = document.getElementById(submitButtonId);
    
    // Clone the submit button to remove all existing event listeners
    const newSubmitButton = submitButton.cloneNode(true);
    submitButton.parentNode.replaceChild(newSubmitButton, submitButton);

    // Add the new event listener
    newSubmitButton.addEventListener('click', function() {

        var userID = sessionStorage.getItem('userID'); 
        const urlParams = new URLSearchParams(window.location.search);
        var reactionId = urlParams.get('reaction_id');

        if (!userID) {
            alert("Please log in.");
            return;
        }
        if (!reactionId && infoType !== 'Gene Info') {
            alert("Please create and save a reaction first.");
            return;
        }
        // Prepare the data to be submitted
        var data = {
            userID: userID,
            infoType: infoType,
            extLinkType: extLinkType,
            refType: refType,
            reactionId: reactionId
        };

        if (infoType !== 'Gene Info') {
            var infoText = document.getElementById(Infotextid).value;
            data.infoText = infoText;
            submitData(data);
        } else {
            const geneInputs = document.getElementById('geneInfoInput');
            var geneinfo = geneInputs.textContent;

            fetch('gene_parsing/', {
                method: 'POST',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken,
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({ geneinfo })
            })
            .then(response => {
                if (!response.ok) {
                    throw new Error(`HTTP error! Status: ${response.status}`);
                }
                return response.json();
            })
            .then(response => {
                if (response.processed_string) {
                    data.infoText = response.processed_string;
                    console.log(data);
                    AddOrganLocation(data);
                } else {
                    if (response.error) {
                        alert(`Error: ${response.error}`);
                    } else {
                        alert('Gene info not added');
                    }            
                }
            })
            .catch((error) => {
                console.error('Error:', error);
                alert('Gene info not added');
            });
            // fetch gene predictions
            fetch(getAIresponse, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken  // Ensure CSRF token is included if needed
                },
                body: JSON.stringify({
                    key: geneinfo,  // Gene data
                    temperature: 0.5  // Temperature value
                })
            }).then(response => {
                response.json().then(response => {
                    if (response.status != 'success') {
                        alert(response.error_message);
                    }
                    else {
                        predictions = response.predictions;
                        displayGenePredictions(predictions);
                    }
                });
            }).catch(error => {
                console.error('Error fetching gene predictions:', error);
            });
        }
    });
}
function displayGenePredictions(predictions) {
    // Get the gene info container div
    const geneInfoContainer = document.getElementsByName('gene_info-div')[0];

    // Create a line separator
    const separator = document.createElement('hr');
    separator.classList.add('separator-line');
    geneInfoContainer.appendChild(separator);

    // Add global title for the predictions section
    const title = document.createElement('h3');
    title.textContent = 'GPT Predicted Reactions';
    title.classList.add('predictions-title');
    geneInfoContainer.appendChild(title);

    // Iterate through each gene in predictions and create toggles
    for (const [gene, reactions] of Object.entries(predictions)) {
        const geneToggleContainer = document.createElement('div');
        geneToggleContainer.classList.add('gene-toggle-container');

        // Create a container for the gene title and arrow toggle
        const geneTitleContainer = document.createElement('div');
        geneTitleContainer.classList.add('gene-title-container');

        // Add a title for each gene
        const geneTitle = document.createElement('h4');
        geneTitle.textContent = `GPT predicted reactions for Gene: ${gene}`;
        geneTitle.classList.add('gene-title');
        geneTitleContainer.appendChild(geneTitle);

        // Create an arrow div that will toggle the visibility of the predictions
        const arrowDown = document.createElement('div');
        arrowDown.classList.add('arrow-down');
        geneTitleContainer.appendChild(arrowDown);

        // Append the gene title container to the gene toggle container
        geneToggleContainer.appendChild(geneTitleContainer);

        // Check if reactions is a string, parse it if necessary
        let reactionListArray = reactions;
        if (typeof reactions === 'string') {
            try {
                reactionListArray = JSON.parse(reactions); // Parse if it's a JSON string
            } catch (error) {
                console.error("Failed to parse reactions: ", error);
                continue;  // Skip this gene if parsing fails
            }
        }

        // Create a div to hold the list of reactions, hidden by default
        const reactionList = document.createElement('ol');  // Ordered list
        reactionList.classList.add('reaction-list');
        reactionList.style.display = 'none';  // Initially hidden

        // Check if the reaction list is empty
        if (reactionListArray.length === 0) {
            const noReactionsMessage = document.createElement('p');
            noReactionsMessage.textContent = `No GPT predicted reactions for Gene: ${gene}`;
            noReactionsMessage.classList.add('no-reactions-message');
            reactionList.appendChild(noReactionsMessage);
        } else {
            // Iterate through the list of reactions for each gene
            reactionListArray.forEach((reaction) => {
                const listItem = document.createElement('li');
                listItem.textContent = reaction;
                reactionList.appendChild(listItem);
            });
        }

        // Append the reaction list to the gene toggle container
        geneToggleContainer.appendChild(reactionList);

        // Add toggle functionality to the arrow
        arrowDown.addEventListener('click', function () {
            this.classList.toggle('active');
            if (reactionList.style.display === 'none') {
                reactionList.style.display = 'block';
            } else {
                reactionList.style.display = 'none';
            }
        });

        // Append the toggle container to the gene info container
        geneInfoContainer.appendChild(geneToggleContainer);
    }
}

    function AddOrganLocation(data) {
        const url = '/gene_details_view/';
        fetch(url, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'Content-Type': 'application/json',
                'X-CSRFToken': 'csrftoken'
            },
            body: JSON.stringify({ infoText: data })
        })
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.json();
        })
        .then(data => {
            if (data.error) {
                console.error('Error from Django:', data.error);
                return;
            }

            const results = data;
            submitData(results);
        })
        .catch(error => {
            console.error('There was a problem with the fetch operation:', error);
        });
    }

    function submitData(data) {
        var infotype = data.infoType;
        if (infotype === 'Reference' || infotype === 'External Link') {
            // assert that both extLinkType and refType cannot be empty string (strip whitespace)
            if (data.extLinkType.trim() === '' && data.refType.trim() === '') {
                alert('Select type from dropdown');
                return;
            } 
        }
        if (data.infoType === 'Reference') {
            foundData = fetchPubMedInfo(data.infoText);
            foundData = foundData.then(function(result) {
                let refFoundData = result.status;
                let refTitle = '';
                let refAbstract = '';
                let refAuthor = '';
                let refFoundMessage = '';
                if (refFoundData === 'success') {
                    refTitle = result.title;
                    refAbstract = result.abstract;
                    refAuthor = result.author;
                }
                else {
                    refFoundMessage = result.message;
                }
            });
        }
        fetch(addInfo2Reaction, {
            method: 'POST',
            body: JSON.stringify(data),
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
                'Content-Type': 'application/json'
            }
        })
        .then(response => response.json().then(data => ({ status: response.status, body: data })))
        .then(obj => {
            if (obj.status === 200) {

                // Refactored to use POST request instead of GET
                fetch(getReactionDetails, {
                    method: 'POST',
                    body: JSON.stringify(data.reactionId),  // Corrected to pass reactionId as an object
                    headers: {
                        'X-Requested-With': 'XMLHttpRequest',
                        'X-CSRFToken': csrfToken,
                        'Content-Type': 'application/json'
                    }
                })
                .then(response => response.json())
                .then(details => {
                    // Conditionally display tab content based on infotype
                    switch (infotype) {
                        case 'Comment':
                            displayTabContent('comments-content', details.comments, data.reactionId);
                            break;
                        case 'Reference':
                            displayTabContent('refs-content', details.references, data.reactionId);
                            break;
                        case 'External Link':
                            displayTabContent('ext-links-content', details.external_links, data.reactionId);
                            break;
                        case 'Gene Info':
                            displayTabContent('gene-info-content', details.gene_info, data.reactionId);
                            break;
                        default:
                            console.error('Unknown infoType:', infotype);
                            break;
                    }
                })
                .catch(error => {
                    console.error('Error fetching reaction details:', error);
                });
    
            } else {
                console.error('Error:', obj.body);
                alert(obj.body.message);
            }
        })

    }
    
    function parseInfo(info) {
        const gprMatch = info.match(/GPR:\s*([^;]+)/);
        const organMatch = info.match(/ORGAN\((.*?)\)/);
        const subcellularMatch = info.match(/SUBCELLULAR\((.*?)\)/);
        
        let gpr = '';
        let organ = [];
        let subcellular = [];
    
        if (gprMatch || organMatch || subcellularMatch) {
            gpr = gprMatch ? gprMatch[1] : '';
            organ = organMatch ? organMatch[1].replace(/_/g, '').split(', ') : [];
            subcellular = subcellularMatch ? subcellularMatch[1].replace(/\[|\]/g, '').split(', ') : [];
        } else {
            // If no GPR, ORGAN, or SUBCELLULAR match, consider the whole info as gpr
            gpr = info;
            organ = ["-"];
            subcellular = ["-"];
        }
    
        return { gpr, organ, subcellular };
    }
    
    
    async function displayTabContent(tableId, reaction_info, reaction_id) {
        const tableElement = document.getElementById(tableId);
        if (!tableElement) {
            console.error('Table not found');
            return;
        }
        // Check if reaction_id is empty, if so, fetch session data
        if (!reaction_id && tableId === 'gene-info-content') {
            try {
                const response = await fetch('/check-session/', {
                    method: 'GET',
                    headers: {
                        'X-CSRFToken': csrfToken,
                        'Content-Type': 'application/json',
                    }
                });
    
                const data = await response.json();
                if (data.status === 'success' && data.gene_info) {
                    reaction_info = data.gene_info; // Use session data
                } else {
                    console.error('No gene info in session or an error occurred.');
                    return;
                }
            } catch (error) {
                console.error('Error fetching session data:', error);
                return;
            }
        }
        // Clear Existing Content
        tableElement.innerHTML = '';
    
        // Create Header Row
        const headerRow = document.createElement('div');
        headerRow.classList.add('header-row');
    
        if (tableId === 'gene-info-content') {
            // Headers for Gene Info content
            const gprHeader = document.createElement('span');
            gprHeader.textContent = 'GPR';
            headerRow.appendChild(gprHeader);
    
            const geneHeader = document.createElement('span');
            geneHeader.textContent = 'Organ';
            headerRow.appendChild(geneHeader);
    
            const organHeader = document.createElement('span');
            organHeader.textContent = 'Subcellular Location';
            headerRow.appendChild(organHeader);
    
            const userHeader = document.createElement('span');
            userHeader.textContent = 'Added by';
            headerRow.appendChild(userHeader);
        } else if (tableId === 'refs-content' || tableId === 'ext-links-content') {
            const typeHeader = document.createElement('span');
            typeHeader.textContent = 'Type';
            headerRow.appendChild(typeHeader);
    
            const identifierHeader = document.createElement('span');
            identifierHeader.textContent = 'Identifier';
            headerRow.appendChild(identifierHeader);
    
            const userHeader = document.createElement('span');
            userHeader.textContent = 'Added by';
            headerRow.appendChild(userHeader);
        } else {
            const commentHeader = document.createElement('span');
            commentHeader.textContent = 'Comment';
            headerRow.appendChild(commentHeader);
    
            const userHeader = document.createElement('span');
            userHeader.textContent = 'Added by';
            headerRow.appendChild(userHeader);
        }
    
        const emptyHeader = document.createElement('span');
        headerRow.appendChild(emptyHeader); // Empty cell for delete button
    
        tableElement.appendChild(headerRow);
    
        for (const item of reaction_info) {
            const row = document.createElement('div');
            row.classList.add('row');
        
            if (tableId === 'gene-info-content') {
                const parsedInfo = parseInfo(item.info);
                const gprField = document.createElement('span');
                gprField.textContent = parsedInfo.gpr;
                row.appendChild(gprField);
        
                const geneField = document.createElement('span');
                geneField.classList.add('gene-field-container');
                row.appendChild(geneField);

                const organField = document.createElement('span');
                row.appendChild(organField);
        
                // // Create buttons for organs
                parsedInfo.organ.forEach((organ) => {
                    const organButton = document.createElement('button');
                    organButton.textContent = organ;
                    organButton.classList.add('organ-btn');
                    organButton.addEventListener('click', () => {
                        applyOrganFunction(organ); // Define what happens when an organ button is clicked
                    });
                    geneField.appendChild(organButton);
                });

                // Create buttons for subcellular locations
                parsedInfo.subcellular.forEach((location) => {
                    const locationButton = document.createElement('button');
                    locationButton.textContent = location;
                    locationButton.classList.add('subcellular-location-btn');
                    locationButton.addEventListener('click', () => {
                        applySubcellularLocation(location);
                    });
                    organField.appendChild(locationButton);
                });
        
            } else if (tableId === 'refs-content') {
                createTooltip(item, row);
            }            
             else if (tableId === 'ext-links-content') {
                const typeField = document.createElement('span');
                typeField.textContent = item.ext_link_type;
                row.appendChild(typeField);
        
                const identifierField = document.createElement('span');
                identifierField.textContent = item.info;
                row.appendChild(identifierField);
        
            } else {
                const commentField = document.createElement('span');
                commentField.textContent = item.info;
                row.appendChild(commentField);
            }
        
            // Add userField and deleteButton for every row
            const userField = document.createElement('span');
            userField.textContent = item.user_name;
            row.appendChild(userField);
            if (reaction_id !== null) {
            const deleteButton = document.createElement('button');
            deleteButton.innerHTML = '<i class="fas fa-trash-alt"></i>';
            deleteButton.onclick = async () => {
                const reactionId = reaction_id;  // Ensure reaction_id is present in your item data
                const tabId = tableId;  // Use the current tabId context
                const itemToDelete = {
                    info: item.info,
                    ref_type: item.ref_type,  // Add only if applicable, based on tabId
                    ext_link_type: item.ext_link_type  // Add only if applicable, based on tabId
                };
        
                await deleteItem(reactionId, tabId, itemToDelete, row);
            };
            row.appendChild(deleteButton);
        
        }
            else{
                const deleteButtonGene = document.createElement('button');
                deleteButtonGene.innerHTML = '<i class="fas fa-trash-alt"></i>';    
                deleteButtonGene.onclick = async () => {
                    const success = await deleteGeneInfoFromSession(item.info);
                    if (success) {
                        row.remove();
                    }
                };  
                row.appendChild(deleteButtonGene);
    
            }


                
            tableElement.appendChild(row);
        }
        
        // Setup tooltips after rows are created
        setupTooltips();
        
    
    }
    
    // Function to apply subcellular location to the relevant fields
    function applySubcellularLocation(location) {
        // Apply the location to all substrate compartments
        const substrateCompsSelects = document.querySelectorAll('#substratesDiv .inputs-group select[name="subs_comps"]');
        substrateCompsSelects.forEach((select) => {
            if (!select.disabled) {
                select.value = location;
            }
        });
    
        // Apply the location to all product compartments
        const productCompsSelects = document.querySelectorAll('#productsDiv .inputs-group select[name="prod_comps"]');
        productCompsSelects.forEach((select) => {
            if (!select.disabled) {
                select.value = location;
            }
        });
    }
    

    async function deleteGeneInfoFromSession(itemToDelete) {
        const confirmation = confirm("Are you sure you want to delete this gene info from the session?");
        if (!confirmation) {
            return;
        }
    
        try {
            const response = await fetch('/delete-gene-info/', {
                method: 'POST',
                headers: {
                    'X-CSRFToken': csrfToken,
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    info_to_delete: itemToDelete,
                }),
            });
    
            const data = await response.json();
            if (data.status === 'success') {
                console.log('Gene info deleted successfully from session');
                return true;  // Indicates successful deletion
            } else {
                console.error('Error deleting gene info from session:', data.message);
            }
        } catch (error) {
            console.error('Error:', error);
        }
        return false;  // Indicates unsuccessful deletion
    }
    

    async function deleteItem(reactionId, tabId, itemToDelete, rowElement) {
        try {
            const response = await fetch('/delete_reaction_info/', {  // URL from your Django path configuration
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'X-CSRFToken': csrfToken
                },
                body: JSON.stringify({
                    reaction_id: reactionId,
                    tab_id: tabId,
                    item_to_delete: itemToDelete
                })
            });
    
            const data = await response.json();
    
            if (response.ok && data.status === 'success') {
                // Remove the row element from the DOM if deletion was successful
                rowElement.remove();
                alert(data.message);
            } else {
                // Handle errors returned by the server
                alert(`Error: ${data.message}`);
            }
        } catch (error) {
            console.error('Error:', error);
            alert('An unexpected error occurred. Please try again.');
        }
    }



    function applyOrganFunction(organ) {
        // Assuming `organField` is the container where organ tags should be added
        const organField = document.querySelector('.tags-input-container'); // Adjust this selector as needed
    
        // Debugging step: Log a warning if organField is not found
        if (!organField) {
            console.warn('organField element not found. Check the selector or HTML structure.');
            return; // Exit the function to avoid further errors
        }
    
        // Function to add a tag to the organField
        function addTag(tagText) {
            // Check if the tag already exists to prevent duplicates
            if ([...organField.children].some(tag => tag.textContent.replace('×', '').trim() === tagText)) return;
    
            // Create the tag element
            const tag = document.createElement('span');
            tag.className = 'tag';
            tag.contentEditable = false; // Make the tag non-editable
            tag.textContent = tagText;
    
            // Create a close button for the tag
            const closeBtn = document.createElement('span');
            closeBtn.className = 'close-btn';
            closeBtn.textContent = '×';
            closeBtn.addEventListener('click', function() {
                organField.removeChild(tag); // Remove the tag when the close button is clicked
            });
    
            // Append the close button to the tag
            tag.appendChild(closeBtn);
            organField.appendChild(tag); // Insert the tag inside the organ field
    
            // Add a space and an empty span to ensure the cursor moves outside the tag
            const spaceNode = document.createTextNode(' ');
            const caretSpan = document.createElement('span');
            caretSpan.innerHTML = '&nbsp;'; // Invisible space to ensure caret positioning
    
            organField.appendChild(spaceNode);
            organField.appendChild(caretSpan);
    
            // Move the cursor to the end of the contenteditable div
            const range = document.createRange();
            const selection = window.getSelection();
            range.setStartAfter(caretSpan);
            range.collapse(true);
            selection.removeAllRanges();
            selection.addRange(range);
    
            organField.focus(); // Refocus the input area
        }
    
        // Add the organ as a tag
        addTag(organ);
    }
    

    function displayreactioninfo(reactionData) {
        // Extract the reactionId from reactionData
        const reactionId = reactionData;
        var user_id = sessionStorage.getItem('userID');
        // Log the reactionId for debugging purposes
    
        if (reactionId === null) {
            // If reactionId is null, fetch session data
            fetch('/check-session/', {
                method: 'GET',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken, // Include CSRF token if needed
                }
            })
                .then(response => response.json())
                .then(data => {
                    if (data.status === 'success') {
                        // If session data is found, display the gene info
                        displayTabContent('gene-info-content', data.gene_info, null);
                    } else {
                        return;
                    }
                })
                .catch(error => {
                    console.error('Error fetching session data:', error);
                });
        } else {
            // If reactionId is not null, first fetch session data
            fetch('/check-session/', {
                method: 'GET',
                headers: {
                    'X-Requested-With': 'XMLHttpRequest',
                    'X-CSRFToken': csrfToken, // Include CSRF token if needed
                }
            })
            .then(sessionResponse => sessionResponse.json())
            .then(sessionData => {
            
                if (sessionData.status === 'success' && sessionData.gene_info) {
                    // Process each item in the session data
                    sessionData.gene_info.forEach(geneInfoItem => {
                        const dataToSubmit = {
                            userID: user_id,  // assuming you have this from your context
                            infoType: 'Gene Info',
                            infoText: geneInfoItem.info,
                            reactionId: reactionId.reaction_id  // Include the reaction ID
                        };
                        // Submit the data to the reaction
                        submitData(dataToSubmit);
                    });
            
                    // After all data has been submitted, clear the session
                    clearSession();
                } else {
                    console.log('No gene info in session or an error occurred.');
                }
            })
            .catch(error => {
                console.error('An error occurred:', error);
                })
                .catch(sessionError => {
                    console.error('Error handling session data:', sessionError);
                })
                .then(() => {
                    // After checking the session, proceed with fetching reaction details
                    fetch(getReactionDetails, {
                        method: 'POST',
                        body: JSON.stringify(reactionId.reaction_id),  // Assuming the POST requires a JSON body
                        headers: {
                            'X-Requested-With': 'XMLHttpRequest',
                            'X-CSRFToken': csrfToken,
                            'Content-Type': 'application/json'  // Ensuring the content type is JSON
                        }
                    })
                        .then(response => response.json())
                        .then(details => {
    
                            // Display content in the appropriate tabs
                            displayTabContent('refs-content', details.references, reactionId.reaction_id);
                            displayTabContent('ext-links-content', details.external_links, reactionId.reaction_id);
                            displayTabContent('gene-info-content', details.gene_info, reactionId.reaction_id);
                            displayTabContent('comments-content', details.comments, reactionId.reaction_id);
                        })
                        .catch(error => {
                            console.error('Error fetching reaction details:', error);
                        });
                });
        }
    }
    // Function to create tooltip content
function createTooltipContent(title, author, abstract) {
    const tooltipContent = document.createElement('div');
    tooltipContent.innerHTML = `
        <strong>Title:</strong> ${title}<br>
        <strong>Author:</strong> ${author}<br>
        <strong>Abstract:</strong> ${abstract}
    `;
    return tooltipContent;
}

// Updated function to handle tooltip creation and data fetching
function createTooltip(item, row) {
    const typeField = document.createElement('span');
    typeField.textContent = item.ref_type;
    row.appendChild(typeField);

    // Create a container for the identifier and tooltip
    const identifierField = document.createElement('span');
    identifierField.textContent = item.info;
    identifierField.classList.add('tooltip-container'); // Add tooltip container class
    // Determine which function to fetch data from
    const fetchPmidorDoi = item.ref_type === 'pmid' ? fetchPubMedInfo(item.info) : fetchDOIInfo(item.info);

    fetchPmidorDoi.then(result => {
        let refFoundData = result.status;
        let refTitle = result.title || 'N/A';
        let refAbstract = result.abstract || 'N/A';
        let refAuthor = result.author || 'N/A';

        if (refFoundData === 'success') {
            // Set tooltip content for found reference
            const tooltipElement = document.createElement('div');
            tooltipElement.classList.add('tooltip-text'); // Add tooltip text class

            const content = createTooltipContent(refTitle, refAuthor, refAbstract);
            tooltipElement.innerHTML = ''; // Clear existing text
            tooltipElement.appendChild(content);
            identifierField.style.color = 'green'; // Color green if found
            identifierField.appendChild(tooltipElement);
        } else {
            identifierField.style.color = 'red'; // Color red if not found
        }
    });

    // Append the tooltip element to the identifier field
    row.appendChild(identifierField);
}
