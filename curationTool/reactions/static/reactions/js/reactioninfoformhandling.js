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
            console.log('Comments Div clicked. infoType set to:', infoType);
            setupSubmitHandler("submitAddInfocom","infoTextInputcom");

        } else if (nameAttr === 'references-div') {
            infoType = 'Reference';
            const identifierType = document.getElementById('identifierType').value;

            extLinkType = '';  // Set external link type to an empty string
            refType = identifierType;      // Set reference type to an empty string
            console.log('References Div clicked. infoType set to:', infoType);
            setupSubmitHandler("submitAddInforef","infoTextInputref");


        } else if (nameAttr === 'extlinks-div') {
            const extType = document.getElementById('extType').value;

            infoType = 'External Link';
            extLinkType = extType;  // Set external link type to the selected value
            refType = '';      // Set reference type to an empty string
            console.log('External Links Div clicked. infoType set to:', infoType);
            setupSubmitHandler("submitAddInfoext","infoTextInputext");

        } else {
            infoType = 'Gene Info';
            extLinkType = '';  // Set external link type to an empty string
            refType = '';      // Set reference type to an empty string
    
            console.log('Gene Info Div clicked. infoType set to:', infoType);
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
        console.log(`${submitButtonId} clicked`);

        var userID = sessionStorage.getItem('userID'); 
        const urlParams = new URLSearchParams(window.location.search);
        var reactionId = urlParams.get('reaction_id');

        if (!userID) {
            alert("Please log in.");
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
            console.log('Info Text:', infoText);
            data.infoText = infoText;
            console.log('Data send:', data);
            submitData(data);
        } else {
            const geneInputs = document.getElementById('geneInfoInput');
            console.log('Gene inputs:', geneInputs);
            var geneinfo = geneInputs.textContent;
            console.log('Gene info:', geneinfo);

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
                console.log('Response:', response);
                if (response.processed_string) {
                    console.log('Processed string:', response.processed_string);
                    data.infoText = response.processed_string;
                    console.log('Gene Data:', data);
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
        }
    });
}


    function AddOrganLocation(data) {
        const url = '/gene_details_view/';
        console.log('Data gene', data);
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
            console.log('Results:', results);
            submitData(results);
        })
        .catch(error => {
            console.error('There was a problem with the fetch operation:', error);
        });
    }

    function submitData(data) {
        var infotype = data.infoType;
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
                console.log("obj", obj.body.message);
                console.log('Data:', data.reactionId);
                
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
                    console.log('Details:', details);
    
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
        console.log('Displaying tab content:', tableId, reaction_info, reaction_id);
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
            console.log('Item:', item);
        
            if (tableId === 'gene-info-content') {
                const parsedInfo = parseInfo(item.info);
                console.log('Parsed Info:', parsedInfo);
                const gprField = document.createElement('span');
                gprField.textContent = parsedInfo.gpr;
                row.appendChild(gprField);
        
                const geneField = document.createElement('span');
                geneField.textContent = parsedInfo.organ.join(', ');
                row.appendChild(geneField);
        
                const organField = document.createElement('span');
                row.appendChild(organField);
        
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
                const typeField = document.createElement('span');
                typeField.textContent = item.ref_type;
                row.appendChild(typeField);
        
                const identifierField = document.createElement('span');
                identifierField.textContent = item.info;
                row.appendChild(identifierField);
        
            } else if (tableId === 'ext-links-content') {
                const typeField = document.createElement('span');
                typeField.textContent = item.ext_link_type;
                console.log('Type:', item.ext_link_type);
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
                console.log('Delete button clicked',item.reactionId);
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