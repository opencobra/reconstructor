let eventListenersAttached = false; // Flag to track if event listeners have been attached

function createSectionHTML(sectionTitle, className, items, reactionId, isExtLink = false, isRef = false, isGene = false) {
    if (isGene) {
        return createGeneSectionHTML(sectionTitle, className, items, reactionId);
    }

    let sectionHTML = `
        <div class="${className}-section">
            <p>${sectionTitle}</p>
            ${items.map((item, index) => `
                <div class="${className}-item" data-reaction-id="${reactionId}" data-index="${index}"> <!-- Ensure this div represents an ext-link-item for external links -->
                    ${isExtLink ? createExtLinkSelect(item, reactionId, index) : ''}
                    ${isRef ? createRefSelect(item, reactionId, index) : ''}
                    <input type="text" class="${className}-input" placeholder="Enter ${sectionTitle.toLowerCase()}" value="${item.info}" data-reaction-id="${reactionId}" data-index="${index}">
                    <button class="remove-${className}" data-reaction-id="${reactionId}" data-index="${index}">Remove</button>
                </div>
            `).join('')}
            <button class="add-${className}" data-reaction-id="${reactionId}">Add ${sectionTitle.slice(0, -1)}</button>
        </div>
    `;
    return sectionHTML;
}
     
function createExtLinkSelect(item, reactionId, index) {
    const linkTypes = ['CHO Models', 'COG', 'EC Number', 'KEGG orthology', 'KEGG reaction', 'MetanetX', 'Rhea', 'SEED', 'Wikipedia'];
    let selectHTML = `
        <select class="ext-link-type-select" data-reaction-id="${reactionId}" data-index="${index}">
            ${linkTypes.map(type => `<option value="${type}" ${item.ext_link_type === type ? 'selected' : ''}>${type}</option>`).join('')}
        </select>
    `;
    return selectHTML;
}
function createRefSelect(item, reactionId, index) {
    const refTypes = ['DOI', 'PMID'];
    let selectHTML = `
        <select class="ref-type-select" data-reaction-id="${reactionId}" data-index="${index}">
            ${refTypes.map(type => `<option value="${type}" ${item.ref_type === type ? 'selected' : ''}>${type}</option>`).join('')}
        </select>
    `;
    return selectHTML;
}
// gene info , gene info, gene_info
// function createSectionHTMLGene(sectionTitle, className, items, reactionId, isExtLink = false, isRef = false, isGene = false) {
//     let sectionHTML = `
//         <div class="${className}-section">
//             <p>${sectionTitle}</p>
//             ${items.map((item, index) => `
//                 <div class="${className}-item">
//                     <input type="text" class="${className}-input" placeholder="Enter ${sectionTitle.toLowerCase()}" value="${item.info}" data-reaction-id="${reactionId}" data-index="${index}">
//                     <button class="remove-${className}" data-reaction-id="${reactionId}" data-index="${index}">Remove</button>
//                 </div>
//             `).join('')}
//         </div>
//     `;
//     return sectionHTML;
// }

function createGeneSectionHTML(sectionTitle, className, items, reactionId) {
    console.log(items);
    console.log(reactionId);
    let sectionHTML = `
        <div class="${className}-section">
            <p>${sectionTitle}</p>
            ${items.map((item, index) => `
                <div class="${className}-item" data-reaction-id="${reactionId}" data-index="${index}">
                <input type="text" class="${className}-input" placeholder="Enter ${sectionTitle.toLowerCase()}" value="${item.info}" data-reaction-id="${reactionId}" data-index="${index}" readonly>
                <button class="remove-${className}" data-reaction-id="${reactionId}" data-index="${index}">Remove</button>
                </div>
            `).join('')}
        </div>
    `;
    return sectionHTML;
}





async function callPrepareAddToVMH(reactionIds) {
    try {
        const response = await fetch(prepAddToVMHUrl, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'Content-Type': 'application/json',
                'X-CSRFToken': csrfToken,
            },
            body: JSON.stringify({
                reactionIds: reactionIds
            })
        });

        if (!response.ok) {
            // Handle HTTP errors
            console.error(`HTTP error! status: ${response.status}`);
            return null;
        }

        const data = await response.json();
        return data;
    } catch (error) {
        // Handle network or other errors
        console.error("Error fetching data: ", error);
        return null; // Or handle error appropriately
    }
}


function toggleDetails(event, headerElement) {
    // Check if the click target is one of the input fields or the select dropdown
    if (event.target.classList.contains('reaction-name-input') || 
        event.target.classList.contains('reaction-abbreviation-input') || 
        event.target.tagName.toLowerCase() === 'select' || 
        event.target.classList.contains("cs-info-button")
    ) {
        return; // Do nothing if clicked inside the input fields or select dropdown
    }

    // Find the toggle icon within the header element
    const toggleIcon = headerElement.querySelector('.toggle-icon');

    // Assume the detailed sections are initially visible; we need to check and set their initial state if undefined
    let initialStateSet = false; // Flag to track if we've set the initial state

    // The detailed sections are immediate siblings of the header, so we toggle their display.
    let currentElement = headerElement.nextElementSibling;
    while (currentElement) {
        // If initial state isn't set, explicitly set it to "none" to ensure consistent behavior
        if (!initialStateSet && (currentElement.style.display === "" || currentElement.style.display === "block")) {
            currentElement.style.display = "none";
            initialStateSet = true; // Mark that we've set the initial state
        } else {
            // Toggle the display style
            if (currentElement.style.display === "none") {
                currentElement.style.display = "block"; // Show the section
            } else {
                currentElement.style.display = "none"; // Hide the section
            }
        }
        // Move to the next sibling
        currentElement = currentElement.nextElementSibling;
    }

    // After setting the initial state or toggling, update the icon accordingly
    if (initialStateSet) {
        // If we just set the initial state, ensure the icon is set for the next action
        toggleIcon.textContent = "+"; // Set to plus since we've just hidden the sections
    } else {
        // Toggle the icon between "+" and "−"
        toggleIcon.textContent = toggleIcon.textContent === "+" ? "−" : "+";
    }
}

document.getElementById('addToVMH').addEventListener('click', async function() {
    var modal = document.getElementById('reactionModal');
    var modalList = document.getElementById('modalReactionsList');
    var modalTitle = document.getElementById('reactionModalTitle');
    modalList.innerHTML = ''; // Clear previous content

    if (checkedReactions.length === 1){
        modalTitle.textContent = "Adding the following reaction:";
    } else {
        modalTitle.textContent = "Adding the following reactions:";
    }

    if (checkedReactions.length === 0) {
        document.getElementById('alertMessage').textContent = "Please choose reactions to be added.";
        document.getElementById('alertModal').style.display = "block";
    } else {
        let vmhResponse;
        modal.style.display = "block"; // Show the modal
        document.getElementById('loadingIndicator').style.display = "flex"; // Show loading indicator
        document.getElementById('loadingText').textContent = 'Gathering Data and (if needed) Generating Abbreviations for Selected Reactions';

        try {
            vmhResponse = await callPrepareAddToVMH(checkedReactions);
        } catch (error) {
            console.error("Error fetching VMH preparation data: ", error);
            document.getElementById('alertMessage').textContent = "An error occurred while fetching data from VMH. Check console for error. Please try again later.";
            document.getElementById('alertModal').style.display = "block";
            modal.style.display = "none";
            return;
        }

        if (vmhResponse.status === 'error') {
            document.getElementById('alertMessage').textContent = vmhResponse.message;
            document.getElementById('alertModal').style.display = "block";
            modal.style.display = "none";
            return;
        }

        document.getElementById('loadingIndicator').style.display = "none";
        document.getElementById('loadingText').textContent = '';

        subsInVMH = vmhResponse.subs_in_vmh;
        prodsInVMH = vmhResponse.prods_in_vmh;
        subsAbbr = vmhResponse.subs_abbr;
        prodsAbbr = vmhResponse.prods_abbr;
        subsNeedNewNames = vmhResponse.subs_need_new_names;
        prodsNeedNewNames = vmhResponse.prods_need_new_names;
        reactionAbbrs = vmhResponse.reaction_abbrs;

        checkedReactions.forEach(function(reactionId) {
            var reactionIdNum = Number(reactionId);
            var reactionIndex = checkedReactions.indexOf(reactionId.toString());
            var reaction = reactions.find(r => r.pk === reactionIdNum);
            if (reaction) {
                let substrates_names = JSON.parse(reaction.fields.substrates_names);
                let products_names = JSON.parse(reaction.fields.products_names);
                let subs_comps = JSON.parse(reaction.fields.subs_comps);
                let prods_comps = JSON.parse(reaction.fields.prods_comps);
                let subs_stoich = JSON.parse(reaction.fields.subs_sch);
                let prods_stoich = JSON.parse(reaction.fields.prods_sch);
                let substrates = JSON.parse(reaction.fields.substrates);
                let products = JSON.parse(reaction.fields.products);
                let substrates_types = JSON.parse(reaction.fields.substrates_types);
                let products_types = JSON.parse(reaction.fields.products_types);
                let subsInVMHForReaction = subsInVMH[reactionIndex];
                let prodsInVMHForReaction = prodsInVMH[reactionIndex];
                let subsAbbrForReaction = subsAbbr[reactionIndex];
                let prodsAbbrForReaction = prodsAbbr[reactionIndex];
                let subsNeedNewNamesForReaction = subsNeedNewNames[reactionIndex];
                let prodsNeedNewNamesForReaction = prodsNeedNewNames[reactionIndex];
                let reactionAbbrForReaction = reactionAbbrs[reactionIndex];

                var listItem = document.createElement('div');
                listItem.className = "modal-reaction-entry";
                listItem.innerHTML = `
                <div class="reaction-header" onclick="toggleDetails(event, this)" style="cursor: pointer;">
                    <span class="toggle-icon">−</span>
                    <text>Description:</text>
                    <input type="text" class="reaction-name-input" placeholder="Reaction Description" value="${reaction.fields.short_name}" data-reaction-id="${reaction.pk}">
                    <text>Abbreviation:</text>
                    <input type="text" class="reaction-abbreviation-input" placeholder="Enter reaction abbreviation" value="" data-reaction-id="${reaction.pk}">
                    <label for="confidencedropdown-${reactionId}">Confidence Score:</label>
                    <select class="confidencedropdown" id="confidencedropdown-${reactionId}" data-reaction-id="${reactionId}">
                        <option value=" ">-</option>
                        <option value="1">1</option>
                        <option value="2">2</option>
                        <option value="3">3</option>
                        <option value="4">4</option>
                    </select>
                    
                    <div class= "cs-info" onclick="toggleInfo()" style="cursor: pointer;">
                    <button class="cs-info-button">i</button>
                    </div>
                </div>

                <div class="reaction-details">
                    <p>Substrates:</p>
                    <table>
                        <thead>
                            <tr>
                                <th>Stoichiometry</th>
                                <th>Comp</th>
                                <th>Name</th>
                                <th>Abbreviation</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${substrates_names.map((name, index) => `
                                <tr class="detail-item" data-reaction-id="${reaction.pk}" data-tooltip-content="${formatTooltipContent(substrates[index], substrates_types[index], subs_comps[index])}">
                                    <td>
                                        <text>${subs_stoich[index]}</text>
                                    </td>
                                    <td>
                                        <text>${subs_comps[index]}</text>
                                    </td>
                                    <td>
                                        <input type="text" name=subsNameInput placeholder="Name" value="${name}" ${subsInVMHForReaction[index]  ? 'readonly' : ''}>
                                        ${subsNeedNewNamesForReaction[index] ? '<span class="info-icon">&#63;</span><div class="tooltip-content">The name ' + name + ' is already assigned in VMH, assign another for this metabolite. Please also check that the metabolite you are adding is not already in VMH.</div>' : ''}
                                    </td>
                                    <td>
                                        <input type="text" name=subsAbbrInput placeholder="Abbreviation" value="${subsAbbrForReaction[index]}" ${subsInVMHForReaction[index]  ? 'readonly' : ''}>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
                <div class="reaction-details">
                    <p>Products:</p>
                    <table>
                        <thead>
                            <tr>
                                <th>Stoichiometry</th>
                                <th>Comp</th>
                                <th>Name</th>
                                <th>Abbreviation</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${products_names.map((name, index) => `
                                <tr class="detail-item" data-reaction-id="${reaction.pk}" data-tooltip-content="${formatTooltipContent(products[index], products_types[index], prods_comps[index])}">
                                    <td>
                                        <text>${prods_stoich[index]}</text>
                                    </td>
                                    <td>
                                        <text>${prods_comps[index]}</text>
                                    </td>
                                    <td>
                                        <input type="text" name=prodsNameInput placeholder="Name" value="${name}" ${prodsInVMHForReaction[index] ? 'readonly' : ''}>
                                        ${prodsNeedNewNamesForReaction[index] ? '<span class="info-icon">&#63;</span><div class="tooltip-content">The name ' + name + ' is already assigned in VMH, assign another for this metabolite. Please also check that the metabolite you are adding is not already in VMH.</div>' : ''}
                                    </td>
                                    <td>
                                        <input type="text" name=prodsAbbrInput placeholder="Abbreviation" value="${prodsAbbrForReaction[index]}" ${prodsInVMHForReaction[index] ? 'readonly' : ''}>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
                `;
                // Append references, external links, comments as before
                let references = reaction.fields.references || [];
                let ext_links = reaction.fields.ext_links || [];
                let comments = reaction.fields.comments || [];
                let gene_info = reaction.fields.gene_info || [];

                gene_info = gene_info.map(item => {
                    if (item.info) {
                        // Split the 'info' field and update it
                        item.info = item.info.split(';')[0];
                    }
                    return item;
                });

                listItem.innerHTML += createSectionHTML('References', 'reference', references, reaction.pk, false,true);
                listItem.innerHTML += createSectionHTML('External Links', 'ext-link', ext_links, reaction.pk, true,false);
                listItem.innerHTML += createSectionHTML('Comments', 'comment', comments, reaction.pk);      
                listItem.innerHTML += createSectionHTML('Gene Info', 'gene-info', gene_info, reaction.pk,false,false,true);
                modalList.appendChild(listItem);
            }
        });
        setupTooltips();
        displayValidationMessage(false); // Clear any previous validation message
        if (!eventListenersAttached) {
            attachDynamicEventListeners();
            eventListenersAttached = true;
        }
    }
});


function toggleInfo() {
    // Create modal div
    var modal = document.createElement('div');
    modal.id = 'myModal';
    modal.className = 'modal';

    // Create modal content div
    var modalContent = document.createElement('div');
    modalContent.className = 'modal-content';

    // Create close button
    // var closeButton = document.createElement('span');
    // closeButton.innerHTML = '&times;';
    // closeButton.className = 'close';
    // closeButton.onclick = function() {
    //     modal.style.display = 'none';
    // };

    // // Append close button to modal content
    // modalContent.appendChild(closeButton);

    // Create modal header
    var modalHeader = document.createElement('h2');
    modalHeader.innerText = 'Confidence score information';
    modalContent.appendChild(modalHeader);

    // Create modal body text
    var modalBody = document.createElement('p');
    modalBody.innerText = 'At the core of the human metabolism resource lies the human metabolic reconstruction, which is amenable to metabolic modeling. Consequently, metabolic functions without genetic evidence but with physiological evidence have also been included in this resource. Many of the reactions are associated with a confidence score indicating the evidence supporting their inclusion into the human metabolic reconstruction.';
    modalContent.appendChild(modalBody);

    // Create table
    var table = document.createElement('table');
    var tableHtml = `
    <table>
    <tr>
        <th>Evidence type</th>
        <th style="text-align: center;">Confidence score</th>
        <th>Examples</th>
    </tr>
    <tr>
        <td>Biochemical data</td>
        <td style="text-align: center;">4</td>
        <td>Direct evidence for gene product function and biochemical reaction: Protein purification, biochemical assays, experimentally solved protein structures, and comparative gene-expression studies.</td>
    </tr>
    <tr>
        <td>Genetic data</td>
        <td style="text-align: center;">3</td>
        <td>Direct and indirect evidence for gene function: Knock-out characterization, knock-in characterization, and over-expression.</td>
    </tr>
    <tr>
        <td>Physiological data</td>
        <td style="text-align: center;">2</td>
        <td>Indirect evidence for biochemical reactions based on physiological data: secretion products or defined medium components serve as evidence for transport and metabolic reactions.</td>
    </tr>
    <tr>
        <td>Sequence data</td>
        <td style="text-align: center;">2</td>
        <td>Evidence for gene function: Genome annotation.</td>
    </tr>
    <tr>
        <td>Modeling data</td>
        <td style="text-align: center;">1</td>
        <td>No evidence is available but reaction is required for modeling. The included function is a hypothesis and needs experimental verification. The reaction mechanism may be different from the included reaction(s).</td>
    </tr>
    <tr>
        <td>Not evaluated</td>
        <td style="text-align: center;">0</td>
        <td></td>
    </tr>
</table>

    `;
    table.innerHTML = tableHtml;
    modalContent.appendChild(table);

    // Create citation
    var citation = document.createElement('p');
    citation.innerText = 'Taken from Thiele, I., Palsson, B. O., "A protocol for generating a high-quality genome-scale metabolic reconstruction.", Nat Protocols, 5(1): 93 - 121 (2010).';
    modalContent.appendChild(citation);

    // Create OK button
    var okButton = document.createElement('button');
    okButton.innerText = 'Ok';
    okButton.onclick = function() {
        modal.style.display = 'none';
    };
    modalContent.appendChild(okButton);

    // Append modal content to modal
    modal.appendChild(modalContent);

    // Append modal to body
    document.body.appendChild(modal);

    // Display the modal
    modal.style.display = 'block';
}

function addToVMH() {
    document.getElementById('loadingIndicator').style.display = "flex";
    let loadingTexts = [
        "Checking if all requirements are met",
        "Initialising Cobra ToolBox",
        "Adding reactions to VMH database"
    ];
    let currentTextIndex = 0;
    document.getElementById('loadingText').textContent = loadingTexts[currentTextIndex];

    // Disable the confirm button to prevent multiple submissions
    setButtonState(true);

    if (!validateInputs()) {
        displayValidationMessage(true, 'Fill in all non-VMH metabolite names');
        document.getElementById('loadingIndicator').style.display = "none"; // Hide loading indicator
        setButtonState(false);
        return;
    } else {
        displayValidationMessage(false);
    }        
    // Function to update the loading text periodically
    const updateLoadingText = () => {
        currentTextIndex++;
        if (currentTextIndex < loadingTexts.length) {
            document.getElementById('loadingText').textContent = loadingTexts[currentTextIndex];
        }
    };
    let textChangeInterval = setInterval(updateLoadingText, 7500);
    const updatedReactions = checkedReactions.map(reactionId => {
        const reactionIdNum = Number(reactionId);
        const reaction = reactions.find(r => r.pk === reactionIdNum);
    
        // Gather Reaction Name and Abbreviation
        const nameInputField = document.querySelector(`.reaction-name-input[data-reaction-id="${reactionId}"]`);
        const abbrInputField = document.querySelector(`.reaction-abbreviation-input[data-reaction-id="${reactionId}"]`);
        if (nameInputField && nameInputField.value.trim() !== '') {
            reaction.fields.short_name = nameInputField.value.trim();
        }
        let reactionAbbr = "";
        if (abbrInputField && abbrInputField.value.trim() !== '') {
            reactionAbbr = abbrInputField.value.trim();
        }
        // Gather Substrates and Products Details
        const detailRows = document.querySelectorAll(`tr.detail-item[data-reaction-id="${reactionId}"]`);
        const subsDetails = [], prodsDetails = [];

        detailRows.forEach(row => {
            const nameInput = row.querySelector('input[type="text"][name*="NameInput"]');
            const abbrInput = row.querySelector('input[type="text"][name*="AbbrInput"]');
            const detail = {
                name: nameInput.value.trim(),
                abbreviation: abbrInput.value.trim()
            };
            
            // Determine if row is a substrate or product by its position
            if (nameInput.name === 'subsNameInput') {
                subsDetails.push(detail);
            } else {
                prodsDetails.push(detail);
            }
        });
            
            
        const inputField = document.querySelector(`.reaction-name-input[data-reaction-id="${reactionId}"]`);
        if (inputField && inputField.value.trim() !== '') {
            reaction.fields.short_name = inputField.value.trim();
        }

        const referencesData = Array.from(document.querySelectorAll(`.reference-input[data-reaction-id="${reactionId}"]`))
        .map(input => {
            const select = input.previousElementSibling;
            return {
                ref_type: select ? select.value : '', // Ensure select is found, otherwise default to empty string
                info: input.value.trim() // Use the input value directly
            };
        });
        // Gathering External Links Data
        const extLinksData = Array.from(document.querySelectorAll(`.ext-link-item[data-reaction-id="${reactionId}"]`))
        .map(container => {
            const select = container.querySelector('.ext-link-type-select');
            const input = container.querySelector('.ext-link-input');
            return {
                ext_link_type: select ? select.value : '', // Ensure select is found, otherwise default to empty string
                info: input.value.trim() // Use the input value directly
            };
        });
        // Gathering Comments Data
        const commentsData = Array.from(document.querySelectorAll(`.comment-input[data-reaction-id="${reactionId}"]`))
                        .map(input => ({ info: input.value.trim() }));
        
        //Gathering Confidence Score Data
        const confidenceDropdown = document.querySelector(`.confidencedropdown[data-reaction-id="${reactionId}"]`);
        const selectedValue = confidenceDropdown.value;
    
        // Returning the updated reaction object
        return {
            pk: reaction.pk,
            short_name: reaction.fields.short_name, 
            abbreviation: reactionAbbr, // Add reaction abbreviation
            substrates_info: JSON.stringify(subsDetails), // Include substrates names and abbreviations
            products_info: JSON.stringify(prodsDetails), // Include products names and abbreviations
            references: JSON.stringify(referencesData),
            ext_links: JSON.stringify(extLinksData),
            comments: JSON.stringify(commentsData),
            confidence_score: JSON.stringify(selectedValue)
        };
        });

    fetch(addToVMHUrl, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrfToken,
        },
        body: JSON.stringify({
            userID: userID,
            reactions: updatedReactions
        })
    })
    .then(response => response.json())
    .then(data => {
        clearInterval(textChangeInterval); // Stop updating the loading text
        document.getElementById('loadingIndicator').style.display = "none"; // Hide loading indicator
        setButtonState(false);
        if (data.status === 'success') {
            const reactions = Object.entries(data.rxn_added_info); // Array of [abbr, [rxn_ID, reaction_formula]]
            const metabolites = Object.entries(data.met_added_info); // Array of [abbr, [met_ID, formula, inchiKey]]
    
            const reactionPart = reactions.length > 0 ? `${reactions.length} ${reactions.length === 1 ? 'reaction' : 'reactions'}` : '';
            const metabolitePart = metabolites.length > 0 ? `${metabolites.length} ${metabolites.length === 1 ? 'metabolite' : 'metabolites'}` : '';
            const conjunction = reactions.length > 0 && metabolites.length > 0 ? 'and' : '';
            const toDatabasePart = 'to VMH database';
    
            document.getElementById('responseModalTitle').textContent = `Successfully added ${reactionPart} ${conjunction} ${metabolitePart} ${toDatabasePart}`.trim();
    
            let responseBody = '';
            if (metabolites.length > 0) {
                responseBody += '<h3>Metabolites Added:</h3><ul>';
                responseBody += '<p>Metabolite Abbreviation (ID in VMH Database):</p><ul>';
                metabolites.forEach(([abbr, [id, formula, inchiKey]]) => {
                    responseBody += `<li><span class="bold">${abbr}</span> (${id})<div class="info-block" id="metInfo_${id}" style="display:none;"><p>Formula: ${formula}</p><p>InChIKey: ${inchiKey}</p></div></li>`;
                });
                responseBody += '</ul>';
            }
    
            if (reactions.length > 0) {
                responseBody += '<h3>Reactions Added:</h3><ul>';
                responseBody += '<p>Reaction Abbreviation (ID in VMH Database):</p><ul>';
                reactions.forEach(([abbr, [id, formula]]) => {
                    responseBody += `<li><span class="bold">${abbr}</span> (${id})<div class="info-block" id="rxnInfo_${id}" style="display:none;"><p>Formula: ${formula}</p></div></li>`;
                });
                responseBody += '</ul>';
            }
    
            document.getElementById('responseModalBody').innerHTML = responseBody;
    
            // Add click event listener to toggle info display
            document.querySelectorAll('#responseModalBody ul li').forEach(item => {
                item.addEventListener('click', function() {
                    const infoBlock = this.querySelector('.info-block');
                    if (infoBlock) {
                        infoBlock.style.display = infoBlock.style.display === 'none' ? 'block' : 'none';
                    }
                });
            });
    
            // Show the modal
            document.getElementById('reactionModal').style.display = "none";
            document.getElementById('responseModal').style.display = "block";
        } else {
            // Handle error case
            document.getElementById('responseModalTitle').textContent = "Error";
            document.getElementById('responseModalBody').textContent = data.message;
            document.getElementById('responseModal').style.display = "block";
            window.scrollTo(0, 0); // Scroll to the top of the page
        }
    })                      
    .catch((error) => {
        clearInterval(textChangeInterval); // Stop updating the loading text
        console.error('Error:', error);
        document.getElementById('loadingIndicator').style.display = "none"; // Hide loading indicator
        setButtonState(false);
        window.scrollTo(0, 0);
    });
}
document.getElementById('confirmAddToVMH').addEventListener('click', addToVMH);

// Close the alert modal when the user clicks on <span> (x)
document.querySelector('.close-alert-btn').addEventListener('click', function() {
    document.getElementById('alertModal').style.display = "none";
});

// Also close the modal if the user clicks anywhere outside of the modal
window.onclick = function(event) {
    if (event.target == document.getElementById('alertModal')) {
        document.getElementById('alertModal').style.display = "none";
    }
};

// Close the modal when the user clicks on <span> (x)
document.getElementsByClassName("close")[0].onclick = function() {
    var modal = document.getElementById('reactionModal');
    modal.style.display = "none";
    window.location.reload();
}
document.getElementById('responseModal').addEventListener('click', function(event) {
    if (event.target.id === "responseModalClose" || event.target.classList.contains('close')) {
        this.style.display = "none";
        // window.location.reload();
    }
});

document.getElementById('backToHome').addEventListener('click', function() {
    window.location.href = '/';
});


