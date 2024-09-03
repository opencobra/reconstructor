function loadChemInfoDiv(reactionData){
    let contentDiv = document.querySelector('div.content-div[name="cheminfo-div"]');
    contentDiv.innerHTML = `
    <div class="div-header">Chemical Information</div>
    `; 
    maindiv = document.createElement('div');
        // Create and append the balanced text
    let balancedText = document.createElement('p');
    let isBalanced = reactionData.balanced_count[0] && reactionData.balanced_charge[0];
    balancedText.textContent = isBalanced ? 'Balanced' : 'Not Balanced';
    balancedText.style.color = isBalanced ? '#4CAF50' : '#F44336'; // Green if balanced, red if not
    balancedText.style.fontSize = '18px'; // Set font size
    balancedText.style.fontWeight = 'bold'; // Make the text bold
    balancedText.style.padding = '10px'; // Add some padding
    balancedText.style.margin = '0'; // Adjust margin if necessary
    balancedText.style.borderRadius = '4px'; // Optional: add border radius for a subtle rounded corner effect
    balancedText.style.display = 'inline-block'; // Adjust display type for better control
    balancedText.style.backgroundColor = isBalanced ? '#E8F5E9' : '#FFEBEE'; // Light background color for contrast
    balancedText.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)'; // Optional: add a subtle shadow for depth
    balancedText.style.marginBottom = '20px'; // Add some margin at the bottom
    balancedText.style.height = '40%'; // Set the height of the element
    balancedText.style.width = '40%'; // Center the text vertically
    maindiv.appendChild(balancedText);
    // Create and append atom types and their colors
    let atomTypes = new Set([...Object.keys(reactionData.subs_atoms[0]), ...Object.keys(reactionData.prods_atoms[0])]);
    let colors = getDistinctColors(atomTypes.size);
    let atomColorMap = {};
    let i = 0;
    atomTypes.forEach(atom => {
        atomColorMap[atom] = colors[i++];
    });

    // Create and append the molecular formula
    let formulaParagraph = document.createElement('p');
    formulaParagraph.textContent = 'Molecular Formula: ';
    let formulaSpan = document.createElement('span');
    formulaSpan.innerHTML = reactionData.molc_formula[0].replace(/([A-Z][a-z]?)/g, match => {
        return atomColorMap[match] ? `<span style="color: ${atomColorMap[match]}">${match}</span>` : match;
    });
    formulaParagraph.appendChild(formulaSpan);
    maindiv.appendChild(formulaParagraph);

    // Create and append mass balance information
    let balancedInfoDiv = document.createElement('div');
    balancedInfoDiv.className = 'balanced-info';
    balancedInfoDiv.textContent = `Mass Balanced? ${reactionData.balanced_count[0] ? 'Yes' : 'No'}`;
    maindiv.appendChild(balancedInfoDiv);

    // Create and append charge balance information
    let chargeBalancedInfoDiv = document.createElement('div');
    chargeBalancedInfoDiv.className = 'balanced-info';
    chargeBalancedInfoDiv.textContent = `Charge Balanced? ${reactionData.balanced_charge[0] ? 'Yes' : 'No'}`;
    maindiv.appendChild(chargeBalancedInfoDiv);
    


// Function to show the loader with a modal background and a specific message
function showLoaderdiv(message) {
    modalBackground.style.display = 'block';
    loader.textContent = message;
    loader.style.display = 'block';
}

// Function to hide the loader and the modal background
function hideLoaderdiv() {
    loader.style.display = 'none';
    modalBackground.style.display = 'none';
}
    // Initialize the container div
    function createSpan(text) {
        let span = document.createElement('span');
        span.textContent = text;
        return span;
    }
    let textElement = document.createElement('div');
    textElement.textContent = 'VMH Formula: ';
    textElement.style.marginRight = '10px';  // Add some space between the text and the div
    textElement.style.display = 'flex';  // Adjust display type for better control
    textElement.style.justifyContent   = 'center';  // Center the text horizontally
    textElement.style.alignItems = "baseline";
    maindiv.appendChild(textElement);

    let formulaDiv = document.createElement('div');
    formulaDiv.style.padding = "10px";
    formulaDiv.style.marginTop = "20px";
    
    let modalBackground = document.createElement('div');
    modalBackground.style.display = 'none'; // Initially hidden
    modalBackground.style.position = 'fixed';
    modalBackground.style.top = '0';
    modalBackground.style.left = '0';
    modalBackground.style.width = '100%';
    modalBackground.style.height = '100%';
    modalBackground.style.backgroundColor = 'rgba(0, 0, 0, 0.75)'; // Semi-transparent black
    modalBackground.style.zIndex = '999'; // Ensure it is below the loader but above other elements
    formulaDiv.appendChild(modalBackground);

    // Create the loader element
    // Create the loader element
    let loader = document.createElement('div');
    loader.style.display = 'none'; // Initially hidden
    loader.style.position = 'absolute'; // Position relative to .content-div
    loader.style.top = '50%';
    loader.style.left = '50%';
    loader.style.transform = 'translate(-50%, -50%)'; // Center the loader within .content-div
    loader.style.padding = '20px';
    loader.style.backgroundColor = 'white';  // Set loader background to white
    loader.style.color = 'black';  // Set text color to black for contrast
    loader.style.borderRadius = '5px';
    loader.style.zIndex = '1000'; // Ensure it is above other elements inside .content-div
    loader.style.textAlign = 'center';  // Center the text inside the loader
    loader.style.boxShadow = '0px 10px 15px rgba(0, 0, 0, 0.3)'; // Optional: Add shadow for better visibility
    loader.textContent = 'Loading...'; // Example text, change as needed

    formulaDiv.appendChild(loader);

    if(reactionData.rxn_formula== null){


    // Helper function to create a span with the given text
    function createSpan(text) {
        let span = document.createElement('span');
        span.textContent = text;
        return span;
    }

    // Helper function to create a span for found metabolites
    function createSpanElement(metaboliteName, compartment) {
        let span = document.createElement('span');
        span.textContent = `${metaboliteName} [${compartment}]`;
        return span;
    }

    // Helper function to create an empty input field with a dropdown and button for metabolites that are not found
    function createEmptyField(metaboliteName, metabolite, mtype, compartment) {
        let input = document.createElement('input');
        input.type = 'text';
        input.placeholder = metaboliteName; // Placeholder only shows metabolite name
        input.disabled = true; // Initially disable the input field

        let dropdown = document.createElement('select');
        let optionGenAbbr = document.createElement('option');
        optionGenAbbr.text = 'Gen-Abbr';
        optionGenAbbr.value = 'gen-abbr';
        dropdown.add(optionGenAbbr);

        let optionWriteOwn = document.createElement('option');
        optionWriteOwn.text = 'Write your own abbr';
        optionWriteOwn.value = 'write-own';
        dropdown.add(optionWriteOwn);

        let button = document.createElement('button');
        button.textContent = 'Apply';

        // Add event listener to the dropdown
        dropdown.addEventListener('change', function() {
            if (dropdown.value === 'write-own') {
                input.disabled = false; // Enable input field if "Write your own abbr" is selected
                button.textContent = 'Save';
            } else {
                input.disabled = true;
                button.textContent = 'Gen-Abbr';
            }
        });

        // Add event listener to the button
        button.addEventListener('click', function() {
            if (dropdown.value === 'gen-abbr') {
                // Create the data object to be sent
                showLoaderdiv('Generating abbreviation...');
                const dataToSend = {
                    metabolite: metabolite,
                    mtype: mtype,
                    metabolite_name: metaboliteName,
                    compartment: compartment
                };

                // Log the data to the console

                // Send a POST request to the Django view with JSON data
                fetch('/create-formula-abbr/', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json'
                    },
                    body: JSON.stringify(dataToSend)
                })
                .then(response => response.json())
                .then(data => {
                    if (data.abbr) {
                        hideLoaderdiv();
                        // Create a span to replace the input, dropdown, and button
                        let resultSpan = document.createElement('span');
                        resultSpan.textContent = `${data.abbr} [${compartment}]`;

                        // Replace the input, dropdown, and button with the span
                        input.replaceWith(resultSpan);
                        dropdown.remove();  // Remove the dropdown
                        button.remove();  // Remove the button
                    } else {
                        console.error('Failed to generate abbreviation');
                    }
                })
                .catch(error => console.error('Error:', error));
            } else if (dropdown.value === 'write-own') {
                // Save the user input as the abbreviation
    
                let userAbbr = input.value;
                if (userAbbr) {
                    let resultSpan = document.createElement('span');
                    resultSpan.textContent = `${userAbbr} [${compartment}]`;

                    // Replace the input, dropdown, and button with the span
                    input.replaceWith(resultSpan);
                    dropdown.remove();  // Remove the dropdown
                    button.remove();  // Remove the button
                }
            }
        });

        let container = document.createElement('div');
        container.style.display = 'inline-block';
        container.style.marginRight = '10px';
        container.appendChild(input);
        container.appendChild(dropdown);
        container.appendChild(button);
        return container;
    }

    // Check if formulas is empty
    if (reactionData.formulas && reactionData.formulas.trim() !== '') {
        // If formulas already exist, display them directly
        let existingFormulaSpan = createSpan(reactionData.formulas);
        formulaDiv.appendChild(existingFormulaSpan);
    } else {
        // If formulas is empty, generate the formula as before
        reactionData.subs_sch.forEach((stoich, index) => {
            if(stoich >1){
            formulaDiv.appendChild(createSpan(parseFloat(stoich).toFixed(1) + " "));}
            if (reactionData.subs_found[index]) {
                formulaDiv.appendChild(createSpanElement(reactionData.substrates[index], reactionData.subs_comps[index]));
            } else {
                formulaDiv.appendChild(createEmptyField(
                    reactionData.substrates_names[index], 
                    reactionData.substrates[index], 
                    reactionData.subs_types[index],
                    reactionData.subs_comps[index]  // Pass the compartment for substrates
                ));
            }
            // Add a plus sign after each substrate except the last one
            if (index < reactionData.subs_sch.length - 1) {
                formulaDiv.appendChild(createSpan(" + "));
            }
        });

                // Add direction symbol (e.g., "->")
        const directionIcon = reactionData.direction === 'forward' ? '->' : reactionData.direction === 'bidirectional' ? '<=> ' : ''; 
        formulaDiv.appendChild(createSpan(directionIcon));


        reactionData.prod_sch.forEach((stoich, index) => {
            if(stoich >1){
                formulaDiv.appendChild(createSpan(parseFloat(stoich).toFixed(1) + " "));}
            if (reactionData.prod_found[index]) {
                formulaDiv.appendChild(createSpanElement(reactionData.products[index], reactionData.prods_comps[index]));
            } else {
                formulaDiv.appendChild(createEmptyField(
                    reactionData.products_names[index], 
                    reactionData.products[index], 
                    reactionData.prods_types[index],
                    reactionData.prods_comps[index]  // Pass the compartment for products
                ));
            }
            // Add a plus sign after each product except the last one
            if (index < reactionData.prod_sch.length - 1) {
                formulaDiv.appendChild(createSpan(" + "));
            }
        });
    }

    // Add a "Save All" button to trigger all "Apply" buttons
    let saveAllButton = document.createElement('button');
    saveAllButton.textContent = 'Save All';
    saveAllButton.style.marginTop = '20px';
    saveAllButton.style.marginLeft = '10px';
    saveAllButton.id = 'save-all-button';


    saveAllButton.addEventListener('click', function() {
        // Find all "Apply" buttons and trigger a click event on each
        let applyButtons = formulaDiv.querySelectorAll('button');
        applyButtons.forEach(button => {
            if (button.textContent === 'Apply' || button.textContent === 'Save') {
                button.click(); // Trigger click event
            }
        });

        // After all inputs are converted to spans, collect the formula string
        let spans = formulaDiv.querySelectorAll('span');
        let formulaString = '';
        spans.forEach(span => {
            formulaString += span.textContent;
        });

        // Assuming you have the reaction ID available
        const reactionId = reactionData.reaction_id;  // Use the actual reaction ID
        showLoaderdiv('Saving formula...');
        // Send the formula string and reaction ID to Django backend
        fetch('/save-formula/', {  // Replace with your Django view URL
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-CSRFToken': csrfToken  // Ensure CSRF token is correctly passed
            },
            body: JSON.stringify({ 
                formula: formulaString,
                reaction_id: reactionId  // Pass the reaction ID here
            })
        })
        .then(response => {
            if (response.ok) {
                
                // Adding a 100-second timeout before hiding the loader
                setTimeout(function() {
                    hideLoaderdiv();
                    document.getElementById('save-all-button').remove(); // Replace 'yourButtonId' with the actual button ID
                }, 1000); // 100 seconds = 100000 milliseconds
            
            }else {
                console.error('Failed to save formula');
            }
        })
        .catch(error => console.error('Error:', error));
    });
    
    formulaDiv.appendChild(saveAllButton);
    textElement.appendChild(formulaDiv);
    maindiv.appendChild(textElement);}
    else{
        let existingFormulaSpan = createSpan(reactionData.rxn_formula);
        formulaDiv.appendChild(existingFormulaSpan);
        textElement.appendChild(formulaDiv);

    }



    // Create and append the atom comparison table
    let comparisonListDiv = document.createElement('div');
    comparisonListDiv.className = 'atoms-comparison';
    populateAtomComparisonTable(
        reactionData.subs_atoms[0],
        reactionData.prods_atoms[0],
        comparisonListDiv,
        { substrates: reactionData.subs_charge[0], products: reactionData.prods_charge[0] },
        reactionData.symb_to_name[0]
    );
    maindiv.appendChild(comparisonListDiv);
    maindiv.style.textAlign = 'center';
    contentDiv.appendChild(maindiv);
    contentDiv.style.display = 'block';
}



function populateAtomComparisonTable(substratesData, productsData, listDiv, charge, symb_to_name) {
    listDiv.style.display = 'flex';
    listDiv.style.justifyContent = 'space-between';
    listDiv.style.alignItems = 'flex-start';
    listDiv.style.padding = '20px';
    listDiv.style.border = '1px solid #ddd';
    listDiv.style.borderRadius = '8px';
    listDiv.style.marginTop = '5%';
    listDiv.style.backgroundColor = '#f9f9f9';
    listDiv.style.fontSize = '16px';  // Replace '14px' with your desired font size

    const allAtoms = { ...substratesData, ...productsData };

    const table = document.createElement('table');
    table.className = 'atom-comparison-table';
    table.style.width = '70%';
    table.style.borderCollapse = 'collapse';

    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');
    ['Atom', 'Substrates', 'Products'].forEach(text => {
        const th = document.createElement('th');
        th.textContent = text;
        th.style.borderBottom = '2px solid #000';
        th.style.padding = '10px';
        th.style.backgroundColor = '#eaeaea';
        th.style.textAlign = 'center';
        th.style.fontWeight = 'bold';
        th.style.fontSize = '19px';  // Replace '16px' with your desired font size
        headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);

    const tbody = document.createElement('tbody');
    Object.keys(allAtoms).forEach(atom => {
        const row = document.createElement('tr');
        const substratesCount = substratesData[atom] || 0;
        const productsCount = productsData[atom] || 0;

        if (substratesCount !== productsCount) {
            row.style.backgroundColor = '#ffcccc'; // Highlight in red
        }

        [symb_to_name[atom] || atom, substratesCount, productsCount].forEach((text, index) => {
            const cell = document.createElement('td');
            cell.textContent = text;
            cell.style.padding = '8px';
            cell.style.borderBottom = '1px solid #ddd';
            cell.style.textAlign = 'center';
            if (index === 0) {
                cell.style.fontWeight = 'bold';
            }
            row.appendChild(cell);
        });
        tbody.appendChild(row);
    });
    table.appendChild(tbody);

    listDiv.appendChild(table);

    const chargeDiv = document.createElement('div');
    chargeDiv.innerHTML = `<strong>Charge Info</strong><br>Substrates: ${charge.substrates}<br>Products: ${charge.products}`;
    chargeDiv.style.marginLeft = '20px';
    chargeDiv.style.padding = '10px';
    chargeDiv.style.border = '1px solid #ccc';
    chargeDiv.style.borderRadius = '8px';
    chargeDiv.style.backgroundColor = '#fff';
    chargeDiv.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
    listDiv.appendChild(chargeDiv);
}

function getDistinctColors(count) {
    // Function to generate distinct colors for atom types
    const colors = [];
    for (let i = 0; i < count; i++) {
        const hue = (i * 137.508) % 360; // Use golden angle approximation
        colors.push(`hsl(${hue}, 100%, 50%)`);
    }
    return colors;
}



