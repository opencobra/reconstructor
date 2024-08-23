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
    listDiv.style.marginTop = '20px';
    listDiv.style.backgroundColor = '#f9f9f9';

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