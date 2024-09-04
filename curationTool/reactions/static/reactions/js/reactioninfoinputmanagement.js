async function fetchPubMedInfo(pmid) {
    try {
        // Correctly construct the URL using the dynamic path from Django
        const url = getPubMedInfo.replace('%s', pmid);
        const response = await fetch(url, {
            method: 'GET',
            headers: {
                'X-CSRFToken': csrfToken,
                'Content-Type': 'application/json',
            },
        });
        const data = await response.json();
        return data;
    } catch (error) {
        return data = { status: 'error', message: 'Failed to fetch PubMed data.' };
    }
}
async function fetchDOIInfo(doi) {
    try {
        const url = getDOIInfo.replace('%s', doi);
        const response = await fetch(url, {
            method: 'GET',
            headers: {
                'X-CSRFToken': csrfToken,
                'Content-Type': 'application/json',
            },
        });
        const data = await response.json();
        return data;
    } catch (error) {
        return { status: 'error', message: 'Failed to fetch DOI data.' };
    }
}
function getHyperlinkUrl(ext_link_type, info) {
    // Placeholder switch case to determine the URL based on ext_link_type
    // Here you can add actual logic to return different URLs based on ext_link_type
    switch (ext_link_type) {
        case 'CHO Models':
            return `https://chomine.boku.ac.at/chomodel/reaction/${info}`;
        case 'MetanetX':
            return `https://www.metanetx.org/equa_info/${info}`;
        case 'Rhea':
            return `https://www.rhea-db.org/rhea/${info}`;
        case 'KEGG reaction':
            return `https://www.kegg.jp/entry/${info}`;
        case 'Wikipedia':
            return `https://en.wikipedia.org/wiki/${info}`;
        case 'SEED':
            return `https://modelseed.org/biochem/reactions/${info}`;
        case 'COG':
            return `ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/static/byCOG/${info}.html`
        case 'EC Number':
            return `https://google.com/search?q=not+implemented+(double+check+EC+number+${info})`;
        case 'KEGG orthology':
            return `https://www.genome.jp/dbget-bin/www_bget?ko:${info}`;
    }
}

// Function to create editable cell
function createEditableCell(text, gene, fieldType, userId, reactionId) {
    const cell = document.createElement('td');
    cell.textContent = text;
    cell.contentEditable = true;
    cell.className = 'editable-cell';

    // Event listener for saving changes on Enter key press
    cell.addEventListener('keypress', async (event) => {
        if (event.key === 'Enter') {
            event.preventDefault(); // Prevent default behavior (new line)
            cell.blur(); // Remove focus to trigger change event

            // Save the updated data
            const updatedValue = cell.textContent.trim();

            // Check if updatedValue is not empty before saving
            if (updatedValue) {
                var userId = sessionStorage.getItem('userID'); 

                await saveUpdatedGeneInfo(gene, fieldType, updatedValue, userId, reactionId);
            } else {
                console.warn('Empty value. No changes saved.');
            }
        }
    });

    return cell;
}

async function saveUpdatedGeneInfo(gene, fieldType, updatedValue, userId, reactionId) {
    try {
        const response = await fetch('/update_gene_info/', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                userID: userId,
                reactionID: reactionId,
                gene: gene,
                fieldType: fieldType,
                updatedValue: updatedValue,
            }),
        });

        const result = await response.json();
        if (response.ok) {
            console.log('Gene information updated successfully:', result);
            // Optionally show a success message to the user
        } else {
            console.error('Error updating gene information:', result);
            // Optionally show an error message to the user
        }
    } catch (error) {
        console.error('Unexpected error:', error);
        // Optionally show an error message to the user
    }
}


// Function to fetch parsed gene info from the Django backend
async function fetchParsedGeneInfo(info) {
    try {
        const response = await fetch(`/parse_gene_info/?info=${encodeURIComponent(info)}`);
        if (response.ok) {
            const data = await response.json();
            return data;
        } else {
            console.error('Failed to fetch parsed gene info:', response.statusText);
            return {};
        }
    } catch (error) {
        console.error('Error occurred while fetching parsed gene info:', error);
        return {};
    }
}


function deleteRow(reactionID, item, tabId, rowElement) {
    const confirmation = confirm("Are you sure you want to delete this item?");
    if (!confirmation) {
        return; // Exit if the user cancels the deletion
    }

    const url = deleteInfo
    const body = JSON.stringify({ reaction_id: reactionID, item_to_delete: item, tab_id: tabId });

    fetch(url, {
        method: 'POST',
        headers: {
            'X-CSRFToken': csrfToken,
            'Content-Type': 'application/json',
        },
        body: body,
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            console.error('Error deleting row:', data.message);
        } else {
            // remove row
            rowElement.remove();
        }
    })
    .catch(error => console.error('Error:', error));
}



