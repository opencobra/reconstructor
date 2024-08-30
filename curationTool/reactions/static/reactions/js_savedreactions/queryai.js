var llm_html = '';

document.getElementById('AIButton').addEventListener('click', function() {
    var dropdown = document.getElementById('dropdown');
    var content = document.querySelector('.content');
    var llmResponseContainer = document.getElementById('llmResponseContainer');
    var errorContainer = document.getElementById('errorContainer');
    if (dropdown.style.display == 'none' || dropdown.style.display == '') {
        dropdown.style.display = 'block';
        errorContainer.style.display = 'block';

        if (llm_html != '') {
            llmResponseContainer.style.display = 'block';
            llmResponseContainer.innerHTML = llm_html;
        }

    } else {
        errorContainer.style.display = 'none';
        dropdown.style.display = 'none';
        llmResponseContainer.style.display = 'none';
    }
});
document.getElementById('submitGene').addEventListener('click', function() {
    var gene = document.getElementById('geneInput').value;
    var temperature = parseFloat(document.getElementById('temperatureSlider').value);  // Get the temperature value
    var submitButton = document.getElementById('submitGene');
    var loadingSpinner = document.getElementById('loadingSpinner');
    var errorContainer = document.getElementById('errorContainer');

    // Clear previous error message
    errorContainer.style.display = 'none';
    errorContainer.innerText = '';

    // Disable the submit button and show the loading spinner
    submitButton.disabled = true;
    loadingSpinner.style.visibility = 'visible';

    fetch('/get_ai_response/', {
        method: 'POST',  // Changed to 'POST'
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': csrftoken  // Ensure CSRF token is included if needed
        },
        body: JSON.stringify({
            key: gene,  // Gene data
            temperature: temperature  // Temperature value
        })
    }).then(response => {
        return response.json().then(data => {
            if (response.ok) {
                return data; // Ensure the server responds with JSON
            } else {
                throw new Error(data.error_message || 'Network response was not ok.');
            }
        });
    })
    .then((data) => {
        if (data.status == 'error') {
            errorContainer.innerText = data['error_message'];
            errorContainer.style.display = 'block';
            return;
        }
        // Update the LLM response container
        var llmResponseContainer = document.getElementById('llmResponseContainer');
        var json_data = data['llm_response_html'];
    
        // Parse the JSON string
        var reactions = JSON.parse(json_data);
    
        // Function to create a link if the metabolite doesn't have an underscore
        function createLink(metabolite) {
            if (!metabolite.includes('_')) {
                return `<a href="https://www.vmh.life/#metabolite/${metabolite}">${metabolite}</a>`;
            } else {
                return metabolite.replace('_', '');
            }
        }
    
        // Function to format the reaction string
        function formatReaction(reaction) {
            
            return reaction.map((metabolite, index) => {
                if (metabolite == '->_') {
                    return '->';
                } else {
                    const linkedMetabolite = createLink(metabolite);
                    if (index == 0 || reaction[index - 1] == '->_') {
                        return `${linkedMetabolite}`;
                    } else {
                        return `+ ${linkedMetabolite}`;
                    }
                }
            }).join(' ');
        }
    
        var formattedReactions = reactions.map((reaction, reactionIndex) => {
            var formattedReaction = formatReaction(reaction);
            return `<p>${reactionIndex + 1}. ${formattedReaction}</p>`;
        }).join('');
    
        llmResponseContainer.innerHTML = formattedReactions;
        llm_html = formattedReactions;
        if (llmResponseContainer.innerHTML == '') {
            errorContainer.innerText = ' Apologies, no reactions predicted by the AI for that gene :('
            errorContainer.style.display = 'block';
            llmResponseContainer.style.display = 'none';
        }
        else {
            errorContainer.style.display = 'none';
            llmResponseContainer.style.display = 'block';
        }
        var content = document.querySelector('.content');
        content.style.paddingTop = '30px';
    }).catch((error) => {
        console.log(error);
        errorContainer.innerText = error.message;
        errorContainer.style.display = 'block';
        console.error('Error:', error);
    }).finally(() => {
        submitButton.disabled = false;
        loadingSpinner.style.visibility = 'hidden';
    });
});

// Update temperature value display
document.getElementById('temperatureSlider').addEventListener('input', function() {
    document.getElementById('temperatureValue').textContent = this.value;
});

function getCookie(name) {
    let cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        const cookies = document.cookie.split(';');
        for (let i = 0; i < cookies.length; i++) {
            const cookie = cookies[i].trim();
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}

const csrftoken = getCookie('csrftoken');