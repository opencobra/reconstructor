// Generates a list of distinct colors in HSL format, avoiding yellow hues, based on the specified count.
function getDistinctColors(count) {
    var colors = [];
    var hueStep = 360 / count;
    var lightness = 40; // Adjust for darker colors, but not too dark

    for (var i = 0; i < 360; i += hueStep) {
        // Skip yellow and near-yellow hues
        if (i >= 50 && i <= 70) {
            continue;
        }

        colors.push('hsl(' + i + ', 100%, ' + lightness + '%)');
    }
    return colors;
}

window.ReactionUtils = window.ReactionUtils || {};

window.ReactionUtils.fetchTemplates = async function (templateList, templatesFetched) {
    if (templatesFetched.value) return;
    try {
        const userID = sessionStorage.getItem('userID');
        console.log(userID);
        const response = await fetch('/list_templates/', {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ userID: userID }),
        });

        if (response.ok) {
            const data = await response.json();
            templateList.length = 0; // Clear existing list
            templateList.push(...data.templates); // Add new templates
            templatesFetched.value = true;
        } else {
            console.error('Failed to fetch templates');
        }
    } catch (error) {
        console.error('Error:', error);
    }
};

// Function to show a toast notification
function showToast(message, color = '#4caf50') {
    // Unified through the Notify toast system (notify.js) so every toast in the
    // app shares one neat style. Legacy callers pass a colour hint; map a red/
    // error hint to an error toast, otherwise treat it as success.
    let type = 'success';
    if (typeof color === 'string') {
        const c = color.toLowerCase();
        if (/error|cc0000|red|#f|#d|#c|#b/.test(c)) {
            type = 'error';
        }
    }
    if (window.Notify) {
        return Notify.toast(message, { type });
    }
    // Minimal fallback if notify.js is unavailable on a given page.
    const toast = document.createElement('div');
    toast.textContent = message;
    toast.style.cssText =
        'position:fixed;bottom:20px;right:20px;background:' + color +
        ';color:#fff;padding:10px 20px;border-radius:8px;box-shadow:0 4px 8px rgba(0,0,0,0.2);font-size:14px;z-index:5000;';
    document.body.appendChild(toast);
    setTimeout(() => toast.remove(), 3000);
}

// Function to save the current visibility state of each div
function updateLocalStorageState() {
    const visibleState = {};
    for (const [buttonId, divName] of Object.entries(button_to_div)) {
        const div = document.getElementsByName(divName + '-div')[0];
        if (div) {
            const isVisible = div.classList.contains('is-visible') || div.style.display === 'block';
            visibleState[divName] = isVisible;
        }
    }
    localStorage.setItem('visibleDivs', JSON.stringify(visibleState));
}

// Function to restore the saved visibility state on page load
function restoreVisibleDivs() {
    const savedState = localStorage.getItem('visibleDivs');
    if (!savedState) return;
    const visibleState = JSON.parse(savedState);
    for (const [divName, isVisible] of Object.entries(visibleState)) {
        const div = document.getElementsByName(divName + '-div')[0];
        if (div) {
            div.style.display = isVisible ? 'block' : 'none';
            div.classList.toggle('is-visible', isVisible);
            // Find the corresponding button and update its active class
            const buttonId = Object.keys(button_to_div).find(key => button_to_div[key] === divName);
            if (buttonId) {
                const button = document.getElementById(buttonId);
                if (button) {
                    if (isVisible) {
                        button.classList.add('active');
                    } else {
                        button.classList.remove('active');
                    }
                }
            }
        }
    }
}
