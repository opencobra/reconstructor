document.addEventListener('DOMContentLoaded', () => {
    const searchInput = document.getElementById('searchInput');
    const reactionListBody = document.querySelector('#reactionList tbody');

    searchInput.addEventListener('input', () => {
        const query = searchInput.value.toLowerCase();
        const rows = reactionListBody.getElementsByTagName('tr'); // Only get rows from tbody

        Array.from(rows).forEach(row => {
            // Check if the row is currently visible
            if (row.style.display !== 'none') {
                const name = row.cells[1].textContent.toLowerCase();
                const subsystem = row.cells[2].textContent.toLowerCase();
                const substrates = row.cells[3].textContent.toLowerCase();
                const products = row.cells[4].textContent.toLowerCase();

                if (name.includes(query) || subsystem.includes(query) || substrates.includes(query) || products.includes(query)) {
                    row.style.display = ''; // Show row
                } else {
                    row.style.display = 'none'; // Hide row
                }
            }
        });
    });
});
