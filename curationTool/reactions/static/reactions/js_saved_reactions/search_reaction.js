document.addEventListener('DOMContentLoaded', () => {
    const searchInput = document.getElementById('searchInput');
    const reactionListBody = document.querySelector('#reactionList tbody');

    searchInput.addEventListener('input', () => {
        const query = searchInput.value.toLowerCase();
        const rows = reactionListBody.getElementsByTagName('tr');

        Array.from(rows).forEach(row => {
            // Column order: Select(0) Name(1) CS(2) Description(3) Subsystem(4) Substrates(5)
            const name = row.cells[1].textContent.toLowerCase();
            const subsystem = row.cells[3].textContent.toLowerCase();
            const substrates = row.cells[4].textContent.toLowerCase();
            const products = row.cells[5].textContent.toLowerCase();

            if (
                name.includes(query) ||
                subsystem.includes(query) ||
                substrates.includes(query) ||
                products.includes(query)
            ) {
                row.style.display = ''; // Show row
            } else {
                row.style.display = 'none'; // Hide row
            }
        });
    });
});
