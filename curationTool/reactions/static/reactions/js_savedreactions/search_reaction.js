document.addEventListener('DOMContentLoaded', () => {
    const searchInput = document.getElementById('searchInput');
    const reactionList = document.getElementById('reactionList');

    searchInput.addEventListener('input', () => {
        const query = searchInput.value.toLowerCase();
        const rows = reactionList.getElementsByTagName('tr');

        Array.from(rows).forEach(row => {
            const name = row.cells[1].textContent.toLowerCase();
            const subsystem = row.cells[2].textContent.toLowerCase();
            const substrates = row.cells[3].textContent.toLowerCase();
            const products = row.cells[4].textContent.toLowerCase();

            if (name.includes(query) || subsystem.includes(query) || substrates.includes(query) || products.includes(query)) {
                row.style.display = '';
            } else {
                row.style.display = 'none';
            }
        });
    });
});
