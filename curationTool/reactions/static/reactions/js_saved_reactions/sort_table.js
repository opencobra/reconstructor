document.addEventListener('DOMContentLoaded', function () {
  const table = document.getElementById('reactionList');
  const tbody = table.querySelector('tbody');
  const originalRows = Array.from(tbody.querySelectorAll('tr'));
  const sortIcon = document.getElementById('sortOrderToggle');

  // false = descending (default: latest first for date, Z-A for others)
  // true  = ascending  (oldest first for date, A-Z for others)
  let sortAscending = false;

  const getValue = (row, key) => {
    switch (key) {
      case 'flag': {
        const icon = row.querySelector('.flag-icon');
        return icon?.getAttribute('data-flag-name')?.toLowerCase() || 'zzz';
      }
      // Column order: Select(0) Name(1) CS(2) Description(3) Subsystem(4)
      // Substrates(5) Direction(6) Products(7) ...
      case 'subsystem':
        return row.children[4]?.textContent.trim().toLowerCase() || 'zzz';
      case 'substrates':
        return row.children[5]?.textContent.trim().toLowerCase() || 'zzz';
      case 'products':
        return row.children[7]?.textContent.trim().toLowerCase() || 'zzz';
      default:
        return '';
    }
  };

  function applySort() {
    const sortKey = document.getElementById('sortBy').value;
    let rows;

    if (!sortKey) {
      rows = sortAscending ? [...originalRows] : [...originalRows].reverse();
    } else {
      rows = [...originalRows].sort((a, b) => {
        const cmp = getValue(a, sortKey).localeCompare(getValue(b, sortKey));
        return sortAscending ? cmp : -cmp;
      });
    }

    rows.forEach((row) => tbody.appendChild(row));

    if (typeof rebindCheckboxListeners === 'function') rebindCheckboxListeners();

    sortIcon.classList.toggle('sort-asc', sortAscending);
  }

  // Default: latest first
  applySort();

  document.getElementById('sortBy').addEventListener('change', () => {
    sortAscending = false;
    applySort();
  });

  sortIcon.addEventListener('click', () => {
    sortAscending = !sortAscending;
    applySort();
  });
});
