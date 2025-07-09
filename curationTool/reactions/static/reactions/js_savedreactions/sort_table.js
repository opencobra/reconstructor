document.addEventListener('DOMContentLoaded', function () {
  const table = document.getElementById('reactionList');
  const tbody = table.querySelector('tbody');
  const originalRows = Array.from(tbody.querySelectorAll('tr')); // Store initial order

  const getValue = (row, key) => {
    switch (key) {
      case 'flag':
        const icon = row.querySelector('.flag-icon');
        return icon?.getAttribute('data-flag-name')?.toLowerCase() || 'zzz'; // 'zzz' pushes nulls to bottom
      case 'subsystem':
        return row.children[3]?.textContent.trim().toLowerCase() || 'zzz';
      case 'substrates':
        return row.children[4]?.textContent.trim().toLowerCase() || 'zzz';
      case 'products':
        return row.children[6]?.textContent.trim().toLowerCase() || 'zzz';
      default:
        return '';
    }
  };

  document.getElementById('sortBy').addEventListener('change', () => {
    const sortKey = document.getElementById('sortBy').value;

    // Reset to original order if no sort selected
    if (!sortKey) {
      tbody.innerHTML = '';
      originalRows.forEach(row => tbody.appendChild(row.cloneNode(true)));
      return;
    }

    const sortedRows = [...tbody.querySelectorAll('tr')].sort((a, b) => {
      const aVal = getValue(a, sortKey);
      const bVal = getValue(b, sortKey);
      return aVal.localeCompare(bVal);
    });

    tbody.innerHTML = '';
    sortedRows.forEach(row => tbody.appendChild(row));
  });
});