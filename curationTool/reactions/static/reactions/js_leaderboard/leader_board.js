function loadLeaderBoardData() {
  fetch('/user-reactions-vmh/')
    .then(response => response.json())
    .then(usersData => {
      // Sort usersData by reactions saved (descending) and limit to top 10
      usersData.sort((a, b) => b.saved - a.saved);
      usersData = usersData.slice(0, 10);

      // Populate the table
      const tableBody = document.getElementById('leaderboard-table-body');
      if (tableBody) {
        tableBody.innerHTML = ''; // Clear existing rows

        usersData.forEach((user, index) => {
          const row = document.createElement('tr');
          row.classList.add('clickable-row');
          row.setAttribute('data-name', user.full_name);
          row.setAttribute('data-saved', user.saved);
          row.setAttribute('data-added', user.added);

          row.innerHTML = `
            <th>${index + 1}</th>
            <td>${user.full_name}</td>
            <td>${user.saved}</td>
            <td>${user.added}</td>
          `;
          tableBody.appendChild(row);
        });

        // Create the chart
        const ctx = document.getElementById('reactionChart');
        if (ctx) {
          const chartContext = ctx.getContext('2d');
          if (chartContext) {
            let reactionChart = new Chart(chartContext, {
              type: 'bar',
              data: {
                labels: usersData.map(user => user.full_name),
                datasets: [
                  {
                    label: 'Reactions Saved',
                    data: usersData.map(user => user.saved),
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 1
                  },
                  {
                    label: 'Reactions Added to VMH',
                    data: usersData.map(user => user.added),
                    backgroundColor: 'rgba(153, 102, 255, 0.2)',
                    borderColor: 'rgba(153, 102, 255, 1)',
                    borderWidth: 1
                  }
                ]
              },
              options: {
                scales: {
                  x: {
                    ticks: {
                      autoSkip: false,
                    }
                  },
                  y: {
                    beginAtZero: true
                  }
                }
              }
            });

            // Event listener for row click
            document.querySelectorAll('.clickable-row').forEach(row => {
              row.addEventListener('click', () => {
                const name = row.getAttribute('data-name');
                const saved = row.getAttribute('data-saved');
                const added = row.getAttribute('data-added');

                reactionChart.data.labels = ['Reactions Saved', 'Reactions Added to VMH'];
                reactionChart.data.datasets = [{
                  label: `${name}'s Stats`,
                  data: [saved, added],
                  backgroundColor: [
                    'rgba(75, 192, 192, 0.2)',
                    'rgba(153, 102, 255, 0.2)',
                  ],
                  borderColor: [
                    'rgba(75, 192, 192, 1)',
                    'rgba(153, 102, 255, 1)',
                  ],
                  borderWidth: 1
                }];
                reactionChart.update();

                document.getElementById('chartTitle').innerText = `${name}'s Stats`;
              });

              // Event listener for row double-click to reset the chart
              row.addEventListener('dblclick', () => {
                document.getElementById('chartTitle').innerText = 'Total Stats';
                reactionChart.data.labels = usersData.map(user => user.full_name);
                reactionChart.data.datasets = [
                  {
                    label: 'Reactions Saved',
                    data: usersData.map(user => user.saved),
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 1
                  },
                  {
                    label: 'Reactions Added to VMH',
                    data: usersData.map(user => user.added),
                    backgroundColor: 'rgba(153, 102, 255, 0.2)',
                    borderColor: 'rgba(153, 102, 255, 1)',
                    borderWidth: 1
                  }
                ];
                reactionChart.update();
              });
            });

            // Initial total stats
            document.getElementById('chartTitle').innerText = 'Total Stats';
            reactionChart.data.labels = usersData.map(user => user.full_name);
            reactionChart.data.datasets = [
              {
                label: 'Reactions Saved',
                data: usersData.map(user => user.saved),
                backgroundColor: 'rgba(75, 192, 192, 0.2)',
                borderColor: 'rgba(75, 192, 192, 1)',
                borderWidth: 1
              },
              {
                label: 'Reactions Added to VMH',
                data: usersData.map(user => user.added),
                backgroundColor: 'rgba(153, 102, 255, 0.2)',
                borderColor: 'rgba(153, 102, 255, 1)',
                borderWidth: 1
              }
            ];
            reactionChart.update();
          }
        }
      } else {
        console.error('Error: The table body element was not found.');
      }
    })
    .catch(error => console.error('Error fetching user data:', error));
}

document.addEventListener('DOMContentLoaded', function() {
  const viewLeaderBoardButton = document.getElementById('view-leader-board');
  
  // Add click event for opening leaderboard if the button is present
  if (viewLeaderBoardButton) {
    viewLeaderBoardButton.addEventListener('click', () => {
      window.open(`${window.location.origin}/stats/`, '_blank');
    });
  }

  // Automatically load leaderboard data if already on the stats page
  if (window.location.pathname === '/stats/') {
    loadLeaderBoardData();
  }
});

