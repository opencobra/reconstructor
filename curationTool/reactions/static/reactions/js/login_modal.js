// Show modal when login button is clicked
document.getElementById('loginButton').addEventListener('click', function () {
  if (this.textContent === 'Log out') {
    // Clear the session storage
    sessionStorage.clear();
    // Reset the user display
    document.getElementById('userDisplay').innerHTML = '<i class="icon user"></i> User : None';
    window.location.href = window.location.origin;
    this.textContent = 'Log in';
    return}

  else{$('#loginModal').modal('show');}
});

// Close modal when close icon is clicked
document.querySelector('.ui.modal .close').addEventListener('click', function () {
  $('#loginModal').modal('hide');
});

// Handle form submission
document.getElementById('loginForm').addEventListener('submit', function (event) {
  event.preventDefault(); // Prevent the form from submitting the traditional way
  var userName = document.getElementById('username').value;
  var password = document.getElementById('password').value;

  if (userName && password) {
      // Prepare the data to be sent in the POST request
      const data = new FormData();
      data.append('username', userName);
      data.append('password', password);
        console.log(data);
      // Adjust the fetch call to use POST
      fetch(`${getUser}`, {
          method: 'POST',
          headers: {
              'X-Requested-With': 'XMLHttpRequest',
              'X-CSRFToken': csrfToken,
          },
          body: data
      })
      .then(response => response.json())
      .then(data => {console.log(data,"success");

          if (data.status === 'success') {
              sessionStorage.setItem('userID', data.userID);
              sessionStorage.setItem('userName', data.userName);
              document.getElementById('userDisplay').innerHTML = `<i class="icon user"></i> User: ${data.userName}`;
                document.getElementById('loginButton').textContent = 'Log out';
                setLoggedInStatusBasedOnUrl();
              // Hide the modal after successful login
              $('#loginModal').modal('hide');
              //todo
            fetchAvailableReactions(data.userID);
          } else {
              alert(data.message);
          }
      })
      .catch(error => {
          console.error('Error:', error);
          alert("Failed to fetch user details.");
      });
  }

});




async function fetchAvailableReactions(userId) {
  try {
      const response = await fetch('available_reactions', {
          method: 'POST',
          headers: {
              'X-Requested-With': 'XMLHttpRequest',
              'X-CSRFToken': csrfToken
          },
          body: JSON.stringify({ user_id: userId })
      });

      if (!response.ok) {
          throw new Error('Network response was not ok');
      }

      const data = await response.json();

      if (data.error) {
          console.error(data.error);
          return;
      }

      handleLoginResponse(data);
  } catch (error) {
      console.error('There was a problem with the fetch operation:', error);
  }
}

function handleLoginResponse(data) {
  const loginModal = document.getElementById('loginModal');
  
  if (!data.available_reaction_ids || data.available_reaction_ids.length === 0) {
      loginModal.style.display = 'none';
  } else {
      const lastReactionId = data.last_index;
      window.location.href = window.location.origin + "/?reaction_id=" + lastReactionId;
  }
}

