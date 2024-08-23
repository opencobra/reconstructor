// Show modal when register link is clicked
document.getElementById('register-link').addEventListener('click', function (event) {
    event.preventDefault(); // Prevent the default link behavior
    $('#loginModal').modal('hide');
    $('#registerModal').modal('show');
});

// Close modal when close icon is clicked
document.querySelectorAll('.ui.modal .close').forEach(button => {
    button.addEventListener('click', function () {
        $(this).closest('.ui.modal').modal('hide');
    });
});

// Handle form submission
document.getElementById('registerForm').addEventListener('submit', function (event) {
    event.preventDefault(); // Prevent the form from submitting the traditional way
    var userName = document.getElementById('username-reg').value;
    var password = document.getElementById('password-reg').value;
    var orchid_id = document.getElementById('orchidid').value;
    var email = document.getElementById('email').value;

    if (userName && password && orchid_id && email) {
        // Prepare the data to be sent in the POST request
        const data = new FormData();
        data.append('username', userName);
        data.append('password', password);
        data.append('orchid_id', orchid_id);
        data.append('email', email);

        // Adjust the fetch call to use POST
        fetch(`${regUser}`, {
            method: 'POST',
            headers: {
                'X-Requested-With': 'XMLHttpRequest',
                'X-CSRFToken': csrfToken,
            },
            body: data
        })
        .then(response => response.json())
        .then(data => {
            if (data.status === 'success') {
                sessionStorage.setItem('userID', data.userID);
                sessionStorage.setItem('userName', data.userName);
                document.getElementById('userDisplay').textContent = `User: ${data.userName}`;
                document.getElementById('loginButton').textContent = 'Log out';
                // Hide the modal after successful registration
                $('#registerModal').modal('hide');
                setLoggedInStatusBasedOnUrl();  
            } else {
                alert(data.message);
            }
        })
        .catch(error => {
            console.error('Error:', error);
            alert("Failed to register user.");
        });
    }
});

// Password validation logic
document.getElementById('password-reg').addEventListener('input', function () {
    const password = this.value;
    const passwordCriteria = {
        length: password.length >= 8,
        uppercase: /[A-Z]/.test(password),
        lowercase: /[a-z]/.test(password),
        number: /[0-9]/.test(password),
        specialChar: /[!@#$%^&*]/.test(password)
    };

    Object.keys(passwordCriteria).forEach(key => {
        const element = document.getElementById(`${key}-criteria`);
        if (passwordCriteria[key]) {
            element.classList.add('valid');
            element.classList.remove('invalid');
        } else {
            element.classList.add('invalid');
            element.classList.remove('valid');
        }
    });
});