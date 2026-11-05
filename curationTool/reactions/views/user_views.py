"""
User authentication and reaction management views.

This module provides API endpoints for user-related actions, including:
    - User authentication (login and registration)
    - Session management
    - Reaction saving and flagging

Functions:
    - get_user: Authenticates users and starts a session.
    - register_user: Handles new user registration.
    - validate_user_ID: Checks if a user ID is valid.
    - set_session_user: Stores a user ID in the session.
    - save_user_reaction: Saves a user's reaction with optional metadata.
"""

import json
from django.http import JsonResponse
from django.shortcuts import render
from reactions.models import ApiKey, User, Reaction, Flag
from reactions.services import api_keys as api_key_service

def get_user(request):
    """
    Authenticate a user and initiate a session.

    Process:
        - Retrieves username and password from request.
        - Validates user credentials.
        - If valid, stores user ID in session.

    Parameters:
        request (HttpRequest): The HTTP request containing POST data with 
                               'username' and 'password'.

    Returns:
        JsonResponse:
            - Success: Returns user details and session ID.
            - Error: If credentials are incorrect or user does not exist.
    """
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')
    try:
        user = User.objects.get(name=username)
        if user.check_password(password):
            # Store the user ID in the session
            request.session['userID'] = user.pk
            return JsonResponse(
                {'status': 'success', 'userName': user.name, 'userID': user.pk})

        return JsonResponse(
            {'status': 'error', 'message': 'Invalid password'})
    except User.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'User does not exist'})


def register_user(request):
    """
    Register a new user.

    Process:
        - Retrieves username, password, email, ORCID ID, and full name from request.
        - Ensures required fields are provided and username is unique.
        - Creates and saves a new user in the database.
        - Stores user ID in session.

    Parameters:
        request (HttpRequest): The HTTP request containing POST data with 
                               'username', 'password', 'email', 'orchid_id', and 'full_name'.

    Returns:
        JsonResponse:
            - Success: Returns user details and session ID.
            - Error: If username already exists or required fields are missing.
    """
    username = request.POST.get('username', '')
    password = request.POST.get('password', '')
    email = request.POST.get('email', '')
    orchid_id = request.POST.get('orchid_id', '')
    full_name = request.POST.get('full_name', '')   # New field
    users = User.objects.all()
    usernames = [user.name for user in users]
    if username in usernames:
        return JsonResponse(
            {'status': 'error', 'message': 'Username already exists'})
    if not username or not password or not email or not full_name:
        return JsonResponse(
            {'status': 'error', 'message': 'Ysername, password, email, and full name are required'})

    try:
        if User.objects.filter(name=username).exists():
            return JsonResponse(
                {'status': 'error', 'message': 'Username already exists'})

        user = User(
            name=username,
            email=email,
            password=password,  # Hash the password before saving
            orchid_id=orchid_id,
            full_name=full_name
        )
        user.save()

        request.session['userID'] = user.pk  # Store the user ID in the session
        return JsonResponse(
            {'status': 'success', 'userName': user.name, 'userID': user.pk})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)})


def validate_user_ID(user_id):
    """
    Validate if a user ID exists in the database.

    Process:
        - Attempts to retrieve a user by primary key.
        - Returns the user object if found.

    Parameters:
        user_id (int): The primary key of the user.

    Returns:
        User or None:
            - User: If the ID is valid.
            - None: If no matching user exists.
    """
    try:
        user = User.objects.get(pk=user_id)
        return user
    except User.DoesNotExist:
        return None


def set_session_user(request):
    """
    Store a user ID in the session.

    Process:
        - Extracts 'userID' from request body.
        - Saves it in the session.
        - Verifies session storage and user validity.

    Parameters:
        request (HttpRequest): The HTTP request containing JSON body with 'userID'.

    Returns:
        JsonResponse:
            - Success: If user ID is successfully set in session.
            - Error: If validation fails.
    """
    req_body = json.loads(request.body)
    user_id = int(req_body.get('userID'))
    request.session['userID'] = user_id
    if request.session.get('userID') == user_id and validate_user_ID(user_id):
        return JsonResponse(
            {'status': 'success', 'message': 'User set in session'})
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'User not set in session'})


def save_user_reaction(request):
    """
    Save a user's reaction with optional metadata.

    Process:
        - Retrieves reaction ID, user ID, and additional details from request.
        - Validates the reaction and user.
        - Updates reaction metadata (short name, description, flags, confidence score).
        - Associates the reaction with the user.

    Parameters:
        request (HttpRequest): The HTTP request containing POST data with 
                               'reaction_id', 'userID', 'short_name', 
                               'description', 'flag_name', 'flag_color', and 'confidence_score'.

    Returns:
        JsonResponse:
            - Success: If the reaction is saved successfully.
            - Error: If user or reaction is invalid or if request is incorrect.
    """
    if request.method == 'POST':
        reaction_id = request.POST.get('reaction_id')
        userID = request.POST.get('userID')
        short_name = request.POST.get('short_name')
        description = request.POST.get('description', '')  # Default to empty string
        flag_name = request.POST.get('flag_name')
        flag_color = request.POST.get('flag_color')
        confidence_score = request.POST.get('confidence_score', None)  # Get confidence score

        try:
            reaction = Reaction.objects.get(pk=reaction_id)
        except Reaction.DoesNotExist:
            return JsonResponse({'status': 'error', 'message': 'Invalid reaction'})

        user = validate_user_ID(userID)

        if user and reaction:
            reaction.short_name = short_name
            reaction.description = description  # Save the description

            # Validate confidence score (allow null)
            if confidence_score not in ["1", "2", "3", "4", None]:
                return JsonResponse({'status': 'error', 'message': 'Invalid confidence score'})
            if flag_name.strip() not in ['None', 'Choose a flag'] and flag_color != 'null':
                flag, _ = Flag.objects.get_or_create(name_flag=flag_name, color=flag_color)
                reaction.flags.add(flag)
            reaction.confidence_score = confidence_score
            reaction.save()
            user.saved_reactions.add(reaction)

            return JsonResponse({'status': 'success'})

        return JsonResponse({'status': 'error', 'message': 'Invalid user or reaction'})

    return JsonResponse({'status': 'error', 'message': 'Invalid request'})


def _session_user(request):
    """Return the signed-in User from the session, or None."""
    user_id = request.session.get('userID')
    if not user_id:
        return None
    try:
        return User.objects.get(pk=user_id)
    except User.DoesNotExist:
        return None


def api_key_page(request):
    """
    Render the API key management page (feature 002, FR-003).

    Shows the key prefix and creation date when a key exists; the raw key is
    only ever returned by api_key_generate, once.
    """
    user = _session_user(request)
    key = None
    if user:
        key = ApiKey.objects.filter(user=user).first()
    return render(request, 'reactions/api_key.html', {
        'session_user': user,
        'key_prefix': key.prefix if key else None,
        'key_created_at': key.created_at if key else None,
    })


def api_key_generate(request):
    """
    Generate (or rotate) the signed-in user's API key.

    Returns the full key exactly once; only its hash is stored. Rotation
    invalidates the previous key immediately.
    """
    if request.method != 'POST':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request'}, status=400)
    user = _session_user(request)
    if not user:
        return JsonResponse(
            {'status': 'error', 'message': 'Sign in to manage your API key.'},
            status=401)
    raw_key = api_key_service.generate_key(user)
    return JsonResponse({'status': 'success', 'api_key': raw_key})
