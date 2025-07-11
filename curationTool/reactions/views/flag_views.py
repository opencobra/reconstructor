"""
This module provides API endpoints for managing user 
flags and associating them with saved reactions.

It includes:
- `get_user_flags`: Retrieves all flags associated with a given user.
- `add_flag`: Allows users to create and assign flags.
- `save_flags_in_saved_reactions`: Saves flags to specific user reactions.
"""

import json

from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import get_object_or_404

from reactions.models import Flag, User, Reaction

def get_user_flags(request, user_id): #  pylint: disable=unused-argument
    """
    Retrieve all flags associated with a specific user.

    Parameters:
        request (HttpRequest): The incoming request object.
        user_id (int): The ID of the user whose flags should be retrieved.

    Returns:
        JsonResponse:
            - Success: A list of flags containing `id`, `name_flag`, and `color`.
            - Error: An error message if the request fails.
    """
    try:
        user = User.objects.get(pk=user_id) #  pylint: disable=no-member
        flags = user.flags.all()
        flags_data = [{'id': flag.id,
                       'name_flag': flag.name_flag,
                       'color': flag.color} for flag in flags]
        return JsonResponse({'status': 'success', 'flags': flags_data})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)


@csrf_exempt
def add_flag(request):
    """
    Create and assign a flag to a user.

    Process:
        - Extracts `user_id`, `name_flag`, and `color` from the request body.
        - Checks if the user exists.
        - Creates a new flag or retrieves an existing one with the same name and color.

    Expects:
        JSON body with:
        - `user_id` (int): The user's ID.
        - `name_flag` (str): The name of the flag.
        - `color` (str): The flag's color.

    Returns:
        JsonResponse:
            - Success: Contains flag details (`id`, `name_flag`, `color`).
            - Error: If user doesn't exist or required data is missing.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('user_id')
            flag_name = data.get('name_flag')
            flag_color = data.get('color')
            user = User.objects.get(pk=user_id) #  pylint: disable=no-member
        except User.DoesNotExist: #  pylint: disable=no-member
            return JsonResponse({'status': 'error', 'message': 'Invalid user'})

        if flag_name and flag_color:
            flag_name = flag_name.strip()
            flag, _ = Flag.objects.get_or_create( #  pylint: disable=no-member
                name_flag=flag_name, color=flag_color, user=user)
            return JsonResponse({
                'status': 'success',
                'message': 'Flag added successfully',
                'flag': {'id': flag.id, 'name_flag': flag.name_flag, 'color': flag.color}
            })
        return JsonResponse({'status': 'error',
                                'message': 'Flag name and color are required'})
    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method'})
    

@csrf_exempt
def remove_flag(request):
    if request.method != 'POST':
        return JsonResponse({'status': 'error',
                             'message': 'Invalid request'}, status=405)

    try:
        data         = json.loads(request.body)
        reaction_id  = int(data['reaction_id'])
        flag_id      = int(data['flag_id'])          # sent by JS
        user_id      = int(data['user_id'])          # still useful for auth

        # 1. make sure the reaction is one the user owns / saved
        user = User.objects.get(pk=user_id)
        reaction = user.saved_reactions.get(pk=reaction_id)  # ← key line
        #   (Use created_by_users.get(pk=…) if that is your ownership rule)

        # 2. fetch the flag regardless of who created it
        flag = Flag.objects.get(pk=flag_id)

        # 3. detach the flag from *this* reaction
        reaction.flags.remove(flag)

        return JsonResponse({'status': 'success'})

    except (User.DoesNotExist, Reaction.DoesNotExist):
        return JsonResponse({'status': 'error',
                             'message': 'Reaction not found for this user'},
                            status=404)

    except Flag.DoesNotExist:
        return JsonResponse({'status': 'error',
                             'message': 'Flag not found'}, status=404)

    except Exception as e:
        return JsonResponse({'status': 'error',
                             'message': str(e)}, status=400)


@csrf_exempt
def save_flags_in_saved_reactions(request):
    """
    Save flags to a list of user reactions.

    Process:
        - Extracts `userID`, `reaction_ids`, `flag_name`, and `flag_color` from the request body.
        - Retrieves the user and associated flag.
        - Associates the flag with the specified reactions.

    Expects:
        JSON body with:
        - `userID` (int): The user's ID.
        - `reaction_ids` (list): List of reaction IDs to associate the flag with.
        - `flag_name` (str): Name of the flag.
        - `flag_color` (str): Color of the flag.

    Returns:
        JsonResponse:
            - Success: If flags are added successfully.
            - Error: If the user, flag, or reactions are not found.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            user_id = data.get('userID')
            reaction_ids = data.get('reaction_ids', [])
            flag_name = data.get('flag_name')
            flag_color = data.get('flag_color')

            user = User.objects.get(pk=user_id) #  pylint: disable=no-member
            flag = Flag.objects.get( #  pylint: disable=no-member
                name_flag=flag_name, color=flag_color, user=user)
            for reaction_id in reaction_ids:
                reaction = get_object_or_404(Reaction, pk=reaction_id)
                reaction.flags.add(flag)

            return JsonResponse({'status': 'success'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)})
    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method'})
