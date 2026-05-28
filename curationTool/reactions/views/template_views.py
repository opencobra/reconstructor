"""
Views for managing reaction templates in the Django application.

This module provides API endpoints for creating, retrieving, listing, 
sharing, updating, and deleting reaction templates. Templates are associated 
with users and can be either user-specific or default templates.

Available Endpoints:
    - create_template: Create a new reaction template.
    - get_rxn_template: Retrieve a reaction template by name.
    - list_templates: List all available reaction templates for a user.
    - share_template: Share a template with another user.
    - update_template: Update an existing template's name and description.
    - delete_template: Delete a user-created reaction template.
"""
import json
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from reactions.models import ReactionTemplate, User
from reactions.utils.reaction_inputs import (
    normalize_reaction_side,
    reaction_has_any_side,
)

@csrf_exempt
def create_template(request):
    """
    Create a new reaction template for a user.

    Process:
        - Validates user authentication.
        - Ensures template name is unique for the user.
        - Extracts form data including substrates, products, and metadata.
        - Saves the new reaction template.

    Parameters:
        request (HttpRequest): The HTTP request containing form data.

    Returns:
        JsonResponse:
            - Success: A success message upon template creation.
            - Error: If the user is not authenticated or template name already exists.
    """
    if request.method == 'POST':
        try:
            # Parse the form data
            form_data = request.POST
            user_id = form_data.get('userID')
            template_name = form_data.get('template_name')
            # Ensure user is authenticated
            if not user_id:
                return JsonResponse(
                    {'status': 'error', 'message': 'User not authenticated.'}, status=403)

            user = User.objects.get(id=user_id)

            # Check if template name already exists
            if (ReactionTemplate.objects.filter(name=template_name, user=user).exists() or
                    ReactionTemplate.objects.filter(is_default=True, name=template_name).exists()):
                return JsonResponse(
                    {'status': 'error',
                    'message': f'You already have a template named {template_name}.'},
                    status=400)

            # Collect the fields for the template. Empty rows are ignored so
            # templates can intentionally represent one-sided reactions.
            substrates_side = normalize_reaction_side(
                form_data.getlist('substrates'),
                form_data.getlist('substrates_type'),
                form_data.getlist('subs_sch'),
                form_data.getlist('subs_comps'),
                uploaded_file_count=len(request.FILES.getlist('substrates')),
            )
            products_side = normalize_reaction_side(
                form_data.getlist('products'),
                form_data.getlist('products_type'),
                form_data.getlist('prod_sch'),
                form_data.getlist('prod_comps'),
                uploaded_file_count=len(request.FILES.getlist('products')),
            )
            substrates = substrates_side['metabolites']
            substrates_types = substrates_side['types']
            subs_comps = substrates_side['compartments']
            subs_sch = substrates_side['stoichiometries']
            products = products_side['metabolites']
            products_types = products_side['types']
            prods_comps = products_side['compartments']
            prods_sch = products_side['stoichiometries']

            if not reaction_has_any_side(substrates, products):
                return JsonResponse(
                    {
                        'status': 'error',
                        'message': 'Enter at least one substrate or one product before creating a template.'
                    },
                    status=400,
                )

            direction = form_data.get('direction', 'forward')
            subsystem = form_data.get('subsystem', 'undefined')
            organs = json.loads(form_data.get('organs', '[]'))

            description = form_data.get('description', '')
            # Create and save the new template
            template = ReactionTemplate.objects.create(
                name=template_name,
                user=user,
                is_default=False,
                description=description,
                substrates=','.join(substrates),
                substrates_types=','.join(substrates_types),
                subs_comps=','.join(subs_comps),
                subs_sch=','.join(subs_sch),
                products=','.join(products),
                products_types=','.join(products_types),
                prods_comps=','.join(prods_comps),
                prods_sch=','.join(prods_sch),
                direction=direction,
                subsystem=subsystem,
                Organs=','.join(organs)
            )

            user.templates.add(template)

            return JsonResponse({'status': 'success',
                                 'message': 'Template created successfully.'})
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=500)
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=400)


@csrf_exempt
def get_rxn_template(request):
    """
    Retrieve a reaction template by name.

    Process:
        - Fetches a user-specific template first.
        - If not found, retrieves the first available template.
        - Returns template details including substrates, products, and metadata.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing JSON data with `reaction_type` and `userID`.

    Returns:
        JsonResponse:
            - Success: Template details if found.
            - Error: If the template is not found.
    """
    if request.method == 'POST':
        try:
            # Parse the request body
            data = json.loads(request.body)
            reaction_type = data.get('reaction_type', '')
            user_id = data.get('userID')

            # Try to fetch a user-specific template
            template = ReactionTemplate.objects.filter(
                name=reaction_type, user=user_id
            ).first()

            # If no user-specific template exists, fetch the first available one
            if not template:
                template = ReactionTemplate.objects.filter(name=reaction_type).first()

            # If no template is found, return an error response
            if not template:
                return JsonResponse({'error': 'Template not found'}, status=404)

            # Prepare the response data
            substrates = template.substrates.split(',') if template.substrates else []
            subs_sch = (
                template.subs_sch.split(',')
                if template.subs_sch
                else ['1'] * len(substrates)
            )
            subs_comps = (
                template.subs_comps.split(',')
                if template.subs_comps
                else ['-'] * len(substrates)
            )
            subs_types = (
                template.substrates_types.split(',')
                if template.substrates_types
                else ['vmh'] * len(substrates)
            )

            products = template.products.split(',') if template.products else []
            prod_sch = (
                template.prods_sch.split(',')
                if template.prods_sch
                else ['1'] * len(products)
            )
            prods_comps = (
                template.prods_comps.split(',')
                if template.prods_comps
                else ['-'] * len(products)
            )
            prods_types = (
                template.products_types.split(',')
                if template.products_types
                else ['vmh'] * len(products)
            )

            result = {
                'substrates': substrates,
                'subs_sch': subs_sch,
                'subs_comps': subs_comps,
                'subs_types': subs_types,
                'products': products,
                'prod_sch': prod_sch,
                'prods_comps': prods_comps,
                'prods_types': prods_types,
                'direction': template.direction,
                'subsystem': template.subsystem,
                'description': template.description,
                'name': template.name,
                'Organs': template.Organs.split(',') if template.Organs else [],
            }

            return JsonResponse(result)

        except ReactionTemplate.DoesNotExist:
            return JsonResponse({'error': 'Template not found'}, status=404)
        except KeyError:
            return JsonResponse({'error': 'Invalid template structure'}, status=400)

    return JsonResponse({'error': 'Invalid request method'}, status=400)

@csrf_exempt
def list_templates(request):
    """
    Retrieve a list of reaction templates available to a user.

    Process:
        - Fetches default templates available to all users.
        - If a user ID is provided, retrieves the user's saved templates.
        - Merges and sorts the list of templates.

    Parameters:
        request (HttpRequest): The HTTP request containing the `userID`.

    Returns:
        JsonResponse:
            - Success: A list of available templates.
            - Error: If the request method is invalid.
    """
    if request.method != 'POST':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=400)

    default_templates = ReactionTemplate.objects.filter(
        is_default=True).values_list('name', flat=True)

    user_id = json.loads(request.body).get('userID')

    if user_id:
        user = User.objects.get(id=user_id)
        user_templates = user.templates.values_list('name', flat=True)
        templates = list(user_templates) + list(default_templates)
    else:
        templates = list(default_templates)
    templates.sort()  # Optionally sort the combined list

    return JsonResponse({'templates': templates})


@csrf_exempt
@csrf_exempt
def share_template(request):
    """
    Share one or more reaction templates with another user.

    Process:
        - Validates the sender, recipient, and provided template names.
        - Copies the shared templates for the recipient user.
        - Ensures unique names for the shared templates.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing `userID`, `share_with_user`, and `template_names`.

    Returns:
        JsonResponse:
            - Success: Confirmation that templates were shared.
            - Error: If the sender, recipient, or templates are invalid.
    """
    if request.method != 'POST':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=400
        )

    data = json.loads(request.body)
    user_id = data.get('userID')
    share_with_username = data.get('share_with_user')
    template_names = data.get('template_names', [])

    if not user_id or not share_with_username or not template_names:
        return JsonResponse(
            {'status': 'error', 'message': 'Missing required fields.'}, status=400
        )

    try:
        sharing_user = User.objects.get(id=user_id)
    except User.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'Sharing user not found.'}, status=404
        )

    try:
        recipient_user = User.objects.get(name=share_with_username)
    except User.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'Recipient user not found.'}, status=404
        )

    # Retrieve templates to share
    templates_to_share = sharing_user.templates.filter(name__in=template_names)

    for template in templates_to_share:
        # Generate a unique name for the recipient’s copy
        new_name = template.name
        if recipient_user.templates.filter(name=new_name).exists():
            base_name = template.name
            counter = 1
            new_name = f"{base_name}{counter}"
            while recipient_user.templates.filter(name=new_name).exists():
                counter += 1
                new_name = f"{base_name}{counter}"

        # Create and save a copy for the recipient
        new_template = ReactionTemplate.objects.create(
            name=new_name,
            user=recipient_user,
            is_default=False,
            substrates=template.substrates,
            products=template.products,
            direction=template.direction,
            substrates_types=template.substrates_types,
            products_types=template.products_types,
            subsystem=template.subsystem,
            subs_comps=template.subs_comps,
            prods_comps=template.prods_comps,
            subs_sch=template.subs_sch,
            prods_sch=template.prods_sch,
            Organs=template.Organs,
            description=template.description,
        )
        recipient_user.templates.add(new_template)
        recipient_user.save()

    return JsonResponse(
        {'status': 'success', 'message': 'Templates shared successfully.'}, status=200
    )


@csrf_exempt
@csrf_exempt
def update_template(request):
    """
    Update the name or description of a reaction template.

    Process:
        - Validates user ownership of the template.
        - Checks for name conflicts.
        - Updates the existing template if owned by the user.
        - Creates a personal copy for the user if the template is a default.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing `old_name`, `new_name`, `userID`, and `description`.

    Returns:
        JsonResponse:
            - Success: Confirmation of template update.
            - Error: If the template is not found or name conflicts exist.
    """
    if request.method != 'POST':
        return JsonResponse({'status': 'error', 'message': 'Invalid request method.'}, status=400)

    try:
        # Parse request body
        data = json.loads(request.body)
        old_name = data.get('old_name')
        new_name = data.get('new_name')
        user_id = data.get('userID')
        description = data.get('description', '')

        if not user_id or not old_name or not new_name:
            return JsonResponse(
                {'status': 'error', 'message': 'Missing required fields.'}, status=400
            )

        user = User.objects.get(id=user_id)

        # Check if the template exists in the user's saved templates
        template = user.templates.filter(name=old_name).first()

        if template and template.user and template.user.id == user.id:
            # Ensure new name is not already used in user's templates or default templates
            if (
                user.templates.filter(name=new_name).exclude(id=template.id).exists()
                or ReactionTemplate.objects.filter(is_default=True, name=new_name).exists()
            ):
                return JsonResponse(
                    {'status': 'error',
                     'message': f'A template named {new_name} already exists.'}, 
                     status=400
                )
            # Update template details
            template.name = new_name
            template.description = description
            template.save()

            return JsonResponse({'status': 'success', 'message': 'Template updated successfully.'})

        # Otherwise, check if the template is a default or belongs to another user
        template = ReactionTemplate.objects.filter(name=old_name).first()

        if template and template.is_default:
            return JsonResponse(
                {'status': 'error', 'message': 'Cannot update default templates.'}, status=403
            )

        if not template:
            return JsonResponse(
                {'status': 'error', 'message': 'Template not found in your list.'}, status=404
            )

        # Ensure the new name does not conflict with other templates
        if (
            user.templates.filter(name=new_name).exists()
            or ReactionTemplate.objects.filter(is_default=True, name=new_name).exists()
        ):
            return JsonResponse(
                {'status': 'error',
                 'message': f'A template named {new_name} already exists.'}, 
                 status=400
            )

        # Create a copy for the user and update it
        template.pk = None  # Create a new record
        template.user = user
        template.name = new_name
        template.description = description
        template.is_default = False
        template.save()
        user.templates.add(template)

        return JsonResponse(
            {'status': 'success', 'message': 'Template copied and updated successfully.'}
        )

    except User.DoesNotExist:
        return JsonResponse({'status': 'error', 'message': 'User not found.'}, status=404)
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)

@csrf_exempt
def delete_template(request):
    """
    Delete a reaction template.

    Process:
        - Validates user authentication.
        - Ensures the template is not a default template.
        - Removes the template from the user's saved list.
        - Deletes the template from the database if it has no other owners.

    Parameters:
        request (HttpRequest): 
            The HTTP request containing `template_name` and `userID`.

    Returns:
        JsonResponse:
            - Success: Confirmation of template deletion.
            - Error: If the template is a default or not found.
    """
    if request.method != 'POST':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=400
        )

    try:
        # Parse request data
        data = json.loads(request.body)
        template_name = data.get('template_name')
        user_id = data.get('userID')

        if not user_id:
            return JsonResponse(
                {'status': 'error', 'message': 'Login required.'}, status=403
            )

        user = User.objects.get(id=user_id)

        # Find the template in the user's saved templates
        template = user.templates.filter(name=template_name).first()

        if not template:
            # Check if the template is a default template
            template = ReactionTemplate.objects.filter(
                name=template_name, is_default=True
            ).first()

            if template:
                return JsonResponse(
                    {'status': 'error', 'message': 'Cannot delete default templates.'},
                    status=403
                )

            return JsonResponse(
                {'status': 'error', 'message': 'Template not found in your list.'},
                status=404
            )

        # Remove the template from the user's saved list
        user.templates.remove(template)

        # Delete the template from the database if no other users are associated with it
        if template.created_by_users.count() == 0:
            template.delete()

        return JsonResponse(
            {'status': 'success', 'message': 'Template deleted successfully.'}
        )

    except User.DoesNotExist:
        return JsonResponse(
            {'status': 'error', 'message': 'User not found.'}, status=404
        )
    except Exception as e:
        return JsonResponse(
            {'status': 'error', 'message': str(e)}, status=500
        )
