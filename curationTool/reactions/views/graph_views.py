from django.http import JsonResponse
from reactions.models import User, Reaction
from reactions.views.user_views import validate_user_ID
import json
def get_graph_info(request):
    """
    View to get graph information for a user's saved reactions. 
    Returns a JSON response with the graph data.
    1. Authenticate the user.
    2. Retrieve the user's saved reactions.
    3 (optional) exclude some reactions based on query parameters.
    4. Return the graph data as JSON:
        - nodes: metabolites (id, name, formula)
        - edges: reactions (id, name, substrates, products)
    """
    userID = request.POST.get('userID')
    user = validate_user_ID(userID)
    if not user:
        return JsonResponse({'status': 'error', 'message': 'Invalid user'})
    reactions = user.saved_reactions.all()
    nodes = {}
    edges = []

    for reaction in reactions:
        substrates = json.loads(reaction.substrates)
        products = json.loads(reaction.products)
        substrate_names = json.loads(reaction.substrate_names)
        product_names = json.loads(reaction.product_names)
        substrates_types = json.loads(reaction.substrates_types)
        products_types = json.loads(reaction.products_types)
        if "Saved" in substrates_types or "Saved" in products_types:
            # deal with saved metabolites here
            pass
        for substrate in substrates:
            if substrate not in nodes:
                nodes[substrate] = {
                    'id': substrate,
                    'vmh_abbreviation': substrate,
                    'name': substrate,
                }
        for product in products:
            if product not in nodes:
                nodes[product] = {
                    'id': product,
                    'vmh_abbreviation': product,
                    'name': product,
                }
        edges.append({
            'from': substrates,
            'to': products,
            'id': reaction.id,
            'name': reaction.short_name,
        })
    return JsonResponse({
        'status': 'success',
        'data': {
            'nodes': nodes,
            'edges': edges,
        }
    })