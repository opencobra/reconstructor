
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
    pass