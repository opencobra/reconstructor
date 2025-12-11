from django.http import JsonResponse
from django.views.decorators.http import require_POST
from reactions.models import User, Reaction, SavedMetabolite
from reactions.views.user_views import validate_user_ID
import json


def _safe_json_loads(value, default=None):
    """Safely parse JSON, returning default on failure."""
    if default is None:
        default = []
    if not value:
        return default
    if isinstance(value, (list, dict)):
        return value
    try:
        return json.loads(value)
    except (json.JSONDecodeError, TypeError):
        return default


def _get_node_id(met_identifier, met_type):
    """
    Generate a standardized node ID based on metabolite type.
    
    - VMH metabolites: "vmh:{abbr}"
    - Saved metabolites: "saved:{id}"
    """
    if met_type.lower() == 'vmh':
        return f"vmh:{met_identifier}"
    elif met_type == 'Saved':
        return f"saved:{met_identifier}"
    else:
        # For other types (chebi, pubchem, etc.) that haven't been saved yet
        # This shouldn't happen in saved reactions, but handle it gracefully
        return f"other:{met_identifier}"


def _build_vmh_node(abbr, name, compartment=None, in_vmh=True, formula=None):
    """Build a node dict for a VMH metabolite."""
    return {
        'id': f"vmh:{abbr}",
        'type': 'vmh',
        'abbr': abbr,
        'name': name,
        'compartment': compartment,
        'in_vmh': in_vmh,
        'formula': formula,
        'vmh_url': f"https://www.vmh.life/#metabolite/{abbr}" if in_vmh else None,
    }


def _build_saved_node(saved_met, compartment=None):
    """Build a node dict for a SavedMetabolite."""
    return {
        'id': f"saved:{saved_met.id}",
        'type': 'saved',
        'saved_id': saved_met.id,
        'abbr': saved_met.vmh_abbr or saved_met.name[:10],
        'name': saved_met.name,
        'compartment': compartment,
        'in_vmh': False,
        'formula': saved_met.mol_formula,
        'source_type': saved_met.source_type,
        'inchi_key': saved_met.inchi_key,
    }


def _build_edge(reaction, substrate_node_ids, product_node_ids):
    """Build an edge dict for a reaction (hyperedge connecting multiple nodes)."""
    direction = reaction.direction or 'forward'
    
    # Get flags for styling
    flags = [{'id': f.pk, 'name': f.name_flag, 'color': f.color} 
             for f in reaction.flags.all()]
    
    return {
        'id': f"rxn_{reaction.id}",
        'reaction_id': reaction.id,
        'name': reaction.short_name or f"Reaction {reaction.id}",
        'description': reaction.description,
        'substrates': substrate_node_ids,
        'products': product_node_ids,
        'direction': direction,
        'reversible': direction.lower() != 'forward',
        'subsystem': reaction.subsystem,
        'balanced': _safe_json_loads(reaction.balanced_count, [None])[0],
        'confidence_score': reaction.confidence_score,
        'flags': flags,
        'vmh_found': reaction.vmh_found,
    }


@require_POST
def get_graph_info(request):
    """
    Generate graph data for visualizing a user's saved reactions as a hypergraph.
    
    Metabolites are nodes, reactions are hyperedges connecting substrate nodes to product nodes.
    
    Node standardization:
    - VMH metabolites: ID = "vmh:{abbr}"
    - Saved metabolites: ID = "saved:{id}"
    
    Returns JSON:
    {
        "status": "success",
        "data": {
            "nodes": { node_id: {node_data}, ... },
            "edges": [ {edge_data}, ... ],
            "stats": { summary statistics }
        }
    }
    """
    # ─────────────────────────────────────────────────────────────────────────
    # 1. Authenticate user
    # ─────────────────────────────────────────────────────────────────────────
    user_id = request.POST.get('userID')
    
    # Handle empty or missing userID
    if not user_id or user_id.strip() == '':
        return JsonResponse({'status': 'error', 'message': 'User ID is required'}, status=400)
    
    user = validate_user_ID(user_id)
    if not user:
        return JsonResponse({'status': 'error', 'message': 'Invalid user'}, status=401)

    # ─────────────────────────────────────────────────────────────────────────
    # 2. Fetch user's saved reactions (optionally filtered by selected IDs)
    # ─────────────────────────────────────────────────────────────────────────
    reactions = user.saved_reactions.prefetch_related('flags').all()
    
    # Check if specific reaction IDs were provided (for filtering to checked reactions)
    reaction_ids_param = request.POST.get('reactionIDs')
    selected_reaction_ids = None
    if reaction_ids_param:
        try:
            selected_reaction_ids = json.loads(reaction_ids_param)
            # Convert to integers for filtering
            selected_reaction_ids = [int(rid) for rid in selected_reaction_ids if rid]
        except (json.JSONDecodeError, ValueError, TypeError):
            selected_reaction_ids = None
    
    # If specific reactions are selected, filter to only those
    if selected_reaction_ids:
        reactions = reactions.filter(id__in=selected_reaction_ids)
    
    if not reactions.exists():
        return JsonResponse({
            'status': 'success',
            'data': {'nodes': {}, 'edges': [], 'stats': {'node_count': 0, 'edge_count': 0}}
        })

    # ─────────────────────────────────────────────────────────────────────────
    # 3. Collect all SavedMetabolite IDs to batch-fetch them
    # ─────────────────────────────────────────────────────────────────────────
    saved_met_ids = set()
    
    for reaction in reactions:
        substrates = _safe_json_loads(reaction.substrates, [])
        products = _safe_json_loads(reaction.products, [])
        subs_types = _safe_json_loads(reaction.substrates_types, [])
        prods_types = _safe_json_loads(reaction.products_types, [])
        
        for met_id, met_type in zip(substrates, subs_types):
            if met_type == 'Saved':
                try:
                    saved_met_ids.add(int(met_id))
                except (ValueError, TypeError):
                    pass
        
        for met_id, met_type in zip(products, prods_types):
            if met_type == 'Saved':
                try:
                    saved_met_ids.add(int(met_id))
                except (ValueError, TypeError):
                    pass

    # Batch fetch all SavedMetabolites
    saved_mets_lookup = {}
    if saved_met_ids:
        saved_mets = SavedMetabolite.objects.filter(id__in=saved_met_ids)
        saved_mets_lookup = {sm.id: sm for sm in saved_mets}

    # ─────────────────────────────────────────────────────────────────────────
    # 4. Build nodes and edges
    # ─────────────────────────────────────────────────────────────────────────
    nodes = {}
    edges = []

    for reaction in reactions:
        # Parse reaction data
        substrates = _safe_json_loads(reaction.substrates, [])
        products = _safe_json_loads(reaction.products, [])
        subs_types = _safe_json_loads(reaction.substrates_types, [])
        prods_types = _safe_json_loads(reaction.products_types, [])
        subs_names = _safe_json_loads(reaction.substrates_names, [])
        prods_names = _safe_json_loads(reaction.products_names, [])
        subs_comps = _safe_json_loads(reaction.subs_comps, [])
        prods_comps = _safe_json_loads(reaction.prods_comps, [])
        subs_found = _safe_json_loads(reaction.subs_found, [])
        prods_found = _safe_json_loads(reaction.prod_found, [])
        
        # Combined formulas (substrates + products)
        all_formulas = _safe_json_loads(reaction.metabolite_formulas, [])
        n_subs = len(substrates)

        substrate_node_ids = []
        product_node_ids = []

        # Process substrates
        for i, (met_id, met_type) in enumerate(zip(substrates, subs_types)):
            node_id = _get_node_id(met_id, met_type)
            substrate_node_ids.append(node_id)
            
            # Only add node if not already present
            if node_id not in nodes:
                name = subs_names[i] if i < len(subs_names) else str(met_id)
                comp = subs_comps[i] if i < len(subs_comps) else None
                in_vmh = subs_found[i] if i < len(subs_found) else False
                formula = all_formulas[i] if i < len(all_formulas) else None
                
                if met_type.lower() == 'vmh':
                    nodes[node_id] = _build_vmh_node(met_id, name, comp, in_vmh, formula)
                elif met_type == 'Saved':
                    saved_met = saved_mets_lookup.get(int(met_id))
                    if saved_met:
                        nodes[node_id] = _build_saved_node(saved_met, comp)
                    else:
                        # Fallback if SavedMetabolite not found
                        nodes[node_id] = {
                            'id': node_id,
                            'type': 'saved',
                            'name': name,
                            'compartment': comp,
                            'in_vmh': False,
                            'formula': formula,
                        }

        # Process products
        for i, (met_id, met_type) in enumerate(zip(products, prods_types)):
            node_id = _get_node_id(met_id, met_type)
            product_node_ids.append(node_id)
            
            if node_id not in nodes:
                name = prods_names[i] if i < len(prods_names) else str(met_id)
                comp = prods_comps[i] if i < len(prods_comps) else None
                in_vmh = prods_found[i] if i < len(prods_found) else False
                formula_idx = n_subs + i
                formula = all_formulas[formula_idx] if formula_idx < len(all_formulas) else None
                
                if met_type.lower() == 'vmh':
                    nodes[node_id] = _build_vmh_node(met_id, name, comp, in_vmh, formula)
                elif met_type == 'Saved':
                    saved_met = saved_mets_lookup.get(int(met_id))
                    if saved_met:
                        nodes[node_id] = _build_saved_node(saved_met, comp)
                    else:
                        nodes[node_id] = {
                            'id': node_id,
                            'type': 'saved',
                            'name': name,
                            'compartment': comp,
                            'in_vmh': False,
                            'formula': formula,
                        }

        # Build edge for this reaction
        edge = _build_edge(reaction, substrate_node_ids, product_node_ids)
        edges.append(edge)

    # ─────────────────────────────────────────────────────────────────────────
    # 5. Compute statistics
    # ─────────────────────────────────────────────────────────────────────────
    vmh_node_count = sum(1 for n in nodes.values() if n.get('type') == 'vmh')
    saved_node_count = sum(1 for n in nodes.values() if n.get('type') == 'saved')
    
    stats = {
        'node_count': len(nodes),
        'edge_count': len(edges),
        'vmh_metabolites': vmh_node_count,
        'saved_metabolites': saved_node_count,
        'reversible_reactions': sum(1 for e in edges if e.get('reversible')),
    }

    # ─────────────────────────────────────────────────────────────────────────
    # 6. Return response
    # ─────────────────────────────────────────────────────────────────────────
    return JsonResponse({
        'status': 'success',
        'data': {
            'nodes': nodes,
            'edges': edges,
            'stats': stats,
        }
    })