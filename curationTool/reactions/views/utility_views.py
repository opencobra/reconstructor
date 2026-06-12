"""
This module contains utility views for handling session data, 
parsing reaction formulas, and fetching external data such as 
PubMed and DOI information.
"""

import json
import re
import requests

from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt

from reactions.models import Reaction
from reactions.utils.utils import parse_xml


def check_session_data(request):
    """
    Delete a specific gene information entry from the session.

    Expects:
        JSON body with 'info_to_delete' key.

    Returns:
        JsonResponse: Success message if deleted, or an error message if failed.
    """
    # Check if 'gene_info' is in session
    gene_info = request.session.get('gene_info', None)

    if gene_info:
        return JsonResponse({'status': 'success', 'gene_info': gene_info})

    return JsonResponse(
        {'status': 'error', 'message': 'No gene info in session.'})

@csrf_exempt
def clear_session(request):
    """
    Clear all session data.

    Returns:
        JsonResponse: Success message if cleared, or an error message if failed.
    """

    if request.method == 'POST':
        try:
            # Clear the entire session
            request.session.flush()  # This will clear all session data

            return JsonResponse({'status': 'success',
                                 'message': 'Session cleared successfully.'})
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=405)


@csrf_exempt
def clear_gene_info_session(request):
    """
    Clear ONLY the pending gene-info bucket from the session, leaving the rest
    of the session (e.g. the logged-in userID) intact.

    The gene panel stages gene info in `session['gene_info']` while creating a
    brand-new reaction, and flushes it into the reaction once one exists. That
    bucket used to be sticky (its JS clear was a no-op), so it leaked into every
    reaction subsequently opened. Switching the open reaction in the tabbed
    workspace (001-multi-reaction-tabbed) calls this first so no stale gene info
    bleeds across reactions (FR-003 / US4 isolation).
    """
    if request.method == 'POST':
        request.session.pop('gene_info', None)
        request.session.modified = True
        return JsonResponse({'status': 'success'})
    return JsonResponse(
        {'status': 'error', 'message': 'Invalid request method.'}, status=405)
@csrf_exempt
def delete_gene_info_from_session(request):
    """
    Delete a specific gene information entry from the session.

    Expects:
        JSON body with the key 'info_to_delete', specifying the gene info to remove.

    Process:
        - Retrieves the existing 'gene_info' list from the session.
        - Filters out the specified entry.
        - Updates the session with the modified list.

    Returns:
        JsonResponse:
            - Success message if the gene info is deleted.
            - Error message if the request method is invalid or an exception occurs.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            info_to_delete = data.get('info_to_delete')

            # Get the gene info list from the session
            gene_info = request.session.get('gene_info', [])

            # Filter out the item to delete
            updated_gene_info = [
                info for info in gene_info if info['info'] != info_to_delete]

            # Update the session with the filtered list
            request.session['gene_info'] = updated_gene_info
            request.session.modified = True

            return JsonResponse({'status': 'success',
                                 'message': 'Gene info deleted from session.'})
        except Exception as e:
            return JsonResponse(
                {'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method.'}, status=405)

@csrf_exempt
def parse_formula_with_compartments(request):
    """
    Parse and format chemical reaction formulas with specified compartments.

    Expects:
        JSON body with 'formulas', 'subs_comps', and 'prods_comps'.

    Returns:
        JsonResponse: Formatted reaction formulas with compartment details.
    """
    if request.method == 'POST':
        try:
            body = json.loads(request.body)
        except json.JSONDecodeError as e:
            return JsonResponse(
                {'error': f'JSON decode error: {str(e)}'}, status=400)

        formulas = body.get('formulas', [])
        subs_comps = body.get('subs_comps', [])
        prods_comps = body.get('prods_comps', [])

        # Log received data for debugging

        detailed_formulas = []

        def extract_components(comp_str):
            """Extracts quantity and name from a component string."""
            match = re.match(r'(\d*)\s*(\w+)', comp_str.strip())
            if match:
                quantity = match.group(1) if match.group(1) else '1'
                name = match.group(2)
                return quantity, name
            return None, None

        for formula in formulas:
            # Normalize formula by ensuring spaces around the arrows
            formula = formula.replace('->', ' -> ').replace('<=>', ' <=> ')
            formula = formula.replace('-> ', ' -> ').replace('<=> ', ' <=> ')
            formula = formula.replace(' ->', ' -> ').replace(' <=>', ' <=> ')

            if ' -> ' in formula:
                direction = 'forward'
                substrates, products = formula.split(' -> ')
            elif ' <=> ' in formula:
                direction = 'reversible'
                substrates, products = formula.split(' <=> ')
            else:
                return JsonResponse(
                    {'error': f'Invalid formula format: {formula}'}, status=400)

            subs_list = [comp.strip() for comp in substrates.split(
                ' + ')] if substrates else []
            prods_list = [comp.strip()
                          for comp in products.split(' + ')] if products else []

            # Remove empty components
            subs_list = [comp for comp in subs_list if comp]
            prods_list = [comp for comp in prods_list if comp]

            # Check if the lengths of subs_comps and prods_comps match the
            # lengths of subs_list and prods_list
            if len(subs_list) != len(subs_comps):
                return JsonResponse(
                    {
                        'error': (
                            "Number of substrate compartments does not "
                            "match number of substrates."
                        )
                    },
                    status=400
                )
            if len(prods_list) != len(prods_comps):
                return JsonResponse(
                    {
                        'error': (
                            "Number of product compartments does not "
                            "match number of products."
                        )
                    },
                    status=400
                )
            subs_detailed = []
            for j, comp in enumerate(subs_list):
                quantity, name = extract_components(comp)
                subs_detailed.append(f"{quantity} {name}[{subs_comps[j]}]")

            prods_detailed = []
            for j, comp in enumerate(prods_list):
                quantity, name = extract_components(comp)
                prods_detailed.append(f"{quantity} {name}[{prods_comps[j]}]")

            if direction == 'forward':
                detailed_formula = ' + '.join(subs_detailed) + \
                    ' -> ' + ' + '.join(prods_detailed)
            else:  # reversible
                detailed_formula = ' + '.join(subs_detailed) + \
                    ' <=> ' + ' + '.join(prods_detailed)

            detailed_formulas.append(detailed_formula)

        return JsonResponse({'detailed_formulas': detailed_formulas})

    return JsonResponse({'error': 'Invalid request method'}, status=400)


@csrf_exempt  # Use csrf_exempt if CSRF token isn't being managed
def save_formula(request):
    """
    Save a given reaction formula to the database.

    Expects:
        JSON body with 'formula' and 'reaction_id'.

    Returns:
        JsonResponse: Success message if saved, or an error if reaction not found.
    """
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            formula = data.get('formula', '')
            reaction_id = data.get('reaction_id', '')

            # Ensure the formula is saved with quotes
            quoted_formula = f'"{formula}"'

            # Find the reaction by the provided reaction ID
            try:
                reaction = Reaction.objects.get(id=reaction_id)
            except Reaction.DoesNotExist:
                return JsonResponse(
                    {'status': 'failed', 'error': 'Reaction not found'}, status=404
                )
            # Update the formulas field with the provided formula string with
            # quotes
            reaction.rxn_formula = quoted_formula
            reaction.save()

            return JsonResponse({'status': 'success'})
        except Exception as e:
            return JsonResponse(
                {'status': 'failed', 'error': str(e)}, status=500)
    return JsonResponse(
        {'status': 'failed', 'error': 'Invalid request method'}, status=400)


def get_pubmed_info(request, pmid):
    """
    Fetch PubMed article details using a given PubMed ID (PMID).

    Parameters:
        pmid (str): The PubMed ID of the article.

    Returns:
        JsonResponse: Contains the article title, authors, and abstract, 
                      or an error message if unsuccessful.
    """
    # Base URL for PubMed API
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'xml'
    }
    response = requests.get(base_url, params=params,timeout=10)

    if response.status_code == 200:
        # Parsing XML response is necessary here to extract needed information
        # This is a placeholder to indicate where parsing should occur
        # You might use libraries like xml.etree.ElementTree or lxml to parse
        # the XML

        pubmed_info = parse_xml(response.content)  # Implement this function
        if pubmed_info['title'] is None:
            return JsonResponse(
                {'status': 'error', 'message': 'PubMed API returned an empty response'}, status=500)
        return JsonResponse({
            'status': 'success',
            'author': pubmed_info['authors'],
            'title': pubmed_info['title'],
            'abstract': pubmed_info['abstract']
        })

    return JsonResponse({'status': 'error',
                            'message': f"PubMed API returned error {response.status_code}"},
                        status=500)


def get_doi_info(request, doi):
    """
    Retrieve metadata of a publication using a DOI from CrossRef API.

    Parameters:
        doi (str): The DOI of the publication.

    Returns:
        JsonResponse: Contains the title, authors, and abstract, 
                      or an error message if unsuccessful.
    """
    doi = doi.replace('DOI:', '')
    base_url = f'https://api.crossref.org/works/{doi}'
    try:
        response = requests.get(base_url, timeout=10)
        if response.status_code == 200:
            data = response.json()['message']

            # Extract needed information
            authors = [
                f"{author['given']} {author['family']}" for author in data.get(
                    'author', [])]
            title = data.get('title', [None])[0]
            abstract = data.get('abstract', None)
            return JsonResponse({
                'status': 'success',
                'author': authors,
                'title': title,
                'abstract': abstract
            })

        return JsonResponse({'status': 'error',
                                'message': 'CrossRef API returned an error'},
                            status=response.status_code)
    except Exception as e:
        return JsonResponse(
            {'status': 'error', 'message': f'Failed to fetch DOI data - {str(e)}'}, status=500)


def chemdoodle_sketcher(request):
    """
    Render the ChemDoodle Sketcher page.

    Returns:
        HttpResponse: Renders the 'chemdoodle_sketcher.html' template.
    """
    return render(request, 'reactions/chemdoodle_sketcher.html')
