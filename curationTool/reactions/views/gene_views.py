"""
This module provides various functions related to gene data retrieval, parsing, 
and mapping to metabolic reactions.

It includes:
- `suggest_genes`: Returns type-ahead gene suggestions (HGNC / Entrez).
- `get_gene_info`: Retrieves gene information from external databases (Entrez, VMH, HGNC).
- `gene_parsing`: Parses gene expressions and validates logical statements.
- `gene_details_view`: Fetches detailed gene-related data including organ and subcellular locations.
- `parse_gene_info`: Extracts organ and subcellular location data from gene annotations.
- `parse_genes`: Extracts individual gene names from logical statements.
- `extract_unique_elements`: Processes lists to extract unique elements.
- `map_locations_to_wbm`: Maps subcellular locations to the WBM categories.
"""

import json
import os
import re
import requests
import pandas as pd
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_GET
from reactions.organ_data import ORGAN_MAPPING, location_mapping
from reactions.utils.utils import fetch_and_map_gene_expression, get_subcellular_locations, _vmh_gene_exists, _entrez_exists, _hgnc_symbol_exists_exact
from reactions.utils.vmh_api import gene_data_new, gene_rows_old, gene_symbol_from_row

from reactions.models import Gene
from django.db.models import Q

@require_GET
def suggest_genes(request):
    """
    Provide type-ahead gene suggestions from the **local Gene database** (instead of querying HGNC API).  
    This allows near real-time response with no external latency.

    Process:
        - Read query token `q` from the query string; ignore if shorter than 2 chars unless it is all digits.
        - Search local Gene table for matches:
            * HGNC symbol (case-insensitive, startswith)
            * Aliases (substring match, case-insensitive)
            * Full gene name (substring match, case-insensitive)
            * Entrez ID (if query is numeric)
        - Return up to 8 results ordered alphabetically by symbol.
        - Each result is always marked `"present": True` since it exists in our DB.

    Parameters:
        request (HttpRequest): Django GET request with query parameter `q` (partial gene token).

    Returns:
        JsonResponse: An object with key `"items"` containing up to 8 dicts:
            {
                "symbol": str,        # HGNC-approved symbol
                "name": str,          # Full gene name
                "hgnc_id": str,       # Stable HGNC identifier
                "entrez_id": str,     # Entrez Gene ID as string (may be "")
                "present": True       # True if present in VMH; False if new
            }
        If no match is found, returns {"items": []}.
    """
    q = (request.GET.get("q") or "").strip()
    if len(q) < 2 and not q.isdigit():
        return JsonResponse({"items": []})

    # Search locally (symbol, aliases, or name)
    results = (
        Gene.objects.filter(
            Q(symbol__istartswith=q) |
            (Q(entrez_id__startswith=q) if q.isdigit() else Q())
        )
        .order_by("symbol")[:8]
    )

    items = []
    for g in results:
        present = _vmh_gene_exists(symbol=g.symbol, entrez_id=g.entrez_id or "")
        items.append({
            "symbol": g.symbol,
            "name": g.name,
            "hgnc_id": g.hgnc_id,
            "entrez_id": g.entrez_id or "",
            "present": bool(present),
        })
    
    return JsonResponse({"items": items})

def get_gene_info(request):
    """
    Retrieve gene information based on the user's input.

    Process:
        - Fetches gene information from Entrez, VMH, or HGNC based on the type of identifier.
        - Returns the gene symbol if found.
        - Handles cases where the gene is not found in the respective databases.

    Parameters:
        request (HttpRequest): The HTTP request object containing 'gene' and 'type'.

    Returns:
        JsonResponse:
            - Success: JSON containing `symbol` or `hgnc_id`.
            - Error: If the gene is not found or an invalid type is provided.
    """
    gene_input = request.POST.get('gene')
    type_input = request.POST.get('type')

    if gene_input.strip() == '':
        return JsonResponse({'error': True, 'message': 'No input provided'})

    if type_input == 'Entrez ID':
        # Check Entrez ID
        entrez_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        entrez_params = {
            "db": "gene",
            "id": gene_input,
            "retmode": "json"
        }

        entrez_response = requests.get(entrez_base_url, params=entrez_params,timeout=10)

        if entrez_response.status_code == 200:
            entrez_data = entrez_response.json()
            if "result" in entrez_data and gene_input in entrez_data["result"]:
                gene_data = entrez_data["result"][gene_input]
                if "name" in gene_data:
                    return JsonResponse(
                        {'error': False, 'symbol': gene_data["name"]})

        # If not found in Entrez, check VMH (new API first, then old fallback).
        vmh_rows = gene_data_new(gene_input)
        if not vmh_rows and "." not in gene_input:
            vmh_rows = gene_data_new(f"{gene_input}.1")

        if not vmh_rows:
            vmh_rows = gene_rows_old(gene_number=gene_input)
        if not vmh_rows and "." not in gene_input:
            vmh_rows = gene_rows_old(gene_number=f"{gene_input}.1")

        if not vmh_rows:
            return JsonResponse(
                {'error': True, 'message': f'Gene number `{gene_input}` not found in VMH'}, status=404) # pylint: disable=line-too-long

        symbol = gene_symbol_from_row(vmh_rows[0]) or vmh_rows[0].get('symbol', '')
        return JsonResponse({'error': False, 'symbol': symbol})

    elif type_input == 'HGNC Symbol':
        # Check HGNC Symbol
        hgnc_base_url = 'https://rest.genenames.org/search/symbol/'
        hgnc_endpoint = f"{hgnc_base_url}{gene_input}"
        hgnc_headers = {'Accept': 'application/json'}
        hgnc_response = requests.get(hgnc_endpoint, headers=hgnc_headers,timeout=10)

        if hgnc_response.status_code != 200:
            return JsonResponse(
                {
                    'error': True,
                    'message': f'HGNC API returned error {hgnc_response.status_code}' 
                    f'for symbol `{gene_input}`'},
                status=500)

        hgnc_data = hgnc_response.json()
        num_found = hgnc_data['response']['numFound']

        if num_found == 0:
            return JsonResponse(
                {'error': True, 'message': f'Gene symbol `{gene_input}` not found in HGNC'}, status=404) # pylint: disable=line-too-long
        if num_found > 1:
            genes = hgnc_data['response']['docs'][:10]
            gene_symbols_and_ids = [
                {'symbol': gene['symbol'], 'hgnc_id': gene['hgnc_id']} for gene in genes]
            message = f"Multiple genes found for symbol `{gene_input}`. Please specify. Found genes: " + ", ".join( # pylint: disable=line-too-long
                [f"{gene['symbol']} (HGNC ID: {gene['hgnc_id']})" for gene in gene_symbols_and_ids])
            return JsonResponse(
                {'error': True, 'message': message}, status=400)
        gene = hgnc_data['response']['docs'][0]
        hgnc_id = gene.get('hgnc_id', '')
        symbol = gene.get('symbol', '')
        return JsonResponse(
            {'error': False, 'hgnc_id': hgnc_id, 'symbol': symbol})
    else:
        return JsonResponse({'error': True, 'message': 'Invalid type input'})


@csrf_exempt
def gene_parsing(request):
    """
    Validate and normalize a user-provided GPR (gene–protein–reaction) expression.

    Process:
        - Ensures the request is POST and parses JSON body for the `geneinfo` string.
        - Rejects illegal characters (anything other than letters, digits, spaces, and parentheses).
        - Tokenizes the string into genes, AND/OR operators, and parentheses.
        - Validates expression grammar (balanced parentheses, correct operator/operand order).
        - Verifies each gene token:
            • If numeric → checks existence as an Entrez Gene ID via NCBI.
            • If alphabetic/alphanumeric → checks exact HGNC symbol existence.
        - On success, pretty-prints spacing around operators and returns the normalized expression.
        - On failure, returns a clear error message and (when applicable) the list of invalid genes.

    Parameters:
        request (HttpRequest): POST request with JSON body containing:
            {
                "geneinfo": "<GPR expression, e.g., 'AOC3 AND (AOC1 OR AOC2)'>"
            }

    Returns:
        JsonResponse:
            - Success (HTTP 200):
                {
                  "processed_string": "<normalized GPR>",
                  "error": None
                }
            - Client error (HTTP 400/405) with details:
                {
                  "processed_string": None,
                  "error": "<reason>",
                  "invalid_genes": ["<gene1>", "<gene2>", ...]  # present only when gene validation fails
                }
    """
    if request.method != 'POST':
        return JsonResponse({'error': 'Invalid request method'}, status=405)

    try:
        data = json.loads(request.body.decode('utf-8'))
    except Exception:
        return JsonResponse({'processed_string': None, 'error': 'Invalid JSON'}, status=400)

    statement = (data.get('geneinfo') or '').strip()
    if not statement:
        return JsonResponse({'processed_string': None, 'error': 'Empty expression'}, status=400)

    # reject illegal characters
    if re.search(r'[^A-Za-z0-9()\s]', statement):
        return JsonResponse({
            'processed_string': None,
            'error': ('The statement contains invalid characters. '
                      'Only letters, numbers, spaces, and parentheses are allowed. '
                      'Use AND/OR as operators.')
        }, status=400)

    # tokenize
    tokens = re.findall(r'AND|OR|\(|\)|[A-Za-z0-9]+', statement, flags=re.IGNORECASE)
    if not tokens:
        return JsonResponse({'processed_string': None, 'error': 'No tokens found'}, status=400)

    # grammar: Expr -> Term ( (AND|OR) Term )*
    stack, expect_operand = [], True
    normalized = []

    for tok in tokens:
        up = tok.upper()
        if up in ('AND', 'OR'):
            if expect_operand:
                return JsonResponse({'processed_string': None,
                                     'error': 'Operator found where a gene or "(" was expected.'}, status=400)
            expect_operand = True
            normalized.append(up)
        elif tok == '(':
            if not expect_operand:
                return JsonResponse({'processed_string': None,
                                     'error': 'Missing operator before "(".'}, status=400)
            stack.append('(')
            expect_operand = True
            normalized.append('(')
        elif tok == ')':
            if expect_operand:
                return JsonResponse({'processed_string': None,
                                     'error': '")" found where a gene was expected.'}, status=400)
            if not stack:
                return JsonResponse({'processed_string': None, 'error': 'Unmatched ")".'}, status=400)
            stack.pop()
            expect_operand = False
            normalized.append(')')
        else:
            # gene token (leave as typed)
            if not expect_operand:
                return JsonResponse({'processed_string': None,
                                     'error': 'Missing operator between genes.'}, status=400)
            expect_operand = False
            normalized.append(tok)

    if stack:
        return JsonResponse({'processed_string': None, 'error': 'Unmatched "(".'}, status=400)
    if expect_operand:
        return JsonResponse({'processed_string': None, 'error': 'Expression ends with an operator.'}, status=400)

    # VALIDATE tokens against HGNC/NCBI WITHOUT altering them
    genes = [t for t in normalized if t not in ('AND', 'OR', '(', ')')]
    invalid = []
    for g in genes:
        ok = _entrez_exists(g) if g.isdigit() else _hgnc_symbol_exists_exact(g)
        if not ok:
            invalid.append(g)

    if invalid:
        unique = sorted(set(invalid), key=str.upper)
        return JsonResponse({
            'processed_string': None,
            'error': 'These genes could not be found in HGNC/NCBI. Please correct them: ' + ', '.join(unique),
            'invalid_genes': unique,
        }, status=400)

    # pretty spacing, keep original gene tokens intact
    pretty = []
    for t in normalized:
        if t in ('(', ')'):
            pretty.append(t)
        elif t in ('AND', 'OR'):
            pretty.append(f' {t} ')
        else:
            pretty.append(t)
    pretty_str = re.sub(r'\s+', ' ', ''.join(pretty)).strip()
    pretty_str = pretty_str.replace('( ', '(').replace(' )', ')')

    return JsonResponse({'processed_string': pretty_str, 'error': None})


@csrf_exempt
def gene_details_view(request):
    """
    Retrieve and process detailed gene-related information.

    Process:
        - Reads a configuration file for data file paths.
        - Loads gene expression data from a CSV file.
        - Extracts and maps gene expression data to organs and subcellular locations.
        - Constructs a structured JSON response with formatted gene details.

    Parameters:
        request (HttpRequest): The HTTP request containing gene-related query data.

    Returns:
        JsonResponse:
            - Success: JSON containing `infoText` with gene-related details.
            - Error: If the data processing fails or an invalid request method is used.
    """
    base_dir = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            '..',
            '..',
            '..'))
    
    # Get file path from environment variable
    file_path = os.getenv('FILE_PATH', 'curationTool/reactions/normal_tissue.tsv')
    
    # Construct the full path to the file
    full_file_path = os.path.join(base_dir, file_path)

    # Load the CSV file using the full file path
    df = pd.read_csv(full_file_path, sep='\t')

    if request.method != 'POST':
        return JsonResponse({"error": "Invalid request method"}, status=405)

    try:
        # Log the raw request body

        data = json.loads(request.body)
        # Extract the nested infoText from the data object
        data_info = data.get("infoText", {})
        data_string = data_info.get("infoText", "")
        if not isinstance(data_string, str):
            raise TypeError("infoText must be a string")

    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)
    except TypeError as e:
        return JsonResponse({"error": str(e)}, status=400)

    genes = parse_genes(data_string)

    gene_details = []
    info_text_parts = [f"GPR: {data_string}"]

    all_organs = set()
    all_subcellular_locations = set()

    for gene in genes:
        if not gene.isalnum():  # Skip if gene is not alphanumeric
            continue

        gene_info = {}
        # Fetch and map gene expression data
        unique_organs, error = fetch_and_map_gene_expression(
            gene, df, ORGAN_MAPPING)
        if error:
            gene_info["ORGAN"] = "Error fetching organs"
        else:
            gene_info["ORGAN"] = unique_organs
            all_organs.add(unique_organs)

        # Get subcellular locations
        subcellular_locations = get_subcellular_locations(gene)
        if subcellular_locations:
            mapped_locations = map_locations_to_wbm(subcellular_locations)
            gene_info["SUBCELLULAR LOCATION"] = mapped_locations
            mapped_locations = extract_unique_elements(mapped_locations)
            all_subcellular_locations.update(mapped_locations)
        else:
            gene_info["SUBCELLULAR LOCATION"] = "Subcellular locations not found"

        gene_details.append(gene_info)

    all_subcellular_locations = list(set(all_subcellular_locations))

    # Combine and deduplicate organ and subcellular location information
    combined_organs_text = ", ".join(all_organs)
    combined_subcellular_text = ", ".join(all_subcellular_locations)

    # Format the final infoText part
    info_text_parts.append(
        f"ORGAN({combined_organs_text}), SUBCELLULAR({combined_subcellular_text})")

    organized_result = {
        "userID": data_info.get("userID", ""),
        "infoType": data_info.get("infoType", ""),
        "extLinkType": data_info.get("extLinkType", ""),
        "refType": data_info.get("refType", ""),
        "reactionId": data_info.get("reactionId", ""),
        "infoText": "; ".join(info_text_parts)
    }

    return JsonResponse(organized_result)


def parse_gene_info(request):
    """
    Parse gene-related information from a formatted string.

    Process:
        - Splits the input into gene sections.
        - Extracts organ and subcellular location data for each gene.
        - Organizes and returns the extracted data in JSON format.

    Parameters:
        request (HttpRequest): The HTTP request containing the `info` parameter.

    Returns:
        JsonResponse:
            - Success: JSON containing gene-organ-subcellular mappings.
            - Error: JSON with an error message if parsing fails.
    """
    info = request.GET.get('info', '')

    if not info:
        return JsonResponse({'error': 'No info provided'}, status=400)

    response_data = {}

    try:
        # Split the input by `; GENE:` to separate each gene's information
        gene_sections = re.split(r';\s*GENE:', info)

        # The first part contains GPR, we'll ignore it as we focus on genes
        for section in gene_sections:
            gene_match = re.match(
                r'([^;]+); ORGAN\(([^)]+)\), SUBCELLULAR\(([^)]+)\)', section)

            if not gene_match:
                continue

            gene = gene_match.group(1).strip()
            organs = [organ.strip()
                      for organ in gene_match.group(2).split(',')]
            subcellular_locations = [loc.strip()
                                     for loc in gene_match.group(3).split(',')]

            if gene not in response_data:
                response_data[gene] = {
                    'Organs': [], 'SubcellularLocations': []}

            response_data[gene]['Organs'].extend(organs)
            response_data[gene]['SubcellularLocations'].extend(
                subcellular_locations)

    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)
    return JsonResponse(response_data)


def parse_genes(data_string):
    """
    Extract individual gene names from a logical expression.

    Process:
        - Splits a string using `AND`, `OR`, or spaces.
        - Filters out non-alphanumeric strings.

    Parameters:
        data_string (str): A logical gene expression string.

    Returns:
        list: A list of extracted gene names.
    """
    # Split by spaces and logical operators, retain only alphanumeric strings
    genes = re.split(r'\s+(?:AND|OR)\s+|\s+', data_string)
    genes = [gene for gene in genes if gene.isalnum()]
    return genes


def extract_unique_elements(input_set):
    """
    Extract unique elements from a set by splitting comma-separated values.

    Process:
        - Iterates through each element and splits it into individual components.
        - Returns a set of unique elements.

    Parameters:
        input_set (set): A set containing comma-separated elements.

    Returns:
        set: A set of unique elements.
    """
    unique_elements = set()
    for item in input_set:
        elements = item.split(',')
        unique_elements.update(elements)
    return unique_elements


def map_locations_to_wbm(subcellular_locations):
    """
    Map UniProt subcellular locations to WBM categories.

    Process:
        - Matches each location with a predefined mapping.
        - Returns the corresponding WBM category if found.

    Parameters:
        subcellular_locations (list): A list of subcellular locations.

    Returns:
        list: A list of mapped WBM categories.
    """
    mapped_locations = []
    for location in subcellular_locations:
        if location in location_mapping:
            mapped_locations.append(location_mapping[location])
    return mapped_locations
