"""
This module provides various functions related to gene data retrieval, parsing, 
and mapping to metabolic reactions.

It includes:
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
from reactions.organ_data import ORGAN_MAPPING, location_mapping
from reactions.utils.utils import fetch_and_map_gene_expression, get_subcellular_locations

from django.conf import settings

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

        # If not found in Entrez, check VMH
        vmh_base_url = settings.BASE_URL
        vmh_endpoint = f"{vmh_base_url}_api/genes/?gene_number={gene_input}"
        vmh_response = requests.get(vmh_endpoint, verify=False,timeout=10)

        if vmh_response.status_code != 200:
            return JsonResponse(
                {
                    'error': True,
                    'message': f'VMH API returned error {vmh_response.status_code}'
                     f'for gene number `{gene_input}`'},
                status=500)

        vmh_data = vmh_response.json()
        if vmh_data['count'] == 0:
            return JsonResponse(
                {'error': True, 'message': f'Gene number `{gene_input}` not found in VMH'}, status=404) # pylint: disable=line-too-long

        gene = vmh_data['results'][0]
        symbol = gene.get('symbol', '')
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
    Parse a gene logical expression and validate its format.

    Process:
        - Ensures the input follows logical statement rules (e.g., `GENE1 AND GENE2`).
        - Identifies errors in the logical structure.
        - Returns the processed string or an error message.

    Parameters:
        request (HttpRequest): The HTTP request containing a JSON with `geneinfo`.

    Returns:
        JsonResponse:
            - Success: JSON with the processed logical statement.
            - Error: JSON with an error message if the format is incorrect.
    """
    if request.method == 'POST':
        data = json.loads(request.body)
        statement = data.get('geneinfo', '')

        # Patterns for various parts of the logical statement
        alphanumeric_pattern = r'[A-Za-z0-9]+'
        operator_pattern = r'(AND|OR)'

        # Full pattern combining the subpatterns
        full_pattern = fr'^{alphanumeric_pattern}\s*({operator_pattern}\s*{alphanumeric_pattern}\s*)*$' # pylint: disable=line-too-long

        # Compile the regular expression
        pattern = re.compile(full_pattern)

        # Match the statement against the pattern
        match = pattern.match(statement)

        response_data = {
            'processed_string': statement if match else None,
            'error': None
        }

        # If there's no match, determine why
        if not match:
            if not re.match(alphanumeric_pattern, statement):
                response_data['error'] = 'The statement must start with a Gene.'
            elif not re.search(fr'\s*{operator_pattern}\s*', statement):
                response_data['error'] = (
                    "The statement must contain at least one AND/OR operator "
                    "after an alphanumeric string."
                    )
            elif re.search(r'[^A-Za-z0-9\s\(\)ANDOR]', statement):
                response_data['error'] = (
                    "The statement contains invalid characters. Only alphanumeric characters, "
                    "spaces, parentheses, and the words AND/OR are allowed."
                )
            else:
                response_data['error'] = (
                    "The statement does not match the required logical pattern. "
                    "Ensure it follows the structure: alphanumeric (AND/OR alphanumeric)."
                )


        return JsonResponse(response_data)
    return JsonResponse({'error': 'Invalid request method'}, status=405)


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
    # Construct the full path to the config.json file
    config_path = os.path.join(base_dir, 'config.json')
    # Load the config file
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)

    # Extract the file_path from the config
    file_path = config.get('file_path')

    # Construct the full path to the file using base_dir and file_path from
    # the config
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
