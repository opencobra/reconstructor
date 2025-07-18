from collections import OrderedDict
from urllib.parse import quote, unquote
import uuid
import os
import re
import xml.etree.ElementTree as ET
import requests
import json
from reactions.models import SavedMetabolite
from django.core import serializers 


def capitalize_first_letter(s):
    return re.sub(
        r'^[^a-zA-Z]*([a-zA-Z])',
        lambda match: match.group(0).upper(),
        s.lower())


def gen_replace_dict(substrates, products, subs_types, prods_types):
    """
    Generates a replacement dictionary mapping each unique substrate and product to a unique key that RDT will generate.

    Input:
    - substrates (list): List of substrate names or identifiers.
    - products (list): List of product names or identifiers.

    Output:
    - (OrderedDict): An ordered dictionary where each unique substrate and product is mapped to a unique key.
    """
    substrates = [unquote(substrate).split('\n')[0].split(' ')[1] if sub_type == 'Draw' else substrate
                  for substrate, sub_type in zip(substrates, subs_types)]

    products = [unquote(product).split('\n')[0].split(' ')[1] if prod_type ==
                'Draw' else product for product, prod_type in zip(products, prods_types)]

    def generate_key(index):
        return f"M{index:05d}"
    # Create the replacement dictionary
    replacement_dict = OrderedDict()
    index = 1
    # Process substrates list
    for substrate in substrates:
        if substrate not in replacement_dict:
            replacement_dict[substrate] = generate_key(index)
            index += 1
    # Process products list
    for product in products:
        if product not in replacement_dict:
            replacement_dict[product] = generate_key(index)
            index += 1
    replacement_dict = {
        value: key for (
            key,
            value) in replacement_dict.items()}
    return replacement_dict


def get_fields(request, mols_list, mols_types, MEDIA_ROOT, MEDIA_URL, side):
    fields = []
    file_idx = 0
    for i in range(len(mols_list)):
        if mols_types[i] != 'MDL Mol file':
            fields.append(mols_list[i])
        else:
            file_input = request.FILES.getlist(side)[file_idx]
            unique_filename = f"{uuid.uuid4()}"
            # Define the permanent file save path
            save_path = os.path.join(MEDIA_ROOT, 'mol_files', unique_filename)

            with open(save_path, 'wb+') as permanent_file:
                # Write the contents of the uploaded file to the permanent file
                for chunk in file_input.chunks():
                    permanent_file.write(chunk)
            save_path = os.path.join(MEDIA_URL, 'mol_files', unique_filename)
            fields.append(save_path)
            file_idx += 1
    return fields


def seperate_metab_names(names_dict):
    # Initialize lists for substrates and products names
    substrates_names = []
    products_names = []

    # Separate and sort the substrates and products based on their keys
    for key, value in sorted(names_dict.items(), key=lambda x: (x[0], int(
            x[0][len('substrate' if 'substrate' in x[0] else 'product'):]))):
        if key.startswith('substrate'):
            substrates_names.append(value)
        elif key.startswith('product'):
            products_names.append(value)
    return substrates_names, products_names


def clean_dict_keys(input_dict):
    cleaned_dict = {}
    for key, value in input_dict.items():
        if key.endswith('edited'):
            new_key = key.replace('edited', '')
        else:
            new_key = key
        cleaned_dict[new_key] = value
    return cleaned_dict


def check_edited_keys(names_dict, key_type):
    """
    Check if keys in the dictionary contain 'edited' for the specified key type.

    Parameters:
    names_dict (dict): Dictionary with keys to check.
    key_type (str): The type of keys to check ('substrate' or 'product').

    Returns:
    list: List of boolean values indicating presence of 'edited' in the keys.
    """
    result = []
    for key in names_dict:
        if key.startswith(key_type):
            result.append('edited' in key)
    return result


def parse_xml(xml_data):
    # Parse the XML data
    root = ET.fromstring(xml_data)

    # Initialize a dictionary to hold the extracted information
    pubmed_info = {
        'authors': [],
        'title': None,
        'abstract': None
    }

    # Extract the title
    title_element = root.find('.//ArticleTitle')
    if title_element is not None:
        pubmed_info['title'] = title_element.text

    # Extract the abstract
    abstract_element = root.find('.//Abstract/AbstractText')
    if abstract_element is not None:
        pubmed_info['abstract'] = abstract_element.text

    # Extract authors
    for author in root.findall('.//AuthorList/Author'):
        lastname = author.find('LastName')
        forename = author.find('ForeName')
        initials = author.find('Initials')
        # Construct a full name, considering which elements are present
        full_name = ' '.join(filter(None, [
            forename.text if forename is not None else '',
            lastname.text if lastname is not None else '',
            f"({initials.text})" if initials is not None else ''
        ]))
        if full_name:  # If a full name was constructed, add it to the list
            pubmed_info['authors'].append(full_name)
    return pubmed_info


def fetch_and_map_gene_expression(gene_name, df, mapping):
    filtered_df = df[(df['Gene name'] == gene_name) &
                     (df['Level'].isin(['Medium', 'High']))]
    unique_tissues = filtered_df['Tissue'].unique()
    mapped_tissues = {mapping[tissue]
                      for tissue in unique_tissues if tissue in mapping and mapping[tissue] != "/"}
    formatted_result = ', '.join(
        f'{tissue}_' for tissue in sorted(mapped_tissues))
    return formatted_result, None


def get_subcellular_locations(gene_name):
    """
    Function to get subcellular locations from UniProt for a given gene name.
    """
    base_url = 'https://rest.uniprot.org/uniprotkb/search'
    params = {
        'query': f'gene_exact:{gene_name}',
        'format': 'json'
    }

    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        data = response.json()
        subcellular_locations = set()
        for entry in data.get('results', []):
            for comment in entry.get('comments', []):
                if comment.get('commentType') == 'SUBCELLULAR LOCATION':
                    for subcell_location in comment.get(
                            'subcellularLocations', []):
                        location_value = subcell_location.get(
                            'location', {}).get('value')
                        if location_value:
                            subcellular_locations.add(location_value)
        return list(subcellular_locations)
    else:
        print(f"Error: {response.status_code}")
        return None


def safe_json_loads(data):
    return json.loads(data) if data is not None else None

def get_external_ids(mols,types):
    external_ids_list = []
    all_external_fields = ['keggId', 'pubChemId', 'cheBlId', 'hmdb', 'metanetx']
    for mol, met_type in zip(mols, types):
        if met_type == 'Saved':
            saved_metab = SavedMetabolite.objects.get(pk=mol)
            external_ids = {
                'keggId': saved_metab.keggId or '',
                'pubChemId': saved_metab.pubChemId or '',
                'cheBlId': saved_metab.cheBlId or '',
                'hmdb': saved_metab.hmdb or '',
                'metanetx': saved_metab.metanetx or ''
            }
        else:
            external_ids = {field: '' for field in all_external_fields}
        external_ids_list.append(external_ids)
    return external_ids_list

def get_mol_weights(mols, types):
    mol_weights = []
    for mol, met_type in zip(mols, types):
        if met_type == 'Saved':
            saved_metab = SavedMetabolite.objects.get(pk=mol)
            mol_w = saved_metab.mol_w if saved_metab.mol_w else None
            mol_weights.append(mol_w)
        else:
            mol_weights.append(None)
    return mol_weights

def parse_mol_formula(formula):
    """
    Parses a molecular formula string (e.g., "C6H12O6") and returns a dictionary with element counts.
    """
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    counts = {}
    for elem, count in matches:
        counts[elem] = counts.get(elem, 0) + (int(count) if count else 1)
    return counts

def reactions_to_json(qs):
    """
    Serialises a queryset of SavedReaction objects (and any
    related fields you need) into the format consumed by
    handleAdd2VMH.js.
    """
    json_str = serializers.serialize(
        'json',
        qs,
        use_natural_foreign_keys=True,
        # you can control which FK/ManyToMany to follow:
        # follow = ('flags', …)
    )
    return json_str


