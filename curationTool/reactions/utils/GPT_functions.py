
# IMPORTS 
import os
import pandas as pd
import numpy as np
import re
import pubchempy as pcp
import coreapi
import re
from openai import OpenAI
import time
import pubchempy as pcp
import requests
from itertools import permutations
from IPython.display import display, HTML


# create a new openai api key
#os.environ["OPENAI_API_KEY"] = "sk-proj-KDaeymwJ58UELawwkqwVT3BlbkFJ3raLweuqMjdC0HqzDHzE"

# set up openai api key
openai_api_key = "sk-proj-KDaeymwJ58UELawwkqwVT3BlbkFJ3raLweuqMjdC0HqzDHzE"#os.environ.get('OPENAI_API_KEY')


def get_gene_name(gene_id):
    client = coreapi.Client()
    
    schema = client.get("https://www.vmh.life/_api/genes/?gene_number=" + str(gene_id))
    return schema['results'][0]['symbol']

def get_ncbi_gene_id(gene_name):
    client = coreapi.Client()
    
    schema = client.get("https://www.vmh.life/_api/genes/?symbol=" + gene_name)
    return schema['results'][0]['gene_number']

def get_vmh_synonyms(abbreviation):
    client = coreapi.Client()
    try:
        schema = client.get("https://www.vmh.life/_api/metabolites/?abbreviation=" + abbreviation + "&format=json")
        return schema['results'][0]['synonyms'].split('***')
    except: 
        schema = abbreviation
        return [schema]
    
def get_vmh_met_from_inchi(inchi, met):
    client = coreapi.Client()
    try:
        schema = client.get("https://www.vmh.life/_api/metabolites/?inchiString=" + inchi + "&format=json")
        return schema['results'][0]['abbreviation']
    except: 
        
        return met
    
    
def get_gene_reactions(ncbi_id):
    
    reactions = []
    # Initialize a client & load the schema document
    client = coreapi.Client()
    schema = client.get("https://www.vmh.life/_api/genereactions/" + ncbi_id)

    for item in schema['results']:
        reactions.append(item['formula'])

    return reactions

def extract_compounds(chemical_formula):
    # Split the formula by spaces and then filter out the '+' and '->' symbols
    elements = [x.strip() for x in chemical_formula.split() if x not in ['+', '->', '<=>']]
    
    # Extract compounds by removing the part inside brackets, if present
    # and removing any leading numbers
    compounds = [re.sub(r'^\.?\s*', '', elem.split('[')[0]) for elem in elements]
    
    # Filter out empty strings
    compounds = [comp for comp in compounds if comp]
    
    return compounds

def vmh_to_normal(formula):
    new_formula = formula
    client = coreapi.Client()
    
    compounds = extract_compounds(formula)
    
    for compound in compounds: 
        
        try:
            schema = client.get("https://www.vmh.life/_api/metabolites/?abbreviation=" +compound)
            normal_name = schema["results"][0]['fullName']    
#             print(compound, normal_name)

            new_formula = new_formula.replace(compound + "[", normal_name + "[")
        except:
            new_formula = new_formula.replace(compound, compound)

        
    
    return new_formula

def get_cid_vmh_api(name):
    cid = 0
    try:
        client = coreapi.Client()
        
        # getting pubchem ID from VMH API
        schema = client.get("https://www.vmh.life/_api/metabolites/?abbreviation=" + name)
        cid = schema["results"][0]['pubChemId']
        
        
    except:
        pass
    return 0 if cid == '' else int(cid)


def get_cid_api(name):
    
#     c = pcp.get_compounds(name, 'name')[0]
#     cid = c.to_dict()['cid']
#     return cid
    try:
        c = pcp.get_compounds(name, 'name')[0]
        cid = c.to_dict()['cid']
    except:
        cid = name
    return cid


def get_metanetx_id(vmh_name):
    api_url = "https://beta.metanetx.org/cgi-bin/mnxweb/id-mapper"
    params = {
        'query_index': 'chem',
        'output_format': 'JSON',
        'query_list': 'vmhM:' + vmh_name
    }
    
    response = requests.get(api_url, params=params)
    
    if response.status_code == 200:
        
        data = response.json()
        # Extracting mnx_id from the response
        mnx_id = data.get('vmhM:' + vmh_name, {}).get('mnx_id')
        return mnx_id
    else:
        return "Error: " + str(response.status_code)
    
    
def get_metanetx_id_by_name(name):
    api_url = "https://beta.metanetx.org/cgi-bin/mnxweb/search"
    params = {
        'format': 'json',
        'db': 'chem',
        'query': name
    }
    
    response = requests.get(api_url, params=params)
    
    if response.status_code == 200:
        
        data = response.json()
        # Extracting mnx_id from the response
#         mnx_id = data.get('vmhM:' + vmh_name, {}).get('mnx_id')
        return data[0]['mnx_id']
    else:
        return "Error: " + str(response.status_code)
    
def parse_metabolic_reactions(text):
    
    # Remove numbers
    equation_no_numbers = re.sub(r'\b\d+\.?\d*\b\s*', '', text)
    
    # Remove brackets and their contents
    equation_no_brackets = re.sub(r'\[\w+\]', '', equation_no_numbers)
    
    # Split the cleaned text into separate reactions based on newlines
    reaction_lines = equation_no_brackets.split('\n')
    
    parsed_reactions = []
    for reaction_line in reaction_lines:
        # Split individual reactions into reactants and products using updated regular expression
        reactions = re.split(r' <=> | = | <--> | → | -> | ⇌ | ↔ | <-> ', reaction_line)
        
        # Remove any leading or trailing spaces from each reaction
        reactions = [reaction.strip() for reaction in reactions if reaction.strip()]
        
        # Combine reactants and products into a single list and break down into its constituent metabolites
        combined_reaction = ' + '.join(reactions)
        metabolites = combined_reaction.split('+')
        metabolites = [metabolite.strip() for metabolite in metabolites if metabolite.strip()]
        
        parsed_reactions.append(metabolites)
    
    return parsed_reactions


def parse_metabolic_reactions_gpt(text):
    # Remove reaction numbers and replace multiple new lines with a single one
    clean_text = re.sub(r'\d+\.', '', text)
    clean_text = re.sub(r'\n+', '\n', clean_text).strip()
    
    # Split the cleaned text into separate reactions based on newlines
    reaction_lines = clean_text.split('\n')
    
    parsed_reactions = []
    for reaction_line in reaction_lines:
        # Split individual reactions into reactants and products using updated regular expression
        reactions = re.split(r' <=> | => | = | --> | → | -> | ⇌ | ↔ ', reaction_line)
        
        # Remove any leading or trailing spaces from each reaction
        reactions = [reaction.strip() for reaction in reactions if reaction.strip()]
        if len(reactions)>1:    
            part1 = reactions[0].split('+')
            part1 = [metabolite.strip() for metabolite in part1 if metabolite.strip()]
            part2 = reactions[1].split('+')
            part2 = [metabolite.strip() for metabolite in part2 if metabolite.strip()]

            metabolites = part1 + ['->'] + part2

            # Using list comprehension to replace '.' and ':' in each string
            metabolites = [s.replace('.', '').replace(':', '') for s in metabolites]       
            parsed_reactions.append(metabolites)
    
    return parsed_reactions


def help_format_answer_with_gpt(raw_answer):
    #print("raw answer: ", raw_answer)
    question = "Given the answer: " + raw_answer + ". Only give the chemical reaction(s) and nothing else i.e. no text. Start with '1.' indicating the first reaction. Only text not latex and use normal letters e.g. H2O instead of H₂O. Also remove electron charge e.g. NAD instead of NAD⁺ "
    message=[{"role": "user", "content": question}]
    client = OpenAI()
    response = client.chat.completions.create(
    model="gpt-4o",
    max_tokens=150,
    temperature=0,
    messages = message)
    
    answer = response.choices[0].message.content 
    return answer

def askGPT4(geneName, temperature):
    import json
    # what is the protein associated with gene X
    # metabolic reactions catalyzed by the protein Y
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    config_path = os.path.join(base_dir, 'config.json')
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)
        os.environ["OPENAI_API_KEY"] = config['API_KEY']
    question = "what are all the metabolic reaction(s) catalyzed by the gene: " + geneName + "? only give the chemical reactions and no other information. Note there may be more than 1 reaction."
    message=[{"role": "user", "content": question}]
    client = OpenAI()
    response = client.chat.completions.create(
        model="gpt-4o",
        max_tokens=200,
        temperature=temperature,  # Use the provided temperature
        messages=message
    )
    answer = response.choices[0].message.content
    return help_format_answer_with_gpt(answer)


def metanetx_to_inchi(metabolite):
    x = metabolite
    api_url = "https://beta.metanetx.org/sparql/?default-graph-uri=https%3A%2F%2Frdf.metanetx.org.beta%2F&query=PREFIX+mnx%3A+%3Chttps%3A%2F%2Frdf.metanetx.org%2Fschema%2F%3E%0D%0APREFIX+owl%3A+%3Chttp%3A%2F%2Fwww.w3.org%2F2002%2F07%2Fowl%23%3E%0D%0APREFIX+rdf%3A+%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23%3E%0D%0APREFIX+rdfs%3A+%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23%3E%0D%0APREFIX+chebi%3A+%3Chttp%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FCHEBI_%3E%0D%0ASELECT+%3Finchi%0D%0AWHERE+%7B%0D%0A++++%3Fmetabolite+a+mnx%3ACHEM+.%0D%0A++++%3Fmetabolite+rdfs%3Alabel+%3Flabel+.%0D%0A++++%3Fmetabolite+rdfs%3Acomment+%27" + x + "%27+.%0D%0A++++%3Fmetabolite+mnx%3AchemRefer+%3Freference%0D%0A++++%0D%0A++++OPTIONAL+%7B+%3Fmetabolite+mnx%3Ainchi++++%3Finchi+%7D%0D%0A++++%0D%0A%7D%0D%0A%0D%0A%0D%0A&format=application%2Fsparql-results%2Bjson"
    params = {
        'format': 'json',
        'db': 'chem',       
        
    }
    
    response = requests.get(api_url)#, params=params)
    
    if response.status_code == 200:
        
        data = response.json()
        # Extracting mnx_id from the response
#         mnx_id = data.get('vmhM:' + vmh_name, {}).get('mnx_id')
        #return data[0]
        try:
            return data['results']['bindings'][0]['inchi']['value']
        except:
            return metabolite
    else:
        return "Error: " + str(response.status_code)

def pubchem_similarity_search(identifier_type, identifier, format="JSON", **kwargs):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    full_url = f"{base_url}/fastsimilarity_2d/{identifier_type}/{identifier}/cids/{format}"
    if kwargs:
        params = "?"
        for key, value in kwargs.items():
            params += f"{key}={value}&"
        full_url += params[:-1]
    
    response = requests.get(full_url)
    response.raise_for_status()
    
    if format == "JSON":
        return response.json()
    else:
        return response.text
    
    
def jaccard_similarity(set1, true_len_set2, set2):
    
    intersection = len(set1.intersection(set2))
    
    union = true_len_set2 + len(set1) - intersection
    return intersection / union if union != 0 else 0

def flatten_extend(matrix):
    flat_list = []
    for row in matrix:
        flat_list.extend(row)
    return flat_list

def evaluate_predictions(predicted_reactions, ground_truth_names, ground_truth_reactions, tracking=False):
    max_jaccard_scores = []
    tracker = []
    
    # Step 1: Jaccard Similarity
    
    for ix, pred in enumerate(predicted_reactions):
        #print(ix)
        tracker_temp = 0
        max_jaccard = 0
        for idx, truth in enumerate(ground_truth_reactions):  
            
            jaccard = jaccard_similarity(set(pred), len(ground_truth_names[idx]), set(flatten_extend(truth)))
          
            if jaccard > max_jaccard:
                tracker_temp = idx
                max_jaccard = jaccard
                
        max_jaccard_scores.append(max_jaccard)
        tracker.append(tracker_temp)
    
    # Step 2 & 3: Over-prediction and Under-prediction
    num_predicted = len(predicted_reactions)
    num_truth = len(ground_truth_reactions)
    over_prediction_penalty = abs(num_predicted - num_truth)
    
    # Step 4: Composite Score
    avg_jaccard = sum(max_jaccard_scores) / num_predicted if num_predicted > 0 else 0
    
    if tracking:
        return {"Average Jaccard": avg_jaccard, "Over-prediction Penalty": over_prediction_penalty}, avg_jaccard, tracker

    else:
        return {"Average Jaccard": avg_jaccard, "Over-prediction Penalty": over_prediction_penalty}, avg_jaccard
    
def check_match(ix, idx, mid, tracker):
   
    vmh_reaction = ground_truths[ix][tracker[idx]]  
    
    for metabolite in vmh_reaction:
        if mid in metabolite:
            return True
    
    return False

def print_output(entrez_id, ix):
    print("Entrez ID: ", entrez_id)
    gene_name = get_gene_name(entrez_id)
    print("Gene name: ", gene_name)
    print_statement, score, tracker = evaluate_predictions(gpt_predictions[ix], ground_truths_names[ix], ground_truths[ix], True)
    print("Jaccard score: ", score, "\n")
    
    #------------------------------------------------------------------------------------
    print("Ground truth reactions:")
    for idx, item in enumerate(ground_truths_names[ix]):
        print(item)
        html = ''
        for i, cid in enumerate(ground_truths[ix][idx]):
            
            print(item[i])
            html += '<a href="https://pubchem.ncbi.nlm.nih.gov/compound/' + str(cid) + '" >' + str(cid) + '</a>' + " "
            
            display(HTML(html))
            print("\n")
        
        print("\n")
        
    # ---------------------------------------------------------
    
    print("Raw GPT output:")
    print(gpt_raw_output[ix], "\n")
    
    #------------------------------------------------------
    
    print("GPT CIDS: \n")
    
    for idx, item in enumerate(gpt_predictions_mid_and_metabolite[ix]):
        
        print("[REACTION ", idx+1, "] :")
        
        for y in item:
            x = y.copy()
            html = x[0] + ' : '
            mid = x[1]
            match = check_match(ix, idx, mid, tracker)
            
            if match:
                html+= '<a style="color: green" id="'+ mid + gene + str(idx) + ' " href="https://www.metanetx.org/chem_info/' + mid + '" >' + mid + '</a>' + ', '

            else:
                html+= '<a style="color: red" id="'+ mid + gene + str(idx) + ' " href="https://www.metanetx.org/chem_info/' + mid + '" >' + mid + ' </a>' + ', '

            display(HTML(html))

        print("\n")        
   
    print("\n------------------------------------\n\n\n")
