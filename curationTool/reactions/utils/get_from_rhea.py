import requests
import re

def get_from_rhea(reaction_abbreviation):
    BASE_URL = 'https://www.rhea-db.org/rhea/?'
    params = {
        "query": reaction_abbreviation,
        "columns": "rhea-id,equation,chebi-id",
        "format": 'tsv',
        "limit": 10
    }

    response = requests.get(BASE_URL, params=params)

    if response.status_code != 200:
        return {"error": f"Failed to fetch data from API. Status code: {response.status_code}"}

    data = response.text

    if not data:
        return {"error": "No reaction found with this abbreviation"}

    lines = data.splitlines()
    headers = lines[0].split("\t")
    results = [line.split("\t") for line in lines[1:]]

    rhea_data = []
    for result in results:
        rhea_dict = {headers[i]: result[i] for i in range(len(headers))}
        rhea_data.append(rhea_dict)

    if not rhea_data:
        return {"error": "No reaction found with this abbreviation"}

    reaction_info = rhea_data[0]
    equation = reaction_info.get("Equation", "")
    chebi_ids = reaction_info.get("ChEBI identifier", "")

    if not equation or not chebi_ids:
        return {"error": "Equation or ChEBI IDs not found"}

    def parse_equation(equation, chebi_ids):
        chebi_ids_list = chebi_ids.split(';')
        if '=' in equation:
            substrates, products = equation.split('=')
            direction = "forward"
        elif '<=>' in equation:
            substrates, products = equation.split('<=>')
            direction = "reversible"
        elif '=>' in equation:
            substrates, products = equation.split('=>')
            direction = "forward"
        elif '<=' in equation:
            substrates, products = equation.split('<=')
            direction = "reverse"
        else:
            return {"error": "Unknown equation format"}

        substrates = [s.strip() for s in substrates.split('+')]
        products = [p.strip() for p in products.split('+')]

        def parse_components(components):
            parsed_components = []
            stoichiometry_pattern = re.compile(r"^(\d+)\s+(.*)$")
            for component in components:
                match = stoichiometry_pattern.match(component)
                if match:
                    stoichiometry = int(match.group(1))
                    component_name = match.group(2)
                else:
                    stoichiometry = 1
                    component_name = component
                parsed_components.append((stoichiometry, component_name))
            return parsed_components

        parsed_substrates = parse_components(substrates)
        parsed_products = parse_components(products)

        substrate_chebis = chebi_ids_list[:len(parsed_substrates)]
        product_chebis = chebi_ids_list[len(parsed_substrates):]

        return {
            "substrates": parsed_substrates,
            "products": parsed_products,
            "substrate_chebis": substrate_chebis,
            "product_chebis": product_chebis,
            "direction": direction
        }

    parsed_equation = parse_equation(equation, chebi_ids)

    result = {
        "direction": parsed_equation["direction"],
        "subs_sch": [1] * len(parsed_equation["substrate_chebis"]),
        "prod_sch": [1] * len(parsed_equation["product_chebis"]),
        "subs_types": ["ChEBI ID"] * len(parsed_equation["substrate_chebis"]),
        "prods_types": ["ChEBI ID"] * len(parsed_equation["product_chebis"]),
        "substrates": parsed_equation["substrate_chebis"],
        "products": parsed_equation["product_chebis"],
        "subs_comps": ["-"] * len(parsed_equation["substrate_chebis"]),
        "prods_comps": ["-"] * len(parsed_equation["product_chebis"]),
        "subsystem": ""
    }

    return result
