def _list_value(values, index, default=''):
    if index < len(values):
        return values[index]
    return default


def _clean_text(value):
    if value is None:
        return ''
    return str(value).strip()


def split_reaction_names(names_dict):
    substrates_names = []
    products_names = []

    for key, value in (names_dict or {}).items():
        if 'substrate' in key:
            substrates_names.append(value)
        elif 'product' in key:
            products_names.append(value)
        else:
            raise ValueError(f"Invalid key: {key}")

    return substrates_names, products_names


def normalize_reaction_side(
        metabolites,
        metabolite_types,
        stoichiometries,
        compartments,
        names=None,
        uploaded_file_count=0):
    normalized = {
        'metabolites': [],
        'types': [],
        'stoichiometries': [],
        'compartments': [],
        'names': [],
    }
    names = names or []
    max_len = max(
        len(metabolites),
        len(metabolite_types),
        len(stoichiometries),
        len(compartments),
        len(names),
    )
    consumed_files = 0

    for index in range(max_len):
        metabolite = _clean_text(_list_value(metabolites, index))
        metabolite_type = _clean_text(_list_value(metabolite_types, index, 'VMH')) or 'VMH'
        has_uploaded_file = (
            metabolite_type == 'MDL Mol file'
            and consumed_files < uploaded_file_count
        )

        if not metabolite and not has_uploaded_file:
            continue

        if metabolite_type == 'MDL Mol file':
            consumed_files += 1

        normalized['metabolites'].append(metabolite)
        normalized['types'].append(metabolite_type)
        normalized['stoichiometries'].append(
            _clean_text(_list_value(stoichiometries, index, '1')) or '1'
        )
        normalized['compartments'].append(
            _clean_text(_list_value(compartments, index, 'c')) or 'c'
        )
        normalized['names'].append(_clean_text(_list_value(names, index)))

    return normalized


def reaction_has_any_side(substrates, products):
    return bool(substrates) or bool(products)
