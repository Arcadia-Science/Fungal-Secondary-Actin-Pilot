#!/usr/bin/env python3
import requests
import json

def get_ott_id_for_order(order):
    response = requests.post(
        "https://api.opentreeoflife.org/v3/tnrs/match_names",
        headers={"Content-Type": "application/json"},
        json={"names": [order]}
    )
    if response.status_code == 200:
        data = json.loads(response.text)
        try:
            return data['results'][0]['matches'][0]['taxon']['ott_id']
        except (IndexError, KeyError):
            return None
    else:
        return None

def get_representative_species(ott_id):
    payload = {"node_id": f"ott{ott_id}", "num_descendants": 1}
    response = requests.post(
        "https://api.opentreeoflife.org/v3/tree_of_life/subtree",
        headers={"Content-Type": "application/json"},
        json=payload
    )
    if response.status_code == 200:
        data = json.loads(response.text)
        if 'newick' in data:
            try:
                newick_string = data['newick']
                newick_string = newick_string.lstrip('(')
                first_species = newick_string.split(')', 1)[0].split(',', 1)[0]
                return first_species
            except IndexError:
                return None
        else:
            return None
    else:
        return None

def main():
    output_list = []
    with open('fungal_orders.txt', 'r') as f:
        orders = f.read().strip().split('\n')

    for order in orders:
        ott_id = get_ott_id_for_order(order)
        if ott_id:
            rep_species = get_representative_species(ott_id)
            if rep_species:
                output_list.append(f"{order},{rep_species}")

    print('\n'.join(output_list))

if __name__ == "__main__":
    main()
