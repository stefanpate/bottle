import requests
import json
import logging
import urllib
from tqdm import tqdm

up_ids_path = '../data/mapping/all_uniprot_ids_from_known_rxns.txt'
jni_uniprot = '../data/mapping/jni_uniprot.json'
base_url = 'https://rest.uniprot.org/uniprotkb/'
batch_size = 500

get_sequence = lambda r: r["sequence"]["value"]
get_organism = lambda x: x["organism"]["scientificName"]
get_existence = lambda x: x["proteinExistence"]
get_reviewed = lambda x: x['entryType']
get_id = lambda x: x["primaryAccession"]

def get_reactions(x):
    rhea_ids = []
    for comment in x["comments"]:
        if comment["commentType"] == 'CATALYTIC ACTIVITY':
            if "physiologicalReactions" in comment['reaction']:
                for elt in comment['reaction']['physiologicalReactions']:
                    rhea_ids.append(elt['reactionCrossReference']['id'])

            elif "reactionCrossReferences" in comment['reaction']:
                for elt in comment['reaction']['reactionCrossReferences']:
                    if elt['database'] == 'Rhea':
                        rhea_ids.append(elt['id'])

    return rhea_ids

def scrape_loop(uniprot_ids, batch_size, res, base_url=base_url, jni_uniprot=jni_uniprot):
    n_batches = len(uniprot_ids) // batch_size
    remainder = len(uniprot_ids) % batch_size
    bad_ids = []
    todo = []

    for i in tqdm(range(n_batches + 1)):
        if i == n_batches:
            id_batch = uniprot_ids[-remainder:]
        else:
            id_batch = uniprot_ids[i * batch_size : (i + 1) * batch_size]

        url = urllib.parse.urljoin(base_url, f"search?query=%28accession%3A{('%29+OR+%28accession%3A').join(id_batch)}%29&size=500")

        response = requests.get(url)
        if response.status_code != 200 or "results" not in response.json():
            logging.warning(f"Error with request for batch {i}")
            message = " ".join(response.json()['messages']).split(' ')
            bad_ids += [elt.strip("'") for elt in message if elt.strip("'") in id_batch]
            todo += [elt for elt in id_batch if elt not in bad_ids]
            continue

        resp = response.json()["results"]
        for r in resp:
            try:
                id = get_id(r)
            except:
                logging.warning(f"Could not find primary accession")
                continue

            res[id] = {'reviewed':None, 'organism':None, 'sequence':None, 'existence':None, 'rhea_ids':[]}

            try:
                res[id]['organism'] = get_organism(r)
            except:
                pass

            try:
                res[id]['reviewed'] = get_reviewed(r)
            except:
                pass

            try:
                res[id]['sequence'] = get_sequence(r)
            except:
                pass

            try:
                res[id]['existence'] = get_existence(r)
            except:
                pass

            try:
                res[id]['rhea_ids'] = get_reactions(r)
            except:
                pass

        save_json(res, jni_uniprot)
        
        with open('../data/mapping/bad_jni_uniprot_ids.txt', 'w') as f:
            for elt in bad_ids:
                f.write(elt + '\n')

    return res, bad_ids, todo


if __name__ == "__main__":
    from src.utils import save_json, load_json
    import os

    if os.path.exists(jni_uniprot):
        res = load_json(jni_uniprot)
        acquired_ids = list(res.keys())
    else:
        res = {}
        acquired_ids = []

    print("Reading JNi UniProt IDs...")
    with open(up_ids_path, "r") as f:
        uniprot_ids = [line.strip() for line in f]

    if os.path.exists('../data/mapping/bad_jni_uniprot_ids.txt'):
        with open('../data/mapping/bad_jni_uniprot_ids.txt', "r") as f:
            bad_ids = [line.strip() for line in f]
    else:
        bad_ids = []

    
    uniprot_ids = [elt for elt in uniprot_ids if elt not in bad_ids + acquired_ids]
    res, bad, todo = scrape_loop(uniprot_ids=uniprot_ids, batch_size=batch_size, res=res)
    bad_ids += bad

    while todo:
        res, bad, todo = scrape_loop(uniprot_ids=todo, batch_size=batch_size, res=res)
        bad_ids += bad
    
    print("Writing to file.")
    save_json(res, jni_uniprot)

    with open('../data/mapping/bad_jni_uniprot_ids.txt', 'w') as f:
        for elt in bad_ids:
            f.write(elt + '\n')
