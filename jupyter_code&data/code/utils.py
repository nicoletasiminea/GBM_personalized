import pandas as pd
from pathlib import Path
from collections import defaultdict
import os
import json
import requests

# read the mutations, the up-regulated genes and the down-regulated genes for each case
m_table = pd.read_excel("../data/Mutations.xlsx")

mutation_t = (
    m_table
    .groupby(m_table.columns[0])[m_table.columns[1]]
    .agg(lambda x: ",".join(sorted(set(x))))
    .to_dict()
)

for patient in mutation_t:
    mutation_t[patient]=mutation_t[patient].split(",")

def load_diff(path_input):

    path=Path.cwd()/path_input
    if not path.exists():
        raise FileNotFoundError(f"Path not found: {path.resolve()}")

    diff_dict = defaultdict(list)
    for file in path.glob("*.txt"):  
        key = file.stem 
        with open(file, "r") as f:
            for line in f:
                clean_line = line.strip()
                diff_dict[key].append(clean_line)

    diff_dict = dict(diff_dict)
    return diff_dict

down_dict = load_diff("../data/DownProteins")
up_dict = load_diff("../data/UpProteins")

# read patients lists (primary, secondary, and all for primary)
patient_ids = pd.read_csv("../data/PatientIds.txt", header=None).iloc[:, 0].tolist()


#define the lists for primary tumors and for recurrent tumors (we exluded those with a high number of mutations)
list_recurrent=["0317a370-6a1e-44bd-bfd0-81fd06ec56fb",
                "c12052d7-f5ae-4fd5-a36b-ac0ce21179a6",
                "f8e7cbd2-b54a-4f30-8fc1-a5c80bfad8b6",
                "57c3b272-7dcb-4937-83c1-02a3694e701f",
                "ddd4c4a7-6057-4244-aa4d-ccffe6e1c95a",
                "928160c7-f711-4b23-832c-8c12aea9fe85",
                "f0f7d061-ace0-4dbd-90ec-df8f5985a7a2",
                "d44d370d-d86e-4006-9530-ab3442cb2848",
                "6dfd7934-3161-4384-be4d-93b14080ece8",
                "fe179a60-5e81-48c1-b437-7c38ed910ba8",
                "1142d18f-9cd9-41b4-a9f2-975730cb713a",
                "b4919075-90a9-4c8a-82a0-a0b5f9944ad3",
                "8223f42b-434d-4966-a123-3b8a947884bc"]

list_all_excluded=["all","ac22ae42-0d53-4b2e-ad6c-f4c98f9be040","d638d8cf-7276-49aa-b2f7-81eb10e7b33d","f5bb4817-1727-41f9-b7f6-3a4e1818c2df"]
list_primary_recurrent=[x for x in patient_ids if x not in list_all_excluded]
list_primary = [x for x in list_primary_recurrent if x not in list_recurrent]

# access EnrichR and take the results from enrichment analysis (after Enrichir website)
def enrichent_results(list_genes, filename, gene_library):
    str_genes = ""
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    str_genes = '\n'.join(list_genes)
    description = 'Example 1'
    payload = {
        'list': (None, str_genes),
        'description': (None, description)
        }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    user_list_id = data['userListId'] 
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    filename = filename
    gene_library = gene_library 

    url = ENRICHR_URL + query_string % (user_list_id, filename, gene_library)
    response = requests.get(url, stream=True)

    with open(filename + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024): 
            if chunk:
                f.write(chunk)


