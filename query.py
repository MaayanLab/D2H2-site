import json
import requests


def query_enricher(gene: str):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join([
        gene
    ])

    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    gene_set_library = 'GWAS_Catalog_2019'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    return data



########################## QUERY KOMP API ######################################



def query_komp(gene: str):

    gene =  gene[0].upper() + (gene[1:]).lower()
    KOMP_URL = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:" + gene
    response = requests.get(KOMP_URL)
    if not response.ok:
        raise Exception('Error analyzing retrieving information')
    data = json.loads(response.text)
    return data

