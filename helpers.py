import json
import requests
from intermine.webservice import Service
import pandas as pd
import urllib.request as urllib


######################### UPDATE GENE LIST ############################

def get_gene_json():

    url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    data = requests.get(url).text

    gene_list = data.split("\n")[1:-1]

    with open('static/data/genes.json', 'w') as f:
        json.dump(gene_list, f)


########################## QUERY ENRICHER ###############################

def query_enricher(genes: list):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)

    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error adding gene list')

    data = json.loads(response.text)

    return data['userListId']

########################## QUERY KOMP API ###############################



def query_komp(gene: str):

    gene =  gene[0].upper() + (gene[1:]).lower()
    KOMP_URL = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:" + gene
    response = requests.get(KOMP_URL)
    if not response.ok:
        raise Exception('Error analyzing retrieving information')
    data = json.loads(response.text)
    return data



######## QUERY MGI ##############

def query_mgi(gene: str):
    service = Service("https://www.mousemine.org/mousemine/service")

    # Returns the phenotypes (MP terms) associated with the specified  mouse genes
    # or other features.

    template = service.get_template('_Feature_Phenotype')


    rows = template.rows(
    B = {"op": "LOOKUP", "value": gene}
    )
    dict_data = {'data': [dict(row) for row in rows]}

    seen = set()
    res = []
    for d in dict_data['data']:
        t = (d['OntologyAnnotation.ontologyTerm.name'], d['OntologyAnnotation.evidence.publications.pubMedId'])
        if t not in seen:
            seen.add(t)
            res.append(d)
    
    return {'data': res}


######## QUERY GWAS ##############

def query_gwas(gene: str):
    gene = gene.upper()
    url = "https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:" + gene + "&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true"

    df = pd.read_csv(url, sep='\t')

    summerized_data = pd.DataFrame(columns=['gene', 'trait', 'mapped_trait_link', 'count'])

    print(df.columns.values)
    i = 0
    seen_traits = set()
    for index, row in df.iterrows():
        trait = row['MAPPED_TRAIT'].strip()
        if not trait in seen_traits:
            
            if "," in row['MAPPED_TRAIT_URI']:
                uri = row['MAPPED_TRAIT_URI'].split(",")[0]
            else:
                uri = row['MAPPED_TRAIT_URI']
            summerized_data.loc[i] = {'gene': gene, 'trait': trait, 'mapped_trait_link': uri, 'count': 1}
            seen_traits.add(trait)
            i+=1
        else:
            idx = summerized_data.index[summerized_data['trait'] == trait].values[0]
            summerized_data.at[idx,'count'] +=1


    summerized_data = summerized_data.sort_values(by=['count'],ascending=False)
    trait_list = list(summerized_data.T.to_dict().values())

    return {'GWAS_Catalog': trait_list}


########## QUERY SigCom LINCS ###############


def sigcom_up_down_genes(up_list, down_list):
    METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"
    DATA_API = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1/"

    input_gene_set = {
        "up_genes": up_list,
        "down_genes": down_list 
    }

    all_genes = input_gene_set["up_genes"] + input_gene_set["down_genes"]

    payload = {
        "filter": {
            "where": {
                "meta.symbol": {
                    "inq": all_genes
                }
            },
            "fields": ["id", "meta.symbol"]
        }
    }
    res = requests.post(METADATA_API + "entities/find", json=payload)
    entities = res.json()

    for_enrichment = {
        "up_entities": [],
        "down_entities": []
    }

    for e in entities:
        symbol = e["meta"]["symbol"]
        if symbol in input_gene_set["up_genes"]:
            for_enrichment["up_entities"].append(e["id"])
        elif symbol in input_gene_set["down_genes"]:
            for_enrichment["down_entities"].append(e["id"])
    

    payload = {
    "meta": {
        "$validator": "/dcic/signature-commons-schema/v6/meta/user_input/user_input.json",
        **for_enrichment
    },
    "type": "signature"
    }
    res = requests.post(METADATA_API + "user_input", json=payload)
    persistent_id = res.json()["id"]

    return ("https://maayanlab.cloud/sigcom-lincs#/SignatureSearch/Rank/%s"%persistent_id)




def sigcom_gene_set(gene_set):
    METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api/"

    input_gene_set = {
        "genes": gene_set,
    }

    all_genes = input_gene_set["genes"]

    payload = {
        "filter": {
            "where": {
                "meta.symbol": {
                    "inq": all_genes
                }
            },
            "fields": ["id", "meta.symbol"]
        }
    }
    res = requests.post(METADATA_API + "entities/find", json=payload)
    entities = res.json()

    for_enrichment = {
    "entities": [],
    "signatures": [],
    "offset": 0,
    "limit": 0
    }


    for e in entities:
        for_enrichment["entities"].append(e["id"])
    

    payload = {
    "meta": {
        "$validator": "/dcic/signature-commons-schema/v6/meta/user_input/user_input.json",
        **for_enrichment
    },
    "type": "signature"
    }
    res = requests.post(METADATA_API + "user_input", json=payload)

    persistent_id = res.json()["id"]

    return ("https://maayanlab.cloud/sigcom-lincs/#/SignatureSearch/Set/%s"%persistent_id)


