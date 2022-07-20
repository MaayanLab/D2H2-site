import json
import requests
import pandas as pd
from bs4 import BeautifulSoup
from selenium import webdriver






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
    print(data)
    return data



########################## QUERY KOMP API ######################################

g = 'Car4'

def query_komp(gene: str):
    KOMP_URL = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:" + g
    response = requests.get(KOMP_URL)
    print(response)
    if not response.ok:
        raise Exception('Error analyzing retrieving information')
    data = json.loads(response.text)
    print(data)
    return data


#query_enricher('HLA-A')
#query_komp(g)

def query_archs4(gene: str):

    url = 'https://maayanlab.cloud/archs4/search/genepage.php?search=go&gene=' + gene
    df = pd.read_html(url)

    print(df)



#query_archs4('HLA-A')

#downloadAnchorElem

r = requests.get("https://maayanlab.cloud/archs4/gene/A2M").content



driver = webdriver.Chrome()

url = "https://maayanlab.cloud/archs4/gene/A2M"

driver.get(url)

document.querySelector("#downloadAnchorElem")


button = driver.find_element_by_path("/html/body/div[1]/h2/span/a")

soup = BeautifulSoup(r, features="lxml")
print(soup.find("a", { "id" : "downloadAnchorElem" }))

#downloadAnchorElem