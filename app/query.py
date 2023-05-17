import openai
import os
from dotenv import load_dotenv
from functools import lru_cache
import json


load_dotenv()

with open('processes.json') as f:
    processes = json.load(f)

validation = {"[Gene]": list(processes["[Gene]"].keys()),
              "[GeneSet]": list(processes["[GeneSet]"].keys())}


openai.api_key = os.getenv("OPENAI_API_KEY")

@lru_cache()
def find_process(query):
    prompt = f"""
    Based on the query from the user: "{query}"
    Additional context: If a gene symbol or name is included in the input then the input will most likely be [Gene]. A GeneSet is a collection of mulitple genes, thus is the user includes 'genes' the input will likely be [GeneSet]. 

    Pick an input type from the list: [[Gene], [GeneSet]]

    Then, pick a process to be exceuted from below based on the query and the chosen input type:
    [Gene]->[Signatures] - In what diabetes related GEO signatures up or down regulate the expression of my gene?
    [Gene]->[Expression] - In what cells and tissues is my gene expressed?
    [Gene]->[Perturbations] - Under what conditions or perturbations is my gene regulated?
    [Gene]->[TFs] - What are the transcription factors that regulate my gene?
    [Gene]->[Traits] - Is my gene associated with traits in human GWAS?
    [Gene]->[Correlation] - Is my gene correlated with other genes?
    [Gene]->[Knockout] - Is there a knockout mouse for my gene and does it show any phenotypes?
    [GeneSet]->[Enrichment] - In which annotated gene sets is my gene set enriched?
    [GeneSet]->[Signatures] - In what diabetes signatures is my gene set enriched/overlaps with?
    [GeneSet]->[TFs] - What transcription factors regulate my gene set?
    [GeneSet]->[Kinases] - What kinases regulate my gene set?
    [GeneSet]->[L1000] - What are the LINCS L1000 small molecules and genetic perturbations that likely up- or down-regulate the expression of my gene set?

    Your response must strictly follow the following format with no other text, description or reasoning:
    [Input]->[Output]
    """

    tag_line = openai.ChatCompletion.create(
    model="gpt-3.5-turbo",
    messages=[
    {"role": "system", "content": "You are an assitant meant to process a user query and pick a given workflow"},
    {"role": "user", "content": prompt}
        ],
    max_tokens =20,
    temperature=.4,
    )

    response = tag_line['choices'][0]['message']['content']
    print(response)
    try:
        input, output = response.split('->')
        if output in validation[input]:
            return {"response": 0, "input": input, "output": output}
        else:
            return {"response": 1, "error": "Process does not exist"}
    except:
        return {"response": 1, "error": "Parsing error"}
    
def select_option(response, options):
    prompt = f"""
    Based on the reponse from the user: "{response}"
    pick an option from the list of options: "{options}"
    Your response must only include the extact string from the list of options with no other text, description, reasoning or punctuation.
    """

    tag_line = openai.ChatCompletion.create(
    model="gpt-3.5-turbo",
    messages=[
    {"role": "system", "content": "You are an assitant meant to process a user response and pick from a predefined list of options"},
    {"role": "user", "content": prompt}
        ],
    max_tokens =20,
    temperature=.4,
    )

    response = tag_line['choices'][0]['message']['content']
    if '"' in response:
        response = response.split('"')[1]
    response = response.replace(".", "")
    
    if response in options:
        return {'option': response}
    else:
        return {'option': 'error'}