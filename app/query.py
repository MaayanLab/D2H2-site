import openai
import os
from dotenv import load_dotenv
from functools import lru_cache
import json


load_dotenv()

with open('static/searchdata/processes.json') as f:
    processes = json.load(f)

validation = {"[Gene]": list(processes["[Gene]"].keys()),
              "[GeneSet]": list(processes["[GeneSet]"].keys()),
              "[Term]": list(processes["[GeneSet]"].keys())
              }

gene_process_descs = "\n".join(map(lambda output: f'[Gene]->{output} - {processes["[Gene]"][output]["gpt_desc"]}', list(processes["[Gene]"].keys())))
geneset_process_descs = "\n".join(map(lambda output: f'[GeneSet]->{output} - {processes["[GeneSet]"][output]["gpt_desc"]}', list(processes["[GeneSet]"].keys())))


openai.api_key = os.getenv("OPENAI_API_KEY")

def determine_valid(query):
    prompt = f"""
    Based on the query from the user: "{query}"
    Pick an input type from the list: [[Metabolite],[RNA-seq file],[Variant],[Transcript],[Gene],[GeneSet],[Term],[Differential Expression],[Study Metadata],[Other]]
    Additional context: 
    A GeneSet is a collection of mulitple genes.
    If the user is asking for associated or upregulated genes or a gene set as an output, you should respond with [Term]. 
    Also select [Term] if the user appears to be asking a general question about a disease, or in general any biomedical term.
    If the user provides a term and a gene, then you should respond with [Term].
    If the the user's query is not relevant to any type in the list, or in general the question is not relevant in a biomedical context, then please respond with [Other].
    Respond in this format including only 1 type:
    [Type]
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to process a user query and decide what type of input and output the user is specifiying"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=.1,
        )

        response = tag_line['choices'][0]['message']['content']
        print(response)
        if '[Gene]' in response:
            return (True, '[Gene]')
        elif '[GeneSet]' in response:
            return (True, '[GeneSet]')
        elif '[Term]' in response:
            return (True, '[Term]')
        elif '[Study Metadata]' in response:
            return (True, '[Study Metadata]')
        else: 
            return (False, '[None]')
    except:
        return (False, '[None]')
    


@lru_cache()
def find_process(query):
    valid, input_type = determine_valid(query)
    if not valid:
        return {"response": 1, "error": "input"}
    if input_type == '[Gene]':
        processes_descs = gene_process_descs
    elif input_type == '[GeneSet]': 
        processes_descs = geneset_process_descs
    elif input_type == '[Term]':
        res = identify_search_term(query)
        if 'term' in res:
            return {"response": 0, "input": '[Term]', "output": '[Search]', "term": res['term']}
        return {"response": 1, "error": "busy"}
    elif input_type == '[Study Metadata]':
        res = identify_search_term(query)
        if 'term' in res:
            return {"response": 0, "input": '[Study Metadata]', "output": '[Search]', "term": res['term']}
        return {"response": 1, "error": "busy"}

    
    prompt = f"""
    Based on the query from the user: "{query}"
    Pick one of the following processes based on their description's similarity to the user's query:
    {processes_descs}
    Your response must strictly follow the following format with no other text, description or reasoning:
    [Input]->[Output]
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to process a user query and pick the relevant input type from the provided list of options."},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=.1,
        )

        response = tag_line['choices'][0]['message']['content']
        response = response.replace('(|)', '')
        try:
            input, output = response.split('->')
            if output in validation[input]:
                return {"response": 0, "input": input, "output": output}
            else:
                return {"response": 1, "error": "Process does not exist"}
        except:
            return {"response": 1, "error": "Parsing error"}
    except:
        return {"response": 1, "error": "busy"}

def select_option(response, options):
    prompt = f"""
    Based on the reponse from the user: "{response}"
    pick an option from the list of options: "{options}"
    Your response must only include the extact string from the list of options (not from the user response) with no other text, description, reasoning or punctuation.
    If the user response does not match one of the options or if they are asking an additional question then respond: "None"
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to process a user response and pick from a predefined list of options"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=0,
        )

        response = tag_line['choices'][0]['message']['content']
        if '"' in response:
            response = response.split('"')[1]
        response = response.replace(".", "")
        if response in options:
            return {'option': response}
        else:
            return {'option': 'error'}
    except:
        return {'option': 'busy'}



def infer_gene(gene):
    prompt = f"""
    Based on text from the user: "{gene}"
    Respond with three comma-separated valid Entrez gene symbols (either human or mouse) that most closely resemble them the user's original input. Only include the three gene symbols with no other text or reasoning."""
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to predict which gene a user was asking for based on a misspelled or unrecognized entry"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=0,
        )

        response = tag_line['choices'][0]['message']['content']
        print(response)
        gs = response.split(',')
        gs = list(map(lambda x: x.strip(), gs))
        print(gs)
        return {'genes': gs}
    except:
        return {'option': 'busy'}
    

def identify_search_term(query):
    prompt = f"""
    Based on text from the user: "{query}"
    Respond with the biomedical term(s) the user included in their question. Only include the term which should be used to serach with and no other text or reasoning. Genes or geneset are general terms and should not be included."""
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to select the biomedical term from the user's query"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=0,
        )

        response = tag_line['choices'][0]['message']['content'].strip()
        print(response)
        return {'term': response}
    except:
        return {'option': 'busy'}
    

    
def determine_association(user_query): 
    prompt = f"""
    You will be provided with text delimited by triple quotes. 
    A functional term is a gene, biological pathway, biological function, drug, phenotype or disease contained in the
    query.
    If it is a query asking about a gene set related to two or more functional terms, the response should be ['Association', [terms in association]]. Otherwise if the query is
    not asking about a gene set related to two or more functional terms, the output should be [False]. 
        ```{user_query}``` 
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to select the biomedical term from the user's query"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=0,
        )

        response = tag_line['choices'][0]['message']['content'].strip()
        print(response)
        return {'term': response}
    except:
        return {'option': 'busy'}


def select_args(query):
    prompt = f"""
    Process the query delimited by triple backticks and return only one most important term. 
    This term should be the gene, list of genes, 
    biological pathway, biological function, drug, phenotype or disease contained in the
    query. The output should not include any prepositions and only one term in the format: [response]. 
    For example, if the query is 'Provide all the differentially expressed 
    genes in alopecia', the output should be [alopecia].
    Another example, if the query is 'Imagine a geneset of all genes that are mentioned with hair loss', 
    the output is should be [hair loss]
 
    ```{query}```
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-4",
        messages=[
        {"role": "system", "content": "You are an assitant meant to select the biomedical term from the user's query"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=0,
        )

        response = tag_line['choices'][0]['message']['content'].strip()
        print(response)
        return {'term': response}
    except:
        return {'option': 'busy'}