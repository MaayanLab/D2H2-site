import openai
import os
from dotenv import load_dotenv
from functools import lru_cache
import json


load_dotenv()

with open('static/searchdata/processes.json') as f:
    processes = json.load(f)

validation = {"[Gene]": list(processes["[Gene]"].keys()),
              "[GeneSet]": list(processes["[GeneSet]"].keys())}

gene_process_descs = "\n".join(map(lambda output: f'[Gene]->{output} - {processes["[Gene]"][output]["gpt_desc"]}', list(processes["[Gene]"].keys())))
geneset_process_descs = "\n".join(map(lambda output: f'[GeneSet]->{output} - {processes["[GeneSet]"][output]["gpt_desc"]}', list(processes["[GeneSet]"].keys())))

openai.api_key = os.getenv("OPENAI_API_KEY")

def determine_valid(query):
    prompt = f"""
    Based on the query from the user: "{query}"
    Pick an input type from the list: [[Metabolite],[RNA-seq file],[Variant],[Transcript],[Gene],[GeneSet]]
    Additional context: If a gene symbol is included in the input then the input will most likely be [Gene]. A GeneSet is a collection of mulitple genes, thus is the user includes 'genes' or 'gene set' the input will likely be [GeneSet]. If the user provides a query with the word gene in it, they could be asking for a gene as an output which is not valid. 
    If the user is asking for a gene or a gene set as an output, you should respond with [Other].
    Please also provide a confidence score in the range of: [0, 1].
    Respond in this format:
    [confidence score],[Type]
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
        {"role": "system", "content": "You are an assitant meant to process a user query and decide what type of input and output the user is specifiying"},
        {"role": "user", "content": prompt}
            ],
        max_tokens =20,
        temperature=1,
        )

        response = tag_line['choices'][0]['message']['content']
        print(response)
        response = response.split(',')[1]
        if '[Gene]' in response:
            return (True, '[Gene]')
        elif '[GeneSet]' in response:
            return (True, '[GeneSet]')
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
    else: 
        processes_descs = geneset_process_descs
    
    prompt = f"""
    Based on the query from the user: "{query}"
    Additional context: If a gene symbol is included in the input then the input will most likely be [Gene]. A GeneSet is a collection of mulitple genes, thus is the user includes 'genes' or 'gene set' the input will likely be [GeneSet]. If the user provides a query with an input does not match any of the included processes, please use the type: [None]. 

    Pick an input type from the list: [[Gene], [GeneSet]]

    Then, pick a process to be exceuted from below based on the query and the chosen input type:
    {processes_descs}

    Your response must strictly follow the following format with no other text, description or reasoning:
    [Input]->[Output]
    """
    try:
        tag_line = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
        {"role": "system", "content": "You are an assitant meant to process a user query and pick a given workflow"},
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
        model="gpt-3.5-turbo",
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
    
