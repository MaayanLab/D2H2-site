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
processes_descs = f"{gene_process_descs}\n{geneset_process_descs}"

openai.api_key = os.getenv("OPENAI_API_KEY")

@lru_cache()
def find_process(query):
    prompt = f"""
    Based on the query from the user: "{query}"
    Additional context: If a gene symbol or name is included in the input then the input will most likely be [Gene]. A GeneSet is a collection of mulitple genes, thus is the user includes 'genes' or 'gene set' the input will likely be [GeneSet]. 

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
        temperature=.4,
        )

        response = tag_line['choices'][0]['message']['content']
        if '"' in response:
            response = response.split('"')[1]
        response = response.replace(".", "")
        print(response)
        if response in options:
            return {'option': response}
        else:
            return {'option': 'error'}
    except:
        return {'option': 'busy'}
    
    