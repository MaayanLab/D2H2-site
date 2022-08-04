from flask import Flask, render_template, url_for, request, session
import os
import json
import plotly
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import GEOparse
import ftfy
from functools import lru_cache
from helpers import *


app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def home():
	return render_template('home.html')

@app.route("/about", methods=['GET', 'POST'])
def about():
    return render_template("about.html")


@app.route("/singlegene", methods=['GET', 'POST'])
def singlegene():
    return render_template("singlegene.html")

@app.route("/geneset", methods=['GET', 'POST'])
def geneset():
    return render_template("geneset.html")

@app.route('/scg', methods=['GET', 'POST'])
def scg():
	return render_template('scg.html', methods=['GET', 'POST'])

@app.route('/getgwas', methods=['GET','POST'])
def get_gwas():

    gene = request.form['gene']

    if gene == '':
        return {'GWAS_Catalog':[]}

    result = query_gwas(gene)

    return result

@app.route('/getkomp', methods=['GET','POST'])
def get_mgi():

    gene = request.form['gene']
    result = query_mgi(gene)

    return result


@app.route('/getsigcom',  methods=['GET','POST'])
def get_sigcom():

	gene_lists = request.get_json()["genes"]
	if len(gene_lists) == 2:
		res= {'url': sigcom_up_down_genes(gene_lists[0], gene_lists[1])}
	elif len(gene_lists) == 1:
		res= {'url':sigcom_gene_set(gene_lists[0])}
	else: 
		res= {'url': "https://maayanlab.cloud/sigcom-lincs/#/"}
	
	return res

#############################################
########## 2. Data
#############################################ow toow
##### 1. Files #####


def fix_para_text(gse_metadata_dict, attribute):
	"""
	Modifies input dictionary
	"""
	if (attribute in gse_metadata_dict) and (len(gse_metadata_dict.get(attribute, [])) != 0):
		for idx in range(0, len(gse_metadata_dict[attribute])):
			gse_metadata_dict[attribute][idx] = ftfy.fix_text(gse_metadata_dict[attribute][idx])


def get_metadata(geo_accession, organ_folder):
    """
    Gets the metadata for the GEO Series.
    Input: geo_accession, a string representing the GEO Accession ID (possibly with the GPL)
        of the form GSEXXXXX(-GPLXXXX)
    Output: A dictionary representing the metadata parsed by GEOparse, with keys being labels and values being that value of the
        metadata label.
    """
    organisms = {}
    # Case where I have GPL in geo_accession
    if "-" in geo_accession:
        geo_accession_num = geo_accession.split("-", 1)[0]
        gpl_num = geo_accession.split("-", 1)[1]
    else:
        geo_accession_num = geo_accession
    # Get gse from GEO
    gse = GEOparse.get_GEO(geo = geo_accession_num, destdir = f'./static/data/{organ_folder}/{geo_accession}', silent=True)

    if "-" in geo_accession:
        gse.metadata['cur_gpl'] = [gpl_num]

    gse.metadata['all_gpls'] = []

    for gpl in gse.gpls.values():
        gse.metadata['all_gpls'].append(gpl.metadata)

    # Add link to gse (put in a list for consistency)
    gse.metadata['gse_link'] = ['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + geo_accession_num]

    # Fix text
    fix_para_text(gse.metadata, 'title')
    fix_para_text(gse.metadata, 'summary')
    fix_para_text(gse.metadata, 'overall_design')

    # Add Pubmed ID, if exists
    if 'pubmed_id' in gse.metadata and (len(gse.metadata.get('pubmed_id', [])) != 0):
        gse.metadata['pubmed_link'] = ["https://pubmed.ncbi.nlm.nih.gov/" + gse.metadata['pubmed_id'][0]]
    
    return gse.metadata

ignore_list = ['mouse_matrix_v11.h5', 'human_matrix_v11.h5', '.DS_Store', 'allgenes.json']

def organs_to_studies(path):
	organs_mapping = {}
	for organ_name in os.listdir(path):
		if organ_name in ignore_list:
			continue
		organs_mapping[organ_name] = []
		
		for study_name in os.listdir(os.path.join(path, organ_name)):
			if study_name in ignore_list:
				continue
			organs_mapping[organ_name].append(study_name)
	return organs_mapping

def sort_studies(organs_mapping):
	for organ, studies_list in organs_mapping.items():
		
		organs_mapping[organ] = sorted(studies_list, key=lambda x: (int(x.split('-', 1)[0][3:])))
	return organs_mapping


organs_mapping = organs_to_studies('static/data')
print(organs_mapping)

organs_mapping = sort_studies(organs_mapping)


gse_metadata = {}
for organ, studies in organs_mapping.items():
	gse_metadata[organ] = {}
	for study in studies:
		gse_metadata[organ][study] = get_metadata(study, organ)


organ_folder_name = {"human": "human", "mouse": "mouse"}
site_metadata = {"human": "human", "mouse": "mouse"}
organ_url_name = {"human": "human", "mouse": "mouse"}


@app.route("/<organ_name>")
def organ_page(organ_name):
	if organ_name in organ_folder_name:
		organ_folder = organ_folder_name[organ_name]
		gse_list = organs_mapping[organ_folder]
		gse_metadata_organ = gse_metadata[organ_folder]
		return render_template('species.html', gse_list = gse_list, gse_metadata = gse_metadata_organ, organ_name=organ_name, organs_mapping=organs_mapping, organ_url_name=organ_url_name)
	
	else:
		return render_template('error.html', organs_mapping=organs_mapping, organ_url_name=organ_url_name)

		


@app.route('/<organ_name>/<geo_accession>')
def gene_explorer(organ_name, geo_accession):
	if organ_name in organ_folder_name:
		organ_folder = organ_folder_name[organ_name]
		gse_list = organs_mapping[organ_folder]
		gse_metadata_organ = gse_metadata[organ_folder]

		if geo_accession not in gse_list:
			return render_template('error.html', gse_list = gse_list)

		# Data
		metadata_file = 'static/data/' + organ_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.json'
		
		# Groups
		with open(metadata_file) as openfile:
			metadata_dict = json.load(openfile)
		
		# Supplementary Metadata
		geo_meta = gse_metadata_organ[geo_accession]

		# Return
		return render_template('viewer.html', metadata_dict=metadata_dict, os=os, geo_accession=geo_accession, geo_meta = geo_meta, gse_list=gse_list, organ_name=organ_name, organs_mapping=organs_mapping, organ_url_name=organ_url_name)#, sample_dataframe=sample_dataframe, conditions_dict=conditions_dict)
	
	else:
		return render_template('error.html', organs_mapping=organs_mapping, organ_url_name=organ_url_name)



##################################################
########## 2.2 APIs
##################################################

#############################################
########## 1. Genes
#############################################

@app.route('/api/genes/<organ_name>/<geo_accession>')

# this will likely stay the same.
@lru_cache(maxsize=None)
def genes_api(organ_name, geo_accession):
	organ_folder = organ_folder_name[organ_name]

	# Get genes json
	expression_file = 'static/data/' + organ_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'
	expression_dataframe = pd.read_csv(expression_file, index_col = 0, sep='\t')

	genes_json = json.dumps([{'gene_symbol': x} for x in expression_dataframe.index])

	# Return
	return genes_json

#############################################
########## 2. Plot
#############################################

@app.route('/api/plot/<organ_name>/<geo_accession>', methods=['GET', 'POST'])
def plot_api(organ_name, geo_accession):
	"""
	Inputs:
	- expression_dataframe, a tsv file representing a matrix with rows having indices representing gene symbols and columns with indices representing Sample ID's. 
		Each entry represents the expression value for a gene in a sample.
	- metadata_dict, a JSON file with the following format: {group_name:{condition_names:[sample_name, ...]}}
	Outputs: 
	- plotly_json, a serialized JSON formatted string that represents the data for the boxplot to be plotted, used by boxplot() function in scripts.js.
	"""
	organ_folder = organ_folder_name[organ_name]
	assay = gse_metadata[organ_folder][geo_accession].get('type')[0]
	
	expression_file = 'static/data/' + organ_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'
	metadata_file = 'static/data/' + organ_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.json'
	expression_dataframe = pd.read_csv(expression_file, index_col = 0, sep='\t')

	with open(metadata_file) as openfile:
		metadata_dict = json.load(openfile)
	
	# Get samples to conditions mapping
	samp_to_cond = {}

	for conditions_dict in metadata_dict.values():
		for condition, samples_list in conditions_dict.items():
			for sample in samples_list:
				samp_to_cond[sample] = condition


	# Get data
	data = request.json
	gene_symbol = data['gene_symbol']
	# print(gene_symbol)
	conditions = data['conditions']

	# Create a new dataframe that maps each sample to its condition
	if gene_symbol not in expression_dataframe.index.values:
		return "error"
	melted_dataframe = expression_dataframe.loc[gene_symbol].rename('expr_vals').rename_axis('Sample').reset_index()
	melted_dataframe['Condition'] = melted_dataframe['Sample'].map(samp_to_cond)

	# Get plot dataframe
	plot_dataframe = melted_dataframe.groupby('Condition')['expr_vals'].agg([np.mean, np.std, lambda x: list(x)])#.rename(columns={'<lambda>': 'points'})#.reindex(conditions)
	plot_dataframe = plot_dataframe.rename(columns={plot_dataframe.columns[-1]: 'points'})
	# print(plot_dataframe)

	# Initialize figure
	fig = go.Figure()

	for condition in conditions:
		if len(plot_dataframe.loc[condition, 'points']) > 1:
			fig.add_trace(go.Box(name=condition, y=plot_dataframe.loc[condition, 'points'], boxpoints='all', pointpos=0))
		else:
			fig.add_trace(go.Scatter(name=condition, x=[condition], y=plot_dataframe.loc[condition, 'points']))
	
	# Determine y-axis expression string

	# Layout

	if assay == 'Expression profiling by array':
		y_units = 'Expression<br>RMA Normalized'
	else:
		y_units = 'Expression'

	fig.update_layout(
		title = {'text': gene_symbol+' gene expression', 'x': 0.5, 'y': 0.85, 'xanchor': 'center', 'yanchor': 'top'},
		xaxis_title = 'Condition',
		yaxis_title = y_units,
		showlegend = False
	)

	# Return
	return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

#######################################################
#######################################################
########## 3. Run App
#######################################################
#######################################################
if __name__ == "__main__":
	app.run(debug=True, host='0.0.0.0')

