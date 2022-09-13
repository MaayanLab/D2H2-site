from array import array
from flask import Flask, render_template, request
from waitress import serve
import os
import json
import plotly
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import GEOparse
import ftfy
from functools import lru_cache
import pickle
from helpers import *
from twitterauth import update_tweets_table

create_meta = False

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def home():
	update_tweets_table()
	return render_template('home.html')

@app.route("/about", methods=['GET', 'POST'])
def about():
    return render_template("about.html")


@app.route("/singlegene/<gene>", methods=['GET', 'POST'])
def singlegene(gene):
    return render_template("singlegene.html", gene=gene)

@app.route("/singlegene", methods=['GET', 'POST'])
def singlegene_home():
    return render_template("singlegene.html")

@app.route("/geneset", methods=['GET', 'POST'])
def geneset():
    return render_template("geneset.html")

@app.route('/scg', methods=['GET', 'POST'])
def scg():
	return render_template('scg.html')

@app.route('/resources', methods=['GET', 'POST'])
def resources():
	return render_template('resources.html')

@app.route('/downloads', methods=['GET', 'POST'])
def downloads():
	return render_template('downloads.html')

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

@app.route('/getresources',  methods=['GET','POST'])
def resources_api():
	table = get_resources()

	return {'resources': table}

@app.route('/getdownloads',  methods=['GET','POST'])
def downloads_api():
	table = get_downloads()

	return {'downloads': table}


@app.route('/gettweets',  methods=['GET','POST'])
def tweets_api():
	table = get_tweets()
	return {'tweets': table}

@app.route('/getworkflows',  methods=['GET','POST'])
def workflows_api():
	table = get_workflows()

	return {'workflows': table}

@app.route('/gettfs',  methods=['GET','POST'])
def gettfs():
	gene = request.form['gene']

	if gene.strip() == '':
		return {'data': []}
	result = query_enricher(gene)

	return {'data': result}

@app.route('/getexample',  methods=['GET','POST'])
def getexample():

	with open('static/searchdata/example_list.txt') as f:
		text = f.read()

	return {'genes': text, 'description': "GSE136134 Ctrl-vs-Insulin 24hrs Human BulkRNAseq hiPSCs_down"}

@app.route('/getdiabetesenrich',  methods=['GET','POST'])
def getdiabetesenrich():

	genes = request.form['genelist']
	description = request.form['description']

	data = query_enricher_diabetes(genes, description)

	return {'data': data}

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


def get_metadata(geo_accession, species_folder):
	"""
	Gets the metadata for the GEO Series.
	Input: geo_accession, a string representing the GEO Accession ID (possibly with the GPL)
		of the form GSEXXXXX(-GPLXXXX)
	Output: A dictionary representing the metadata parsed by GEOparse, with keys being labels and values being that value of the
		metadata label.
	"""
    # Case where I have GPL in geo_accession
	if "-" in geo_accession:
		geo_accession_num = geo_accession.split("-", 1)[0]
		gpl_num = geo_accession.split("-", 1)[1]
	else:
		geo_accession_num = geo_accession
    # Get gse from GEO
	if f'{geo_accession}_family.soft.gz' not in os.listdir(f'./static/data/{species_folder}/{geo_accession}'):
		gse = GEOparse.get_GEO(geo = geo_accession_num, destdir = f'./static/data/{species_folder}/{geo_accession}', silent=True)
	else:
		gse = GEOparse.get_GEO(filepath = f'./static/data/{species_folder}/{geo_accession}/{geo_accession}_family.soft.gz', silent=True)

	if "-" in geo_accession:
		gse.metadata['cur_gpl'] = [gpl_num]
	else:
		gse.metadata['cur_gpl'] = [gse.metadata['platform_id'][0]]

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
	
	
	metadata_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	metadata_dataframe = pd.read_csv(metadata_file, sep='\t')
	gse.metadata['numsamples'] = metadata_dataframe.shape[0] - 1


	return gse.metadata

ignore_list = ['mouse_matrix_v11.h5', 'human_matrix_v11.h5', '.DS_Store', 'allgenes.json']

def species_to_studies(path):
	species_mapping = {}
	for species_name in os.listdir(path):
		if species_name not in ignore_list:
			species_mapping[species_name] = []
			for study_name in os.listdir(os.path.join(path, species_name)):
				if study_name not in ignore_list:
					species_mapping[species_name].append(study_name)
	return species_mapping

def sort_studies(species_mapping):
	species_url_mapping = {}
	for species_folder, studies_list in species_mapping.items():
		if species_folder not in ignore_list:
			species = folder_to_url[species_folder]
			species_url_mapping[species] = sorted(studies_list, key=lambda x: (int(x.split('-', 1)[0][3:])))
	return species_url_mapping

url_to_folder = {"human": "human", "mouse": "mouse"}
folder_to_url = {folder:url for url, folder in url_to_folder.items()}

species_mapping = species_to_studies('static/data')
species_mapping = sort_studies(species_mapping)
# print(species_mapping)



if create_meta:
	gse_metadata = {}
	for species, geo_accession_ids in species_mapping.items():
		gse_metadata[species] = {}
		for geo_accession in geo_accession_ids:
			gse_metadata[species][geo_accession] = get_metadata(geo_accession, url_to_folder[species])
	with open('static/searchdata/metadata-v1.pickle', 'wb') as f:
		pickle.dump(gse_metadata, f, protocol=pickle.HIGHEST_PROTOCOL)
else:
	with open('static/searchdata/metadata-v1.pickle', 'rb') as f:	
		gse_metadata = pickle.load(f)



study_to_species = {study:species_name for species_name, studies_metadata in gse_metadata.items() for study in studies_metadata.keys()}


@app.route("/<species_or_gse>", methods=['GET', 'POST'])
def species_or_viewerpg(species_or_gse):
	# test if species
	if species_or_gse in gse_metadata:
		return render_template('species.html', species=species_or_gse, gse_metadata=gse_metadata, species_mapping=species_mapping)
	# test if gsea
	elif species_or_gse in study_to_species:
		geo_accession = species_or_gse
		species = study_to_species[geo_accession]
		species_folder = url_to_folder[species]
		metadata_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
		metadata_dataframe = pd.read_csv(metadata_file, sep='\t')
		metadata_dict = metadata_dataframe.groupby('Group')['Condition'].apply(set).to_dict()
		metadata_dict_samples = metadata_dataframe.groupby('Condition')['Sample_geo_accession'].apply(list).to_dict()
		sample_dict = {}
		for key in metadata_dict_samples.keys():
			sample_dict[key] = {'samples':metadata_dict_samples[key], 'count': len(metadata_dict_samples[key])}

		return render_template('viewer.html', metadata_dict=metadata_dict, metadata_dict_samples=sample_dict, geo_accession=geo_accession, gse_metadata=gse_metadata, species=species, species_mapping=species_mapping)
	else:
		return render_template('error.html', gse_metadata=gse_metadata, species_mapping=species_mapping)





#############################################
########## 1. Genes
#############################################


# this will likely stay the same.

@app.route('/api/genes/<geo_accession>')

# this will likely stay the same.
@lru_cache(maxsize=None)
def genes_api(geo_accession):
	if geo_accession == 'human':
		with open('static/searchdata/t2d-human.json', 'r') as f:
			human_genes = json.load(f)
		genes_json = json.dumps([{'gene_symbol': x} for x in human_genes['human_genes']])
	elif geo_accession == 'mouse':
		with open('static/searchdata/t2d-mouse.json', 'r') as f:
			mouse_genes = json.load(f)
		genes_json = json.dumps([{'gene_symbol': x} for x in mouse_genes['mouse_genes']])
	else:
		species = study_to_species[geo_accession]
		species_folder = url_to_folder[species]

		# Get genes json
		expression_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'
		expression_dataframe = pd.read_csv(expression_file, index_col = 0, sep='\t')

		genes_json = json.dumps([{'gene_symbol': x} for x in expression_dataframe.index])

	# Return
	return genes_json

#############################################
########## 2. Plot
#############################################

@app.route('/api/plot/<geo_accession>', methods=['GET', 'POST'])

def plot_api(geo_accession):
	"""
	Inputs:
	- expression_dataframe, a tsv file representing a matrix with rows having indices representing gene symbols and columns with indices representing Sample ID's. 
		Each entry represents the expression value for a gene in a sample.
	- metadata_dict, a JSON file with the following format: {group_name:{condition_names:[sample_name, ...]}}
	Outputs: 
	- plotly_json, a serialized JSON formatted string that represents the data for the boxplot to be plotted, used by boxplot() function in scripts.js.
	"""
	species = study_to_species[geo_accession]
	species_folder = url_to_folder[species]

	assay = gse_metadata[species][geo_accession].get('type')[0]

	expression_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'
	metadata_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	expression_dataframe = pd.read_csv(expression_file, index_col = 0, sep='\t')
	metadata_dataframe = pd.read_csv(metadata_file, sep='\t')
	
	# Get data
	data = request.json
	gene_symbol = data['gene_symbol']

	conditions = data['conditions']
	# print(conditions)

	# Create a new dataframe that maps each sample to its condition
	melted_dataframe = expression_dataframe.loc[gene_symbol].rename('expr_vals').rename_axis('Sample_geo_accession').reset_index().merge(metadata_dataframe, on="Sample_geo_accession")


	# Get plot dataframe
	plot_dataframe = melted_dataframe.groupby('Condition')['expr_vals'].agg([np.mean, np.std, lambda x: list(x)])#.rename(columns={'<lambda>': 'points'})#.reindex(conditions)
	plot_dataframe = plot_dataframe.rename(columns={plot_dataframe.columns[-1]: 'points'})
	# print(plot_dataframe)

	# Initialize figure
	fig = go.Figure()

	# print(plot_dataframe)
	
	# Loop
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
	return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

@app.route('/api/volcano', methods=['GET', 'POST'])
def plot_volcano_api():
	request.form
	gene = request.form["gene"]
	species = request.form["species"]
	json_item_plot = send_plot(species, gene)
	# Return
	return json.dumps(json_item_plot)


########## 3. Conditions ##########

@app.route('/api/conditions/<geo_accession>')

def conditions_api(geo_accession):
	species = study_to_species[geo_accession]
	species_folder = url_to_folder[species]
	metadata_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	metadata_dataframe = pd.read_csv(metadata_file, sep='\t')
	conditions = set(metadata_dataframe['Condition'])
	return json.dumps([{'Condition': x} for x in conditions])

#############################################
########## Conditions to Samples
#############################################
@app.route('/api/samples/<geo_accession>')

def samples_api(geo_accession):
	species = study_to_species[geo_accession]
	species_folder = url_to_folder[species]

	metadata_file = 'static/data/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	metadata_dataframe = pd.read_csv(metadata_file, sep='\t')

	conditions_mapping = metadata_dataframe.groupby('Condition')['Sample_geo_accession'].apply(list).to_dict()

	return json.dumps(conditions_mapping)


@app.route('/api/data',  methods=['GET', 'POST'])
def get_study_data():

	response_json = request.get_json()
	geo_accession = response_json['gse']
	control = response_json['control']
	perturb = response_json['perturb']

	metadata_file = 'static/data/' + 'human' + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	expression_file = 'static/data/' + 'human' + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'

	

	with open(metadata_file, 'r') as f:
		meta_data = f.read()

	with open(expression_file, 'r') as f:
		expression_data = f.read()

	selected_conditions = []
	for i, row in enumerate(meta_data.split('\n')):
		if i == 0:
			selected_conditions.append(row)
		
		entries = row.split('\t')
		if len(entries) != 3:
			break
		if entries[1] == control or entries[1] == perturb:
			selected_conditions.append(row)


	selected_meta = '\n'.join(selected_conditions)

	data_dict = {'meta': selected_meta, 'expression': expression_data}

	return data_dict
	








#######################################################
#######################################################
########## 3. Run App
#######################################################
#######################################################
if __name__ == "__main__":
	
	#serve(app, host="0.0.0.0", port=5000)
	app.run(debug=True, host="0.0.0.0")

