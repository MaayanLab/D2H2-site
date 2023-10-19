from flask import Flask, render_template, request
import os
import ast
import json
import s3fs
import plotly
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import GEOparse
import ftfy
from functools import lru_cache
import pickle
import anndata
from dotenv import load_dotenv
load_dotenv()

from helpers import *
from dge import *
from query import *
from log_chats import *
from log_suggested_study import *

#Added the route for s3 bucket
endpoint = os.environ.get('ENDPOINT', 'https://d2h2.s3.amazonaws.com/')
base_url = os.environ.get('BASE_URL', 'data')
ROOT_PATH = os.environ.get('ROOT_PATH', '/')
BASE_PATH = os.environ.get('BASE_PATH', 'maayanlab.cloud')
DEBUG = os.environ.get('DEBUG', True).lower() in ('true', '1', 't')
UPDATE_STUDIES = os.environ.get('UPDATE_STUDIES', False).lower() in ('true', '1', 't')

s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})


app = Flask(__name__, static_url_path=ROOT_PATH + 'static')


#### PAGES ####
month_dict = {"Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04", "May": "05", "Jun": "06", "Jul": "07", "Aug": "08", "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12"}
@app.route(ROOT_PATH, methods=['GET', 'POST'])
def home():
	return render_template('home.html', base_path=BASE_PATH, gse_metadata=gse_metadata, numstudies=numstudies)

@app.route(f"{ROOT_PATH}/about", methods=['GET', 'POST'])
def about():
    return render_template("about.html", base_path=BASE_PATH, numstudies=numstudies)

@app.route(f"{ROOT_PATH}/mockup", methods=['GET', 'POST'])
def mockup():
    return render_template("mockup.html", base_path=BASE_PATH, numstudies=numstudies)

@app.route(f"{ROOT_PATH}/help", methods=['GET', 'POST'])
def help():
    return render_template("help.html", base_path=BASE_PATH, numstudies=numstudies)

@app.route(f"{ROOT_PATH}/singlegene", methods=['GET', 'POST'])
def singlegene_home():
    return render_template("singlegene.html", base_path=BASE_PATH, numstudies=numstudies)

@app.route(f"{ROOT_PATH}/geneset", methods=['GET', 'POST'])
def geneset_home():
    return render_template("geneset.html", base_path=BASE_PATH, numstudies=numstudies)

@app.route(f'{ROOT_PATH}/scg', methods=['GET', 'POST'])
def scg():
	return render_template('scg.html', base_path=BASE_PATH, numstudies=numstudies)

@app.route(f'{ROOT_PATH}/resources', methods=['GET', 'POST'])
def resources():
	return render_template('resources.html', base_path=BASE_PATH, numstudies=numstudies)

@app.route(f'{ROOT_PATH}/contribute', methods=['GET', 'POST'])
def contribute():
	return render_template('contribute.html', base_path=BASE_PATH, numstudies=numstudies)

@app.route(f'{ROOT_PATH}/downloads', methods=['GET', 'POST'])
def downloads():
	return render_template('downloads.html', base_path=BASE_PATH, gse_metadata=gse_metadata, gse_metadata_single=gse_metadata_single, species_mapping=species_mapping, numstudies=numstudies,  month_dict=month_dict, endpoint=endpoint, version=sigs_version)
####################



########## QUERY APIs ####################
@app.route(f'{ROOT_PATH}/queryexpression', methods=['GET','POST'])
def query_expression():

    gene = request.form['gene']

    result = query_generanger(gene)

    return result

@app.route(f'{ROOT_PATH}/getgwas', methods=['GET','POST'])
def get_gwas():

    gene = request.form['gene']

    if gene == '':
        return {'GWAS_Catalog':[]}

    result = query_gwas(gene)

    return result

@app.route(f'{ROOT_PATH}/getkomp', methods=['GET','POST'])
def get_mgi():
    gene = request.form['gene']
    result = query_mgi(gene)
    return result


@app.route(f'{ROOT_PATH}/getsigcom',  methods=['GET','POST'])
def get_sigcom():
	gene_lists = request.get_json()["genes"]
	if len(gene_lists) == 2:
		res= {'url': sigcom_up_down_genes(gene_lists[0], gene_lists[1])}
	elif len(gene_lists) == 1:
		res= {'url':sigcom_gene_set(gene_lists[0])}
	else: 
		res= {'url': "https://maayanlab.cloud/sigcom-lincs/#/"}
	
	return res

@app.route(f'{ROOT_PATH}/getkea3',  methods=['GET','POST'])
def getkea3():
	geneset = request.form["geneset"]
	geneset = geneset.split(',')
	ADDLIST_URL = 'https://amp.pharm.mssm.edu/kea3/api/enrich/'
	payload = {
        'gene_set': geneset,
        'query_name': ''
    }
	response = requests.post(ADDLIST_URL, data=json.dumps(payload))
	if not response.ok:
		raise Exception('Error analyzing gene list')
	return json.loads(response.text)

@app.route(f'{ROOT_PATH}/getchea3',  methods=['GET','POST'])
def getchea3():
	geneset = request.form["geneset"]
	geneset = geneset.split(',')
	ADDLIST_URL = 'https://maayanlab.cloud/chea3/api/enrich/'
	payload = {
        'gene_set': geneset,
        'query_name': ''
    }
	response = requests.post(ADDLIST_URL, data=json.dumps(payload))
	if not response.ok:
		raise Exception('Error analyzing gene list')

	return json.loads(response.text)

@app.route(f'{ROOT_PATH}/gettfs',  methods=['GET','POST'])
def gettfs():
	gene = request.form['gene']

	if gene.strip() == '':
		return {'data': []}
	result = query_enricher(gene)

	return {'data': result}

@app.route(f'{ROOT_PATH}/getdiabetesenrich',  methods=['GET','POST'])
def getdiabetesenrich():

	genes = request.form['genelist']
	description = request.form['description']

	data = query_enricher_diabetes(genes, description)

	return {'data': data}

@app.route(f'{ROOT_PATH}/enrichrURL',  methods=['GET','POST'])
def enrichrURL():
	ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'

	genes_str = request.form['genelist']
	description = ''
	payload = {
		'list': (None, genes_str),
		'description': (None, description)
	}
	response = requests.post(ENRICHR_URL, files=payload)
	if not response.ok:
		raise Exception('Error analyzing gene list')

	data = json.loads(response.text)
	return {"url": f"https://maayanlab.cloud/Enrichr/enrich?dataset={data['shortId']}"}
	
###########################

##### FILL TABLES ######
@app.route(f'{ROOT_PATH}/getresources',  methods=['GET','POST'])
def resources_api():
	table = get_resources()

	return {'resources': table}

@app.route(f'{ROOT_PATH}/getdownloads',  methods=['GET','POST'])
def downloads_api():
	table = get_downloads()

	return {'downloads': table}


@app.route(f'{ROOT_PATH}/gettweets',  methods=['GET','POST'])
def tweets_api():
	table = get_tweets()
	return {'tweets': table}

@app.route(f'{ROOT_PATH}/getworkflows',  methods=['GET','POST'])
def workflows_api():
	table = get_workflows()

	return {'workflows': table}
###############################


@app.route(f'{ROOT_PATH}/getexample',  methods=['GET','POST'])
def getexample():

	with open('static/searchdata/example_list.txt') as f:
		text = f.read()

	return {'genes': text, 'description': "My gene set"}

@app.route(f'{ROOT_PATH}/getexample2',  methods=['GET','POST'])
def getexample2():

	with open('static/searchdata/example_list2.txt') as f:
		text = f.read()

	return {'genes': text, 'description': "My gene set"}



@app.route(f'{ROOT_PATH}/dgeapi',  methods=['GET','POST'])
def dge():
	response_json = request.get_json()

	perturb = response_json['perturb']
	control = response_json['control']
	method = response_json['method']
	gse = response_json['gse']
	species = response_json['species']
	norms = response_json['norms']
	expr_file = '{base_url}/{species}/{gse}/{gse}_Expression.tsv'.format(species=species, gse=gse, base_url=base_url)
	meta_file = '{base_url}/{species}/{gse}/{gse}_Metadata.tsv'.format(species=species, gse=gse, base_url=base_url)

	data, title = compute_dge(expr_file, meta_file, method, control, perturb, norms['logCPM'], norms['log'], norms['z'], norms['q'])

	jsonplot = make_dge_plot(data,title, method)

	string_data = data.to_string()

	return json.dumps({'table': string_data, 'plot': jsonplot})

@app.route('/dgeapisingle',  methods=['GET','POST'])
def dgesingle():
	response_json = request.get_json()
	method = response_json['method']
	gse = response_json['gse']
	species = response_json['species']
	condition_group = response_json['conditiongroup']
	cluster_group = response_json['diffcluster']

	metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=gse, base_url=base_url),'r')
	metadict = json.load(metajson)
	base_expression_filename = metadict[condition_group]['filename']
	expr_file = '{base_url}/{species}/{gse}/{file}'.format(species=species, gse=gse, base_url=base_url, file=base_expression_filename)
	
	data_dict = compute_dge_single(expr_file, method, 'Cluster', 'leiden',cluster_group, True)

	jsonplot = None
	string_data = None
	description = None
	for key in data_dict:

		jsonplot = make_dge_plot(data_dict[key],key, method)
		string_data = data_dict[key].to_string()
		description = key
		break

	return json.dumps({'table': string_data, 'plot': jsonplot, 'description':description})

@app.route('/api/precomputed_dge',  methods=['GET','POST'])
def fetch_precomputed_dge():
	response_json = request.get_json()
	sig = response_json['sig']
	species = response_json['species']
	dge_tab = get_precomputed_dge(sig, species)
	return dge_tab.to_json(orient='index')


#This function makes the umap tsne and pca plots for the single cell data based off the precomputed coordinates for these plots. 
#It is called in the generate_single_plots within the main.js file. 
@app.route('/singleplots',  methods=['GET','POST'])
def makesingleplots():
	response_json = request.get_json()
	gse = response_json['gse']
	species = response_json['species']
	condition_group = response_json['conditiongroup']
	metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=gse, base_url=base_url),'r')
	metadict = json.load(metajson)
	base_expression_filename = metadict[condition_group]['filename']
	expr_file = '{base_url}/{species}/{gse}/{file}'.format(species=species, gse=gse, base_url=base_url, file=base_expression_filename)
	#The data is originally in genesxcells so transpose to make it cells x genes

	f = read_anndata_h5(expr_file)
	pca_df = pd.DataFrame(data=f['varm/X_pca'][:,:2], columns = ['x', 'y'])
	umap_df = pd.DataFrame(data=f['varm/X_umap'][:], columns = ['x', 'y'])
	tsne_df = pd.DataFrame(data=f['varm/X_tsne'][:,:2], columns = ['x', 'y'])
	clus_numbers = f["var/leiden/codes"][:]
	leiden_values = list(map(lambda x: "Cluster " + str(x), clus_numbers))
	values_dict = {"Cluster": leiden_values}
	category_list_dict = {"Cluster": list(sorted(set(leiden_values)))}
	cells = f['var/column_names'][:].astype(str)
	
	jsonplotumap = make_single_visialization_plot(umap_df, values_dict,'umap', ["Cluster"], cells, "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False)
	jsonplottsne = make_single_visialization_plot(tsne_df, values_dict,'tsne', ["Cluster"], cells, "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False)
	jsonplotpca = make_single_visialization_plot(pca_df, values_dict,'pca', ["Cluster"], cells, "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False)


	return json.dumps({'umapplot': jsonplotumap, 'tsneplot':jsonplottsne, 'pcaplot':jsonplotpca })

#This function gets the different computed leiden clusters from the expression matrix and returns it as json dict for the cluster table on the single viewer page that is called when a new condition-profile is clicked.
@app.route('/getclusterdata', methods=['GET', 'POST'])
def getclusterinfo():
	#The json below holds information about the conditiongroup that we are looking at for this data as well the specific species. 
	response_json = request.get_json()
	gse = response_json['gse']
	species = response_json['species']
	condition_group = response_json['conditiongroup']
	metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=gse, base_url=base_url),'r')
	metadict = json.load(metajson)
	base_expression_filename = metadict[condition_group]['filename']
	expression_file = base_url + '/' + species + '/' + gse + '/' + base_expression_filename
	adata = read_anndata_h5(expression_file)
	#Stores the list of cluster names. 
	leiden_data = adata["var/leiden/categories"][:].astype(str)
	clus_numbers = adata["var/leiden/codes"][:]
	leiden_data_vals = list(map(lambda x: "Cluster " + str(x), clus_numbers))
	classes = sorted(leiden_data)
	classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
	#Stores the number of of cells correlated to each cluster. 
	metadata_dict_counts = pd.Series(leiden_data_vals).value_counts().to_dict()

	return {"classes":classes, "metadict":metadata_dict_counts}


@app.route('/submitcontributionform',  methods=['GET','POST'])
def submit_contribution_form():
	request_json = request.get_json()
	title = request_json['title']
	pmid = request_json['pmid']
	geo = request_json['geo']
	conditions = request_json['conditions']
	model = request_json['model']
	platform = request_json['platform']
	keywords = request_json['keywords']
	authors = request_json['authors']
	email = request_json['email']
	if len(pmid.strip()) == 0:
		pmid = 'N/A'
	if len(geo.strip()) == 0:
		geo = 'N/A'
	if len(email.strip()) == 0:
		email = 'N/A'
	conditions = '|'.join(conditions.strip().split('\n'))
	keywords = ','.join(keywords.strip().split(','))
	authors = '|'.join(authors.strip().split('\n'))
	data_for_sheet = [title,pmid, geo, conditions, model, platform, keywords, authors, email]
	print(data_for_sheet)
	print(request.get_json())
	status = log_suggested_study(data_for_sheet)
	if status == 'success':
		message = "Thank you for suggesting the study: \"{}\"".format(title)
	else:
		message = "There was an error in submitting the data. Please try again".format(title)
	return json.dumps({'response': message, 'status': status})

#############################################
########## 2. Data
#############################################
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
	gse = GEOparse.get_GEO(geo = geo_accession_num, silent=True)

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
	
	metadata_file = f'{base_url}/{species_folder}/{geo_accession}/{geo_accession}_Metadata.tsv'
	
	metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
	gse.metadata['numsamples'] = metadata_dataframe.shape[0]

	return gse.metadata



#### CHECK IF METADATA IS COMPLETE/ IF NEW STUDIES WERE ADDED, ADD THEM TO METADATA

with open(f'static/searchdata/metadata-v2.pickle', 'rb') as f:	
		gse_metadata = pickle.load(f)


url_to_folder = {"human": "human", "mouse": "mouse"}
folder_to_url = {folder:url for url, folder in url_to_folder.items()}

numstudies= [len(gse_metadata['human'].keys()), len(gse_metadata['mouse'].keys())]
study_to_species = {study:species_name for species_name, studies_metadata in gse_metadata.items() for study in studies_metadata.keys()}
species_mapping = {'human': gse_metadata['human'], 'mouse': gse_metadata['mouse']}

with open('static/searchdata/metadatasingle-v2.pickle', 'rb') as f:	
	gse_metadata_single = pickle.load(f)

#For single cell studies
url_to_folder_single = {"human_single": "human_single", "mouse_single": "mouse_single"}
folder_to_url_single = {folder:url for url, folder in url_to_folder_single.items()}
species_mapping_single = {'human_single': gse_metadata_single['human_single'], 'mouse_single': gse_metadata_single['mouse_single']}
#Creating the metadata for the single files only here. Making it separate from the above in case of additions/changes as the project continues. """

#Single cell studies from study to the species name
study_to_species_single = {study:species_name for species_name, studies_metadata in gse_metadata_single.items() for study in studies_metadata.keys()}
numstudies = {'human': len(gse_metadata['human']), 'mouse': len(gse_metadata['mouse']), 'human_single': len(gse_metadata_single['human_single']), 'mouse_single': len(gse_metadata_single['mouse_single'])}
def load_new_studies():
	#Need to use account and pass in order to update the files.
	s3 = s3fs.S3FileSystem(key = os.environ.get('AWS_ACCESS_KEY_ID'), secret = os.environ.get('AWS_SECRET_ACCESS_KEY'))
	#This needs to be the base URL when accessing the s3 bucket through a verified account
	base_url = 'd2h2/data'
	mouse_gses = list(s3.walk(f'{base_url}/mouse', maxdepth=1))[0][1]
	human_gses = list(s3.walk(f'{base_url}/human', maxdepth=1))[0][1]
	#Bulk and microarray study to species name dictionary
	if numstudies['human'] != len(human_gses) or numstudies['mouse'] != len(mouse_gses):
		print("Change in bulk studies")
		for species, geo_accession_ids in species_mapping.items():
			if species not in gse_metadata:
				gse_metadata[species] = {}
			#Adding additional human and mouse studies from the s3 bucket.
			if species =='human':
				for geo_accession in human_gses:
					if geo_accession not in geo_accession_ids:
						gse_metadata[species][geo_accession] = get_metadata(geo_accession, url_to_folder[species])
			if species == 'mouse':
				for geo_accession in mouse_gses:
					if geo_accession not in geo_accession_ids:
						gse_metadata[species][geo_accession] = get_metadata(geo_accession, url_to_folder[species])
		with open('static/searchdata/metadata-v2.pickle', 'wb') as f:
			pickle.dump(gse_metadata, f, protocol=pickle.HIGHEST_PROTOCOL)

	
	mouse_singlegses = list(s3.walk(f'{base_url}/mouse_single', maxdepth=1))[0][1]
	human_singlegses = list(s3.walk(f'{base_url}/human_single', maxdepth=1))[0][1]
	if numstudies['human_single'] != len(human_singlegses) or numstudies['mouse_single'] != len(mouse_singlegses):
		print("Change in single cell studies")
		for species, geo_accession_ids in species_mapping_single.items():
			if species =='human_single':
				for geo_accession in human_singlegses:
					if geo_accession not in geo_accession_ids:
						gse_metadata_single[species][geo_accession] = get_metadata(geo_accession, url_to_folder_single[species])
						metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=geo_accession, base_url=base_url),'r')
						metadict = json.load(metajson)
						dict_to_store_cell_numbers = {}
						total_cells = 0
						for key in metadict:
							#Getting each profile from the single cell data and storing the cell count for each of the different profiles. 
							key_for_dict = key.split(':')[0]
							base_expression_filename = metadict[key]['filename']
							#Add the base path of data here since it uses another modules public s3 credentials 
							expression_file = 'data' + '/' + species + '/' + geo_accession + '/' + base_expression_filename
							adata = read_anndata_h5(expression_file)
							clus_numbers = adata["var/leiden/codes"][:]
							num_cells = len(clus_numbers)
							total_cells += num_cells
							dict_to_store_cell_numbers[key_for_dict] = num_cells
						dict_to_store_cell_numbers['Total'] = total_cells
						gse_metadata_single[species][geo_accession]['cell_count'] = dict_to_store_cell_numbers
						print("Human Single Cell Study")
						print(geo_accession)
						print(dict_to_store_cell_numbers)
			if species == 'mouse_single':
				for geo_accession in mouse_singlegses:
					if geo_accession not in geo_accession_ids:
						gse_metadata_single[species][geo_accession] = get_metadata(geo_accession, url_to_folder_single[species])
						metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=geo_accession, base_url=base_url),'r')
						metadict = json.load(metajson)
						dict_to_store_cell_numbers = {}
						total_cells = 0
						for key in metadict:
							#Getting each profile from the single cell data and storing the cell count for each of the different profiles. 
							key_for_dict = key.split(':')[0]
							base_expression_filename = metadict[key]['filename']
							#Add the base path of data here since it uses another modules public s3 credentials 
							expression_file = 'data' + '/' + species + '/' + geo_accession + '/' + base_expression_filename
							adata = read_anndata_h5(expression_file)
							clus_numbers = adata["var/leiden/codes"][:]
							num_cells = len(clus_numbers)
							total_cells += num_cells
							dict_to_store_cell_numbers[key_for_dict] = num_cells
						dict_to_store_cell_numbers['Total'] = total_cells
						gse_metadata_single[species][geo_accession]['cell_count'] = dict_to_store_cell_numbers
						print("Mouse Single Cell Study")
						print(geo_accession)
						print(dict_to_store_cell_numbers)
		with open('static/searchdata/metadatasingle-v2.pickle', 'wb') as f:
			pickle.dump(gse_metadata_single, f, protocol=pickle.HIGHEST_PROTOCOL)


@app.route(f'{ROOT_PATH}/<studies_or_gse>', methods=['GET', 'POST'])
def species_or_viewerpg(studies_or_gse):
	# test if species
	if studies_or_gse == 'bulkrna':
		num_samples = {'human': sum(map(lambda x: x.get('numsamples'), gse_metadata['human'].values())),'mouse': sum(map(lambda x: x.get('numsamples'), gse_metadata['mouse'].values()))}
		return render_template('species.html', species='human', gse_metadata=gse_metadata, species_mapping=species_mapping, num_samples=num_samples, numstudies=numstudies, month_dict=month_dict)
	#Checking for the single cell studies and loading that summary page
	elif studies_or_gse == 'scrna':
		num_samples = {'human_single': sum(map(lambda x: x.get('numsamples'), gse_metadata_single['human_single'].values())), 'mouse_single': sum(map(lambda x: x.get('numsamples'), gse_metadata_single['mouse_single'].values()))}
		num_cells = {'human_single': sum(map(lambda x: x.get('cell_count').get('Total'), gse_metadata_single['human_single'].values())), 'mouse_single': sum(map(lambda x: x.get('cell_count').get('Total'), gse_metadata_single['mouse_single'].values()))}
		return render_template('single_species.html', species='human_single', gse_metadata_single=gse_metadata_single, species_mapping=species_mapping, gse_metadata=gse_metadata, num_samples=num_samples, numstudies=numstudies,  month_dict=month_dict, num_cells=num_cells)
	
	# test if gse
	elif studies_or_gse in study_to_species:
		geo_accession = studies_or_gse
		species = study_to_species[geo_accession]
		species_folder = url_to_folder[species]
		metadata_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.tsv'
		metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
		gses = list(metadata_dataframe['Sample_geo_accession'])
		metadata_dict = metadata_dataframe.groupby('Group')['Condition'].apply(set).to_dict()
		metadata_dataframe['combined'] = metadata_dataframe['Group'] + ' ' + metadata_dataframe['Condition']
		metadata_dict_samples = metadata_dataframe.groupby('combined')['Sample_geo_accession'].apply(list).to_dict()
		sample_dict = {}
		for key in metadata_dict_samples.keys():
			filtered_samps = list(filter(lambda x: x in gses, metadata_dict_samples[key]))
			sample_dict[key] = {'samples': filtered_samps, 'count': len(filtered_samps)}
		dge_precomputed = get_precomputed_dge_options(geo_accession, species)
		return render_template('viewer.html', metadata_dict=metadata_dict, metadata_dict_samples=sample_dict, geo_accession=geo_accession, gse_metadata=gse_metadata, species=species, species_mapping=species_mapping, numstudies=numstudies,  month_dict=month_dict, dge_precomputed=dge_precomputed)
	
	#Check for the single study individual viewer page
	elif studies_or_gse in study_to_species_single:
		geo_accession = studies_or_gse
		#Will be human_single or mouse_single not just species name
		species = study_to_species_single[geo_accession]
		species_folder = url_to_folder_single[species]
		#The metadata json is in form with key being the condition/profile and the value is a dictionary storing the file name and list of gsms that are part of that condition. 
		meta_file = s3.open(base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_metasep.json', 'r')
		metadata_json = json.load(meta_file)
		list_of_conditions = []
		for key in metadata_json.keys():
			list_of_conditions.append(key)
		default_condition = list_of_conditions[0]
		expression_base_name = metadata_json[default_condition]['filename']
		expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + expression_base_name
		adata = read_anndata_h5(expression_file)
		#Stores the list of cluster names. 
		leiden_data = adata["var/leiden/categories"][:].astype(str)
		clus_numbers = adata["var/leiden/codes"][:]
		leiden_data_vals = list(map(lambda x: "Cluster " + str(x), clus_numbers))
		classes = sorted(leiden_data)
		classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
		#Stores the number of of cells correlated to each cluster. 
		metadata_dict_counts = pd.Series(leiden_data_vals).value_counts().to_dict()
		meta_file.close()
		return render_template('single_viewer.html', study_conditions = list_of_conditions, metadata_dict=classes, metadata_dict_samples=metadata_dict_counts, geo_accession=geo_accession, gse_metadata_single=gse_metadata_single, species=species, species_mapping=species_mapping, gse_metadata=gse_metadata, numstudies=numstudies,  month_dict=month_dict)
	else:
		return render_template('error.html', base_path=BASE_PATH, gse_metadata=gse_metadata, species_mapping=species_mapping, numstudies=numstudies)



#############################################
########## 1. Genes
#############################################

@app.route(f'{ROOT_PATH}/api/genes/<geo_accession>')
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
	elif geo_accession in study_to_species_single:
		species_folder = study_to_species_single[geo_accession]
		expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.h5'
		adata = anndata.read_h5ad(expression_file)
		adata_df = adata.to_df()
		genes_json = json.dumps([{'gene_symbol': x} for x in adata_df.index])
		#go into anndata and get the genes for each human single
	elif geo_accession == 'combined':
		with open('static/searchdata/allgenes-comb.json', 'r') as f:
			human_genes = json.load(f)
		with open('static/searchdata/t2d-mouse.json', 'r') as f:
			mouse_genes = json.load(f)
		all_genes = human_genes + mouse_genes['mouse_genes']	
		genes_json = json.dumps([{'gene_symbol': x} for x in all_genes])
	else:
		species = study_to_species[geo_accession]
		species_folder = url_to_folder[species]

		# Get genes json
		expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.tsv'
		expression_dataframe = pd.read_csv(s3.open(expression_file), index_col = 0, sep='\t')

		genes_json = json.dumps([{'gene_symbol': x} for x in expression_dataframe.index])

	# Return
	return genes_json

#This is the route for the single cell studies with condition as well. It will be called upon initial loading and changing the samples
#to focus on
@app.route('/api/singlegenes/<geo_accession>/<condition>')

# this will likely stay the same.
@lru_cache(maxsize=None)
def genes_api_single(geo_accession, condition):
	species_folder = study_to_species_single[geo_accession]
	meta_file = s3.open(base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_metasep.json', 'r')
	metadata_json = json.load(meta_file)
	expression_base_name = metadata_json[condition]['filename']
	expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + expression_base_name
	f = read_anndata_h5(expression_file)
	genes = np.array(f['obs/gene_symbols'][:].astype(str))
	genes_json = json.dumps([{'gene_symbol': x} for x in genes])
	#go into anndata and get the genes for each human single

	return genes_json
#############################################
########## 2. Plot
#############################################
#SINGLE CELL DATA PLOT
@app.route('/api/plot_single/<geo_accession>/<condition>', methods=['GET', 'POST'])

def plot_api_single(geo_accession, condition):
	species = study_to_species_single[geo_accession]
	species_folder = url_to_folder_single[species]
	assay = gse_metadata_single[species][geo_accession].get('type')[0]
	meta_file = s3.open(base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_metasep.json', 'r')
	metadata_json = json.load(meta_file)
	base_expression_name = metadata_json[condition]['filename']
	expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + base_expression_name
	#expression_adata = read_anndata_raw(expression_file)
	data = request.json
	gene_symbol = data['gene_symbol']

	conditions = data['conditions']

	f = read_anndata_h5(expression_file)

	genes = np.array(f['obs/gene_symbols'][:].astype(str))

	clus_numbers = f["var/leiden/codes"][:]

	leiden_values = list(map(lambda x: "Cluster " + str(x), clus_numbers))

	idx = np.where(genes == gene_symbol)[0][0]

	vals = f['raw/X'][idx,:]

	# Create a new dataframe that maps each sample to its condition
	melted_dataframe = pd.DataFrame(data=[vals, leiden_values]).T
	melted_dataframe.columns = ['expr_vals', 'Condition']
	
	# Get plot dataframe
	plot_dataframe = melted_dataframe.groupby('Condition')['expr_vals'].agg([np.mean, np.std, lambda x: list(x)])#.rename(columns={'<lambda>': 'points'})#.reindex(conditions)
	plot_dataframe = plot_dataframe.rename(columns={plot_dataframe.columns[-1]: 'points'})

	# Initialize figure
	fig = go.Figure()
	
	# Loop
	for condition in conditions:
		if len(plot_dataframe.loc[condition, 'points']) > 1:
			fig.add_trace(go.Box(name=condition, y=plot_dataframe.loc[condition, 'points'], boxpoints='all', pointpos=0))
		else:
			fig.add_trace(go.Scatter(name=condition, x=[condition], y=plot_dataframe.loc[condition, 'points']))

	y_units = 'Expression'

	fig.update_layout(
		title = {'text': gene_symbol+' gene expression', 'x': 0.5, 'y': 0.85, 'xanchor': 'center', 'yanchor': 'top'},
		xaxis_title = 'Condition',
		yaxis_title = y_units,
		showlegend = False
	)
	return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

@app.route(f'{ROOT_PATH}/api/plot/<geo_accession>', methods=['GET', 'POST'])
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

	expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.tsv'
	metadata_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.tsv'
	expression_dataframe = pd.read_csv(s3.open(expression_file), index_col = 0, sep='\t')
	metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
	
	# Get data

	data = request.json
	gene_symbol = data['gene_symbol']

	conditions = data['conditions']



	# Create a new dataframe that maps each sample to its condition

	melted_dataframe = expression_dataframe.loc[gene_symbol].rename('expr_vals').rename_axis('Sample_geo_accession').reset_index().merge(metadata_dataframe, on="Sample_geo_accession")

	# Get plot dataframe
	melted_dataframe['Combined'] = melted_dataframe['Condition'] + ' ' + melted_dataframe['Group']

	plot_dataframe = melted_dataframe.groupby('Combined')['expr_vals'].agg([np.mean, np.std, lambda x: list(x)])#.rename(columns={'<lambda>': 'points'})#.reindex(conditions)
	plot_dataframe = plot_dataframe.rename(columns={plot_dataframe.columns[-1]: 'points'})
	

	# Initialize figure
	fig = go.Figure()

	condition_name = False
	if len(melted_dataframe['Group'].unique()) == 1:
		condition_name = True
		group = melted_dataframe['Group'].values[0]
	for condition in conditions:
		if condition_name:
			displayname = condition.replace(group, '').strip()
		else:
			displayname = condition
		if len(plot_dataframe.loc[condition, 'points']) > 1:
			fig.add_trace(go.Box(name=displayname, y=plot_dataframe.loc[condition, 'points'], boxpoints='all', pointpos=0))
		else:
			fig.add_trace(go.Scatter(name=displayname, x=[condition], y=plot_dataframe.loc[condition, 'points']))
	
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
		showlegend = False,
		xaxis = dict(
		tickangle = 30,
		tickfont = dict(size=10)
    	)
	)
	return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

@app.route(f'{ROOT_PATH}/api/volcano', methods=['GET', 'POST'])
def plot_volcano_api():
	request.form
	gene = request.form["gene"]
	species = request.form["species"]
	json_item_plot = send_plot(species, gene)
	# Return
	return json.dumps(json_item_plot)


########## 3. Conditions ##########

@app.route(f'{ROOT_PATH}/api/conditions/<geo_accession>')
def conditions_api(geo_accession):
	species = study_to_species[geo_accession]
	species_folder = url_to_folder[species]
	metadata_file = base_url+ '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.tsv'
	metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
	if len(set(metadata_dataframe['Group'])) == 1:
		conditions = set(metadata_dataframe['Condition'])
		col = 'Condition'
	else:
		metadata_dataframe['combined'] = metadata_dataframe['Condition'] + ' ' + metadata_dataframe['Group']
		conditions = set(metadata_dataframe['combined'])
		col = 'combined'
	valid_conditions = []
	for condition in conditions:

		if len(metadata_dataframe[metadata_dataframe[col] == condition]) > 1:
			valid_conditions.append(condition)

	return json.dumps([{'Condition': x} for x in valid_conditions])

#############################################
########## Conditions to Samples
#############################################
@app.route(f'{ROOT_PATH}/api/samples/<geo_accession>')

def samples_api(geo_accession):
	species = study_to_species[geo_accession]
	species_folder = url_to_folder[species]

	metadata_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.tsv'
	metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')

	conditions_mapping = metadata_dataframe.groupby('Condition')['Sample_geo_accession'].apply(list).to_dict()

	return json.dumps(conditions_mapping)


@app.route(f'{ROOT_PATH}/api/data',  methods=['GET', 'POST'])
def get_study_data():

	response_json = request.get_json()
	geo_accession = response_json['gse']
	control = response_json['control']
	perturb = response_json['perturb']
	species = response_json['species']

	metadata_file = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Metadata.tsv'
	expression_file = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Expression.tsv'

	with s3.open(metadata_file, 'r') as f:
		meta_data = f.read()

	with s3.open(expression_file, 'r') as f:
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
		elif f'{entries[1]} {entries[2]}' == control or f'{entries[1]} {entries[2]}' == perturb:
			entries[1] = f'{entries[1]} {entries[2]}'
			selected_conditions.append('\t'.join(entries))

	selected_meta = '\n'.join(selected_conditions)

	data_dict = {'meta': selected_meta, 'expression': expression_data}

	return data_dict
	

@app.route(f'{ROOT_PATH}/api/bulksampvis',  methods=['GET', 'POST'])
def visualize_samps():
	response_json = request.get_json()
	geo_accession = response_json['gse']
	species = response_json['species']
	meta_df = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Metadata.tsv'
	
	meta_df = pd.read_csv(s3.open(meta_df), sep='\t', index_col=0)

	if len(meta_df.columns) > 2:
		pca_plot = interactive_circle_plot(meta_df.rename(columns={'pca_x': 'x', 'pca_y': 'y'}), "PC-1", "PC-2", 'Condition', 'pca')
		tsne_plot = interactive_circle_plot(meta_df.rename(columns={'tsne_x': 'x', 'tsne_y': 'y'}), "t-SNE-1", "t-SNE-2", 'Condition', 'tsne')
		umap_plot = interactive_circle_plot(meta_df.rename(columns={'umap_x': 'x', 'umap_y': 'y'}), "UMAP-1", "UMAP-2", 'Condition', 'umap')
	else:
		pca_plot, tsne_plot, umap_plot = 'None', 'None', 'None'

	return json.dumps({'pcaplot': pca_plot, 'tsneplot': tsne_plot, 'umapplot': umap_plot})

@app.route(f'{ROOT_PATH}/api/query_gpt',  methods=['GET', 'POST'])
def query_gpt():
	response_json = request.get_json()
	query = response_json['query']
	res = find_process(query)

	return res

@app.route(f'{ROOT_PATH}/api/query_options',  methods=['GET', 'POST'])
def query_options():
	response_json = request.get_json()
	response = response_json['response']
	options = response_json['options']
	res = select_option(response, options)
	return res

@app.route(f'{ROOT_PATH}/api/query_genes',  methods=['GET', 'POST'])
def query_genes():
	response_json = request.get_json()
	g = response_json['gene']
	res = infer_gene(g)
	return res

@app.route(f'{ROOT_PATH}/api/record_chat',  methods=['GET', 'POST'])
def record_chat():
	try:
		response_json = request.get_json()
		user_chat = response_json['user_chat']
		response = response_json['response']
		userid = response_json['user_id']
	except:
		print('Error logging chat:')
		print('User chat:', 'user_chat' in response_json)
		print('D2H2 response:', 'response' in response_json)
		print('User id:', 'user_id' in response_json)
	if not DEBUG:
		log_chat(user_chat, response, userid)
	return {}

@app.route('/api/run_geneshot', methods=['GET', 'POST'])
def run_geneshot():
	if request.method == "POST":
		term = request.get_json()['term']
		try:
			res = query_geneshot(term) 
			return {'search_term':term, 'genes': '\t'.join(res[0]), 'count': str(res[1])}
		except:
			return {'search_term':term, 'error': True}


	
#######################################################
#######################################################
########## 3. Run App
#######################################################
#######################################################
if __name__ == "__main__":
	if UPDATE_STUDIES:
		load_new_studies()
	
	app.run(debug=DEBUG, host="0.0.0.0", port=5000)

