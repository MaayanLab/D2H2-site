from flask import Flask, render_template, request
from waitress import serve
import os
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
import datetime
from helpers import *
from twitterauth import update_tweets_table
from dge import *
import anndata
import scanpy as sc
from sklearn.preprocessing import StandardScaler

endpoint = os.environ.get('ENDPOINT', 'https://minio.dev.maayanlab.cloud/')
base_url = os.environ.get('BASE_URL', 'd2h2/data')
ROOT_PATH = os.environ.get('ROOT_PATH', '/')
BASE_PATH = os.environ.get('BASE_PATH', 'maayanlab.cloud')

print(endpoint)
s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})

print(s3.ls(''))
app = Flask(__name__, static_url_path=ROOT_PATH + 'static')


@app.route(ROOT_PATH, methods=['GET', 'POST'])
def home():
	#update_tweets_table(datetime.datetime.date)
	return render_template('home.html', base_path=BASE_PATH, gse_metadata=gse_metadata, numstudies=numstudies)

@app.route(f"{ROOT_PATH}/about", methods=['GET', 'POST'])
def about():
    return render_template("about.html", base_path=BASE_PATH, numstudies=numstudies)

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

@app.route(f'{ROOT_PATH}/downloads', methods=['GET', 'POST'])
def downloads():
	return render_template('downloads.html', base_path=BASE_PATH, numstudies=numstudies)

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

@app.route(f'{ROOT_PATH}/gettfs',  methods=['GET','POST'])
def gettfs():
	gene = request.form['gene']

	if gene.strip() == '':
		return {'data': []}
	result = query_enricher(gene)

	return {'data': result}

@app.route(f'{ROOT_PATH}/getexample',  methods=['GET','POST'])
def getexample():

	with open('static/searchdata/example_list.txt') as f:
		text = f.read()

	return {'genes': text, 'description': "GSE136134 Ctrl-vs-Insulin 24hrs Human BulkRNAseq hiPSCs_down"}

@app.route(f'{ROOT_PATH}/getdiabetesenrich',  methods=['GET','POST'])
def getdiabetesenrich():

	genes = request.form['genelist']
	description = request.form['description']

	data = query_enricher_diabetes(genes, description)

	return {'data': data}


@app.route(f'{ROOT_PATH}/dgeapi',  methods=['GET','POST'])
def dge():
	response_json = request.get_json()

	perturb = response_json['perturb']
	control = response_json['control']
	method = response_json['method']
	gse = response_json['gse']
	species = response_json['species']
	norms = response_json['norms']
	expr_file = '{base_url}/{species}/{gse}/{gse}_Expression.txt'.format(species=species, gse=gse, base_url=base_url)
	meta_file = '{base_url}/{species}/{gse}/{gse}_Metadata.txt'.format(species=species, gse=gse, base_url=base_url)
	if method == 'limma' or method == 'edgeR':
		data, title = compute_dge(expr_file, meta_file, method, control, perturb, False, False, False, False)
	else:
		data, title = compute_dge(expr_file, meta_file, method, control, perturb, norms['logCPM'], norms['log'], norms['z'], norms['q'])

	jsonplot = make_dge_plot(data,title, method)

	string_data = data.to_string()

	return json.dumps({'table': string_data, 'plot': jsonplot})
@app.route('/dgeapisingle',  methods=['GET','POST'])
def dgesingle():
	response_json = request.get_json()
	print(response_json)
	method = response_json['method']
	gse = response_json['gse']
	species = response_json['species']
	condition_group = response_json['conditiongroup']
	cluster_group = response_json['diffcluster']
	norms = response_json['norms']
	metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=gse, base_url=base_url),'r')
	metadict = json.load(metajson)
	base_expression_filename = metadict[condition_group]['filename']
	expr_file = '{base_url}/{species}/{gse}/{file}'.format(species=species, gse=gse, base_url=base_url, file=base_expression_filename)
	#compute_dge_single(adata, diff_gex_method, 'Cluster', 'leiden', True)
	#adata = read_anndata_raw(expr_file)
	#adata_raw = adata.raw.to_adata()
	#Passing the gene information as cells x genes for the differential expression method. 
	#adata = adata.T
	#adata.raw = adata_raw.T
	
	#if method == 'limma' or method == 'edgeR':
	data_dict = compute_dge_single(expr_file, method, 'Cluster', 'leiden',cluster_group, True)
	#data_dict = compute_dge_single(expr_file, method, 'Cluster', 'leiden',cluster_group, True)

	jsonplot = None
	string_data = None
	description = None
	for key in data_dict:

		jsonplot = make_dge_plot(data_dict[key],key, method)
		string_data = data_dict[key].to_string()
		description = key
		break

	return json.dumps({'table': string_data, 'plot': jsonplot, 'description':description})


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



	""" time_start = datetime.datetime.now()
	adata = read_anndata_raw(expr_file).T
	time_end = datetime.datetime.now()
	print('TIME TO READ ANNDATA ' + str(time_end-time_start))
	values_dict = dict()
	values_dict["Cluster"] = adata.obs["leiden"].values
	category_list_dict = dict()
	category_list_dict["Cluster"] = list(sorted(adata.obs["leiden"].unique()))
	#Each of the values for the umap, tsne, and pca were precomputed so obtain them from the anndata object. 
	umap_df = pd.DataFrame(adata.obsm['X_umap'])
	umap_df.columns = ['x', 'y']
	# scaler = StandardScaler().fit(adata.obsm['X_pca'][:,:2])
	# X_scaled = scaler.transform(adata.obsm['X_pca'][:,:2])
	pca_df = pd.DataFrame(adata.obsm['X_pca'][:,:2])
	pca_df.columns = ['x', 'y']

	print(pca_df)
	tsne_df = pd.DataFrame(adata.obsm['X_tsne'][:,:2])
	tsne_df.columns = ['x', 'y'] """
	
	jsonplotumap = make_single_visialization_plot(umap_df, values_dict,'umap', ["Cluster"], cells, "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False)
	jsonplottsne = make_single_visialization_plot(tsne_df, values_dict,'tsne', ["Cluster"], cells, "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False)
	jsonplotpca = make_single_visialization_plot(pca_df, values_dict,'pca', ["Cluster"], cells, "Scatter plot of the samples. Each dot represents a sample and it is colored by ", category_list_dict=category_list_dict, category=True, dropdown=False)


	return json.dumps({'umapplot': jsonplotumap, 'tsneplot':jsonplottsne, 'pcaplot':jsonplotpca })

#This function gets the different computed leiden clusters from the expression matrix and returns it as json dict for the cluster table on the single viewer page that is called when a new condition-profile is clicked.
@app.route('/getclusterdata', methods=['GET', 'POST'])
def getclusterinfo():
	print("IN GET CLUSTER DATA")
	#The json below holds information about the conditiongroup that we are looking at for this data as well the specific species. 
	response_json = request.get_json()
	print(response_json)
	gse = response_json['gse']
	species = response_json['species']
	condition_group = response_json['conditiongroup']
	metajson = s3.open('{base_url}/{species}/{gse}/{gse}_metasep.json'.format(species=species, gse=gse, base_url=base_url),'r')
	metadict = json.load(metajson)
	base_expression_filename = metadict[condition_group]['filename']
	expression_file = base_url + '/' + species + '/' + gse + '/' + base_expression_filename
	#Transposing to get the data with cells as rows and genes as columns. 

	adata = read_anndata_h5(expression_file)
	print('after reading adata')
	#Stores the list of cluster names. 
	leiden_data = adata["var/leiden/categories"][:].astype(str)
	clus_numbers = adata["var/leiden/codes"][:]
	leiden_data_vals = list(map(lambda x: "Cluster " + str(x), clus_numbers))
	classes = sorted(leiden_data)
	classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
	#Stores the number of of cells correlated to each cluster. 
	metadata_dict_counts = pd.Series(leiden_data_vals).value_counts().to_dict()

	""" adata = read_anndata_raw(expression_file).T
	classes = sorted(adata.obs["leiden"].unique().tolist())
	classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
	metadata_dict_counts = adata.obs["leiden"].value_counts().to_dict() """


	return {"classes":classes, "metadict":metadata_dict_counts}
	
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
	
	
	metadata_file = f'{base_url}/{species_folder}/{geo_accession}/{geo_accession}_Metadata.txt'
	#Why was it -1
	
	metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
	gse.metadata['numsamples'] = metadata_dataframe.shape[0]

	return gse.metadata



#### CHECK IF METADATA IS COMPLETE/ IF NEW STUDIES WERE ADDED, ADD THEM TO METADATA

with open(f'static/searchdata/metadata-v1.pickle', 'rb') as f:	
		gse_metadata = pickle.load(f)


url_to_folder = {"human": "human", "mouse": "mouse"}
folder_to_url = {folder:url for url, folder in url_to_folder.items()}

""" numstudies= [len(gse_metadata['human'].keys()), len(gse_metadata['mouse'].keys())]
mouse_gses = list(s3.walk(f'{base_url}/mouse', maxdepth=1))[0][1]
human_gses = list(s3.walk(f'{base_url}/human', maxdepth=1))[0][1]






#Bulk and microarray study to species name dictionary

if numstudies[0] == len(human_gses) and numstudies[1] == len(mouse_gses):
	for species, geo_accession_ids in species_mapping.items():
		if species not in gse_metadata:
			gse_metadata[species] = {}
		for geo_accession in geo_accession_ids:
			if geo_accession not in gse_metadata[species]:
				gse_metadata[species][geo_accession] = get_metadata(geo_accession, url_to_folder[species])
	with open('static/searchdata/metadata-v1.pickle', 'wb') as f:
		pickle.dump(gse_metadata, f, protocol=pickle.HIGHEST_PROTOCOL) """
species_mapping = {'human': gse_metadata['human'], 'mouse': gse_metadata['mouse']}

study_to_species = {study:species_name for species_name, studies_metadata in gse_metadata.items() for study in studies_metadata.keys()}

######### SINGLE CELL METADATA CREATION BASED OFF ABOVE FOR BULK AND MICROARRAY#####
#mouse_singlegses = list(s3.walk(f'{base_url}/mouse_single', maxdepth=1))[0][1]
#human_singlegses = list(s3.walk(f'{base_url}/human_single', maxdepth=1))[0][1]

with open('static/searchdata/metadatasingle-v1.pickle', 'rb') as f:	
	gse_metadata_single = pickle.load(f)

#For single cell studies
url_to_folder_single = {"human_single": "human_single", "mouse_single": "mouse_single"}
folder_to_url_single = {folder:url for url, folder in url_to_folder_single.items()}
species_mapping_single = {'human_single': gse_metadata_single['human_single'], 'mouse_single': gse_metadata_single['mouse_single']}
#Creating the metadata for the single files only here. Making it separate from the above in case of additions/changes as the project continues. """




""" numstudies_single= [len(gse_metadata_single['human_single'].keys()), len(gse_metadata_single['mouse_single'].keys())]

if numstudies_single[0] != len(human_singlegses) and numstudies_single[1] != len(mouse_singlegses):
	for species, geo_accession_ids in species_mapping_single.items():
		for geo_accession in geo_accession_ids:
			if geo_accession not in gse_metadata_single[species]:
				gse_metadata_single[species][geo_accession] = get_metadata(geo_accession, url_to_folder_single[species])
	with open('static/searchdata/metadatasingle-v1.pickle', 'wb') as f:
		pickle.dump(gse_metadata_single, f, protocol=pickle.HIGHEST_PROTOCOL)
"""
#Single cell studies from study to the species name
study_to_species_single = {study:species_name for species_name, studies_metadata in gse_metadata_single.items() for study in studies_metadata.keys()}


numstudies = [len(gse_metadata['human']), len(gse_metadata['mouse']), len(gse_metadata_single['human_single']), len(gse_metadata_single['mouse_single'])]


@app.route(f'{ROOT_PATH}/<species_or_gse>', methods=['GET', 'POST'])
def species_or_viewerpg(species_or_gse):
	# test if species
	if species_or_gse in gse_metadata:
		num_samples = sum(map(lambda x: x.get('numsamples'), gse_metadata[species_or_gse].values()))
		return render_template('species.html', species=species_or_gse, gse_metadata=gse_metadata, species_mapping=species_mapping, num_samples=num_samples, numstudies=numstudies)
	#Checking for the single cell studies and loading that summary page
	elif species_or_gse in gse_metadata_single:
		num_samples = sum(map(lambda x: x.get('numsamples'), gse_metadata_single[species_or_gse].values()))
		return render_template('single_species.html', species=species_or_gse, gse_metadata_single=gse_metadata_single, species_mapping=species_mapping, gse_metadata=gse_metadata, num_samples=num_samples, numstudies=numstudies)
	# test if gsea
	elif species_or_gse in study_to_species:
		geo_accession = species_or_gse
		species = study_to_species[geo_accession]
		species_folder = url_to_folder[species]
		metadata_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
		metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
		metadata_dict = metadata_dataframe.groupby('Group')['Condition'].apply(set).to_dict()
		metadata_dict_samples = metadata_dataframe.groupby('Condition')['Sample_geo_accession'].apply(list).to_dict()
		sample_dict = {}
		for key in metadata_dict_samples.keys():
			sample_dict[key] = {'samples':metadata_dict_samples[key], 'count': len(metadata_dict_samples[key])}
		return render_template('viewer.html', metadata_dict=metadata_dict, metadata_dict_samples=sample_dict, geo_accession=geo_accession, gse_metadata=gse_metadata, species=species, species_mapping=species_mapping, numstudies=numstudies)
	#Check for the single study individual viewer page
	elif species_or_gse in study_to_species_single:
		geo_accession = species_or_gse
		#Will be human_single or mouse_single not just species name
		species = study_to_species_single[geo_accession]
		species_folder = url_to_folder_single[species]
		#The metadata json is in form with key being the condition/profile and the value is a dictionary storing the file name and list of gsms that are part of that condition. 
		print('before reading json file')
		meta_file = s3.open(base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_metasep.json', 'r')
		metadata_json = json.load(meta_file)
		print(metadata_json)
		list_of_conditions = []
		for key in metadata_json.keys():
			list_of_conditions.append(key)
		default_condition = list_of_conditions[0]
		expression_base_name = metadata_json[default_condition]['filename']
		expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + expression_base_name
		print(expression_file)
		# expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.h5'
		print('before_reading anndata')
		# adata = anndata.read_h5ad(s3.open(expression_file)).T
		# print('after_reading_anndata')
		# classes = sorted(adata.obs["leiden"].unique().tolist())
		# #Stores the list of cluster names. 
		# classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
		# #Stores the number of of cells correlated to each cluster. 
		# metadata_dict_counts = adata.obs["leiden"].value_counts().to_dict()
		adata = read_anndata_h5(expression_file)
		print('after reading adata')
		#Stores the list of cluster names. 
		leiden_data = adata["var/leiden/categories"][:].astype(str)
		clus_numbers = adata["var/leiden/codes"][:]
		leiden_data_vals = list(map(lambda x: "Cluster " + str(x), clus_numbers))
		classes = sorted(leiden_data)
		classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
		#Stores the number of of cells correlated to each cluster. 
		metadata_dict_counts = pd.Series(leiden_data_vals).value_counts().to_dict()
		meta_file.close()
		return render_template('single_viewer.html', study_conditions = list_of_conditions, metadata_dict=classes, metadata_dict_samples=metadata_dict_counts, geo_accession=geo_accession, gse_metadata_single=gse_metadata_single, species=species, species_mapping=species_mapping, gse_metadata=gse_metadata, numstudies=numstudies)
	else:
		return render_template('error.html', base_path=BASE_PATH, gse_metadata=gse_metadata, species_mapping=species_mapping, numstudies=numstudies)



#############################################
########## 1. Genes
#############################################


# this will likely stay the same.

@app.route(f'{ROOT_PATH}/api/genes/<geo_accession>')

# this will likely stay the same.
@lru_cache(maxsize=None)
def genes_api(geo_accession):
	print(geo_accession)
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
		print([{'gene_symbol': x} for x in adata_df.index][:10])
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
		expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'
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
	print("IN SINGLE API GENES")
	print(geo_accession)
	print(condition)
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
	"""
	Inputs:
	- expression_dataframe, a tsv file representing a matrix with rows having indices representing gene symbols and columns with indices representing Sample ID's. 
		Each entry represents the expression value for a gene in a sample.
	- metadata_dict, a JSON file with the following format: {group_name:{condition_names:[sample_name, ...]}}
	Outputs: 
	- plotly_json, a serialized JSON formatted string that represents the data for the boxplot to be plotted, used by boxplot() function in scripts.js.
	"""
	species = study_to_species_single[geo_accession]
	species_folder = url_to_folder_single[species]
	time_start = datetime.datetime.now()
	print("IN PLOT AP SINGLE")
	# print(condition)
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


	melted_dataframe = pd.DataFrame(data=[vals, leiden_values]).T
	melted_dataframe.columns = ['expr_vals', 'Condition']
	time_now = datetime.datetime.now()

	""" expression_adata = read_anndata_raw(expression_file)
	
	expression_adata = expression_adata.raw.to_adata()
	expression_dataframe = expression_adata.to_df()
	adata = expression_adata.T
	
	print(time_now-time_start)
	leiden_vals = adata.obs["leiden"].tolist()
	cell_names = adata.obs['column_names'].tolist()
	print(leiden_vals)
	print(cell_names)

	df_dict = {'Sample_geo_accession':cell_names, 'Condition':leiden_vals}
	metadata_dataframe = pd.DataFrame.from_dict(df_dict) """
	# Get data
	
	# print(conditions)

	# Create a new dataframe that maps each sample to its condition
	#melted_dataframe = expression_dataframe.loc[gene_symbol].rename('expr_vals').rename_axis('Sample_geo_accession').reset_index().merge(metadata_dataframe, on="Sample_geo_accession")
	
	print(melted_dataframe)
	# Get plot dataframe
	plot_dataframe = melted_dataframe.groupby('Condition')['expr_vals'].agg([np.mean, np.std, lambda x: list(x)])#.rename(columns={'<lambda>': 'points'})#.reindex(conditions)
	print(plot_dataframe.shape)
	plot_dataframe = plot_dataframe.rename(columns={plot_dataframe.columns[-1]: 'points'})
	time_after = datetime.datetime.now()
	print(time_after - time_now)
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

	expression_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'
	metadata_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
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

	# print(plot_dataframe)
	
	# Loop
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
	metadata_file = base_url+ '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	metadata_dataframe = pd.read_csv(s3.open(metadata_file), sep='\t')
	conditions = set(metadata_dataframe['Condition'])
	return json.dumps([{'Condition': x} for x in conditions])

#############################################
########## Conditions to Samples
#############################################
@app.route(f'{ROOT_PATH}/api/samples/<geo_accession>')

def samples_api(geo_accession):
	species = study_to_species[geo_accession]
	species_folder = url_to_folder[species]

	metadata_file = base_url + '/' + species_folder + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
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

	metadata_file = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	expression_file = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'

	

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


	selected_meta = '\n'.join(selected_conditions)

	data_dict = {'meta': selected_meta, 'expression': expression_data}

	return data_dict
	

@app.route(f'{ROOT_PATH}/api/bulksampvis',  methods=['GET', 'POST'])
def visualize_samps():
	response_json = request.get_json()
	geo_accession = response_json['gse']
	species = response_json['species']
	meta_df = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Metadata.txt'
	expr_df = base_url + '/' + species + '/' + geo_accession + '/' + geo_accession + '_Expression.txt'

	pca_df, tsne_df, umap_df = bulk_vis(expr_df, meta_df)

	pca_plot = interactive_circle_plot(pca_df, "PC-1", "PC-2", pca_df.columns[0], 'pca')

	tsne_plot = interactive_circle_plot(tsne_df, "t-SNE-1", "t-SNE-2", tsne_df.columns[0], 'tsne')

	umap_plot = interactive_circle_plot(umap_df, "UMAP-1", "UMAP-2", umap_df.columns[0], 'umap')


	return json.dumps({'pcaplot': pca_plot, 'tsneplot': tsne_plot, 'umapplot': umap_plot})


#######################################################
#######################################################
########## 3. Run App
#######################################################
#######################################################
if __name__ == "__main__":

	#serve(app, host="0.0.0.0", port=5000)
	app.run(debug=True, host="0.0.0.0")


