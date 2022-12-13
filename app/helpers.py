from functools import lru_cache
import json
import requests
from intermine.webservice import Service
import pandas as pd
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize
# bokeh
from bokeh.plotting import figure
from bokeh.embed import json_item
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool, CustomJS, Span, Select, Legend, PreText, Paragraph, LinearColorMapper, ColorBar, CategoricalColorMapper
from bokeh.palettes import Category20
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import anndata
import s3fs
import h5py
import os
import re
import hashlib


base_url = os.environ.get('BASE_URL')

s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': 'https://minio.dev.maayanlab.cloud/'})

########################## QUERY ENRICHER ###############################

@lru_cache()
def query_enricher(gene):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/genemap'
    query_string = '?json=true&setup=true&gene=%s'



    response = requests.get(ENRICHR_URL + query_string % gene)
    if not response.ok:
        raise Exception('Error finding gene')

    data = json.loads(response.text)



    for l in data['categories']:
        if l['name'] == 'Transcription':
            transcription_libs = l['libraries']

    res = []
    for lib in transcription_libs:
        if lib['name'] in data['gene']:
            lib['tfs'] = data["gene"][lib['name']]
            res.append(lib)
    return res


@lru_cache()
def query_enricher_diabetes(genelist, description):
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = genelist

    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    listid = data["userListId"]


    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = listid
    gene_set_library = 'Diabetes_Perturbations_GEO_2022'
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
    )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)


    return data


########################## QUERY KOMP API ###############################


@lru_cache()
def query_komp(gene: str):

    gene =  gene[0].upper() + (gene[1:]).lower()
    KOMP_URL = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select?q=marker_symbol:" + gene
    response = requests.get(KOMP_URL)
    if not response.ok:
        raise Exception('Error analyzing retrieving information')
    data = json.loads(response.text)
    return data



######## QUERY MGI ##############

@lru_cache()
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

@lru_cache()
def query_gwas(gene: str):
    gene = gene.upper()
    url = "https://www.ebi.ac.uk/gwas/api/search/downloads?q=ensemblMappedGenes:" + gene + "&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&genotypingfilter[]=&traitfilter[]=&dateaddedfilter=&facet=association&efo=true"

    df = pd.read_csv(url, sep='\t')

    summerized_data = pd.DataFrame(columns=['gene', 'trait', 'mapped_trait_link', 'count'])
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



############## GET RESOURCES TABLE FROM GOOGLE DRIVE ##############


@lru_cache()
def get_resources():

    table = pd.read_csv('./static/searchdata/resources.csv', index_col=0)
    table.fillna('', inplace=True)
    resources_list = table.values.tolist()
    table_list = [list(table.columns.values)] + resources_list
    return table_list
    

#### UPDATE RESOURCES TABALE ####


def update_resources():
    sheet_url = "https://docs.google.com/spreadsheets/d/1uKuDw8eN3si7QmukL7n4dQsOdn3qN8xUiu-rRapC1Cs/edit#gid=0"
    url_1 = sheet_url.replace('/edit#gid=', '/export?format=csv&gid=')
    table = pd.read_csv(url_1)
    table.to_csv('static/searchdata/resources.csv')

#update_resources()

@lru_cache()
def get_downloads():

    table = pd.read_csv('./static/searchdata/downloads.csv', index_col=0)
    table.fillna('', inplace=True)
    resources_list = table.values.tolist()
    table_list = [list(table.columns.values)] + resources_list
    return table_list
    

#### UPDATE RESOURCES TABALE ####

def update_downloads():
    sheet_url = "https://docs.google.com/spreadsheets/d/17if2nhNAOQMESA6EyO_eicqinMmT-P5Ly3TxbkSZZRQ/edit#gid=0"
    url_1 = sheet_url.replace('/edit#gid=', '/export?format=csv&gid=')
    table = pd.read_csv(url_1)
    table.to_csv('static/searchdata/downloads.csv')


#update_downloads()

def get_workflows():

    table = pd.read_csv('./static/searchdata/workflows.csv')
    resources_list = table.values.tolist()
    table_list = [list(table.columns.values)] + resources_list
    return table_list

def get_tweets():

    table = pd.read_csv('./static/searchdata/tweets.csv', index_col=None)
    resources_list = table.values.tolist()
    table_list = [list(table.columns.values)] + resources_list
    return table_list



red_map = cm.get_cmap('Reds_r')
red_norm = colors.Normalize(vmin=-0.25, vmax=1)
blue_map = cm.get_cmap('Blues_r')
blue_norm = colors.Normalize(vmin=-0.25, vmax=1)

def load_files(species, gene):
    root_path = 'static/searchdata/'
    pval_rna_df = pd.read_feather(f'{root_path}all_{species}_pval.f', columns=['index', gene]).set_index('index')
    # RNA-seq fold change
    fc_rna_df = pd.read_feather(f'{root_path}all_{species}_fc.f', columns=['index', gene]).set_index('index')

    # microarray data
    try:
        has_micro = True 
        pval_micro_df = pd.read_feather(f"{root_path}{species}_pruned_affy_pv.f", columns=['index', gene]).set_index('index')
        fc_micro_df = pd.read_feather(f"{root_path}{species}_pruned_affy_fc.f", columns=['index', gene]).set_index('index')
    except:
        print('no micro data')
        has_micro = False 
        pval_micro_df = pd.DataFrame()
        fc_micro_df = pd.DataFrame()
    inst_df_input = pd.read_csv(f"{root_path}{species}_instances.tsv", sep='\t', index_col=0)

    return  pval_rna_df, fc_rna_df, inst_df_input, pval_micro_df, fc_micro_df, has_micro

def combine_data(pval_df, fc_df, gene, isRNA=False, inst_df=None):
    # extract and combine data for each gene
    comb_df = pd.DataFrame()
    comb_df['sig'] = pval_df.index.tolist()
    comb_df['pval'] = pval_df[gene].tolist()
    comb_df['logpv'] = np.negative(np.log10(comb_df['pval']))
    comb_df['fc'] = fc_df[gene].tolist()
    if isRNA:
        comb_df['inst'] = comb_df['sig'].apply(lambda x: inst_df.loc[x, 'session_id'])
    return comb_df


def map_color(fc, pv):
    if fc < 0:
        return colors.to_hex(red_map(red_norm(pv)))
    elif fc == 0:
        return '#808080'
    else:
        return colors.to_hex(blue_map(blue_norm(pv)))

def make_plot(comb_df, species, gene, micro=False, micro_df=None):
    # create links from Bulk RNA-seq Appyter instance session IDs
    comb_df['inst'] = comb_df['inst'].apply(lambda x: f'https://appyters.maayanlab.cloud/Bulk_RNA_seq/{x}')

    # set color and size for each point on plot
    rna_colors = [map_color(r.fc, r.pval) for r in comb_df.itertuples()]
    rna_sizes = [12 if r.pval < 0.05 else 6 for r in comb_df.itertuples()]

    if micro:
        micro_colors = [map_color(r.fc, r.pval) for r in micro_df.itertuples()]
        micro_sizes = [12 if r.pval < 0.05 else 6 for r in micro_df.itertuples()]

    # generate data source
    data_source = ColumnDataSource(
        data=dict(
            x = comb_df['fc'],
            y = comb_df['logpv'],
            sig = comb_df['sig'],
            pval = comb_df['pval'], 
            fc = comb_df['fc'], 
            colors = rna_colors, 
            sizes = rna_sizes,
            label = ['RNA-seq']*comb_df.shape[0]
        )
    )

    # generate microarray data source if it exists
    if micro:
        micro_data_source = ColumnDataSource(
            data=dict(
                x = micro_df['fc'],
                y = micro_df['logpv'], 
                sig = micro_df['sig'],
                pval = micro_df['pval'], 
                fc = micro_df['fc'],
                colors = micro_colors,
                sizes = micro_sizes,
                label = ['Microarray']*micro_df.shape[0]
            )
        )
    # create hover tooltip
    tools = [
        ("Signature", "@sig"),
        ("P-Value", "@pval"),
        ("Fold Change", "@fc")
    ]
    # generate plot and relevant plot labels
    plot = figure(
        plot_width=700,
        plot_height=500,
        tooltips=tools
    )
    plot.circle(
        'x', 'y', 
        size='sizes',
        alpha=0.7, 
        line_alpha=0,
        line_width=0.01, 
        source=data_source,
        fill_color='colors', 
        name=f'{gene}_t2d_expression_volcano_plot',
        legend_group='label'
    )

    if micro:
        plot.square(
            'x', 'y',
            size='sizes',
            alpha=0.7,
            line_alpha=0,
            line_width=0.01,
            source=micro_data_source,
            fill_color='colors',
            name=f'{gene}_t2d_expression_volcano_plot',
            legend_group='label'
        )

    plot.xaxis.axis_label = 'log2(Fold Change)'
    plot.yaxis.axis_label = '-log10(P-value)'
    plot.title.text = f"Differential Expression of {gene} in {species} Type 2 Diabetes Transcriptomics Signatures"
    plot.title.align = 'center'
    plot.title.text_font_size = '14px'

    return plot

# create download link for table results
def download_link(df, fname, isRNA=False):
    if isRNA:
        df['Link to Bulk RNA-seq Analysis'] = df['Link to Bulk RNA-seq Analysis'].apply(
            lambda x: x.split('href=')[1].split('>')[0].replace('"', '')
        )

    df['Link to GEO Study'] = df['Link to GEO Study'].apply(
        lambda x: x.split('href=')[1].split('>')[0].replace('"', '')
    )

    df['Signature'] = df['Signature'].apply(lambda x: x.replace('* ', ''))
    csv = df.to_csv(fname, sep='\t', index=False)
    link = f'<div>Download full results: <a href="{fname}" target=_blank>{fname}</a></div>'
    return fname

# get GEO links 
@lru_cache()
def geo_link(sig_name, clickable):
    gse_id = sig_name.split('_')[0].replace('* ', '')
    geo_path = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
    if clickable:
        return f'<a target="_blank" href="{geo_path}{gse_id}">{gse_id}</a>'
    else:
        return f'{geo_path}{gse_id}'

@lru_cache()
def appyter_link(sig_name, inst=''):
    text = f'Analysis of {sig_name}'
    return f'<a target="_blank" href="{inst}">{text}</a>'

# create tables of significant results with links to GEO 
def make_tables(comb_df, species, gene, is_upreg, isRNA=False):
    root_path = 'static/searchdata/'
    if isRNA:
        sigranks = pd.read_feather(f"{root_path}all_{species}_fc_sigrank.f", columns=['index', gene]).set_index('index')
    else:
        sigranks = pd.read_feather(f"{root_path}{species}_affy_fc_sigrank.f", columns=['index', gene]).set_index('index')
    dir_df = comb_df[comb_df['fc'] > 0] if is_upreg else comb_df[comb_df['fc'] < 0]
    dir_df = dir_df.drop(columns='logpv').sort_values(by='pval', ascending=True)
    if isRNA:
        dir_df['inst'] = dir_df.apply(lambda row: appyter_link(row.sig, inst=row.inst), axis=1)
    dir_df['rank'] = [sigranks.loc[sig, gene] for sig in dir_df['sig']]
    if dir_df.shape[0] != 0:
        dir_df['sig'] = dir_df.apply(lambda row: f"* {row.sig}" if row.pval < 0.05 else row.sig, axis=1)
        dir_df['pval'] = dir_df['pval'].apply(lambda x: f'{x:.3e}')
        dir_df['fc'] = dir_df['fc'].apply(lambda x: f'{x:.4f}')
    if isRNA:
        dir_df = dir_df.rename(columns={'sig': 'Signature', 'pval': 'P-value', 'fc': 'Log2 Fold Change', 'inst': 'Link to Bulk RNA-seq Analysis', 'rank': 'Gene Rank in Signature'})
    else:
        dir_df = dir_df.rename(columns={'sig': 'Signature', 'pval': 'P-value', 'fc': 'Log2 Fold Change', 'rank': 'Gene Rank in Signature'})
    dir_df['Link to GEO Study'] = dir_df['Signature'].apply(geo_link, clickable=True)
    return dir_df

@lru_cache()
def send_plot(species, gene):
    pval_rna_df, fc_rna_df, inst_df_input, pval_micro_df, fc_micro_df, micro_exists = load_files(species, gene)
    comb_df_input = combine_data(pval_rna_df, fc_rna_df, gene, isRNA=True, inst_df=inst_df_input)
    if micro_exists:
        micro_df_input = combine_data(pval_micro_df, fc_micro_df, gene)
        plot = make_plot(comb_df_input, species, gene, micro=True, micro_df=micro_df_input)
    else:
        plot = make_plot(comb_df_input, species, gene)
    fname='static/searchdata/t2d-files/'
    up_comb_df_input = download_link(make_tables(comb_df_input, species, gene, is_upreg=True, isRNA=True), fname + gene + '_' + species + '_' + 'RNA-upreg.tsv', isRNA=True)
    dn_comb_df_input = download_link(make_tables(comb_df_input, species, gene, is_upreg=False, isRNA=True), fname + gene + '_' + species + '_' + 'RNA-dnreg.tsv', isRNA=True)
    if micro_exists and not(micro_df_input.empty):
        micro_df_input.dropna(inplace=True)
        up_micro_df = make_tables(micro_df_input, species, gene, is_upreg=True)
        up_micro_df_input = download_link(up_micro_df, fname + gene + '_' + species + '_' + 'micro-upreg.tsv')
        dn_micro_df = make_tables(micro_df_input, species, gene, is_upreg=False)
        dn_micro_df_input = download_link(dn_micro_df, fname + gene + '_' + species + '_' + 'micro-dnreg.tsv')
    else:
        micro_exists = False
        up_micro_df_input = ''
        dn_micro_df_input = ''
    return {'plot': json_item(plot, 'volcano-plot'), 'micro': micro_exists, 'tables': [up_comb_df_input, dn_comb_df_input, up_micro_df_input, dn_micro_df_input]}

def make_dge_plot(data, title, method, id_plot='dge-plot'):

    # set color and size for each point on plot

    if method == 'limma':
        colors = [map_color(r[1]['logFC'], r[1]['P.Value']) for r in data.iterrows()]
        sizes = [12 if r[1]['P.Value'] < 0.05 else 6 for r in data.iterrows()]
        data['logp'] = data['P.Value'].apply(lambda x: -np.log10(x))
        data['gene'] = data.index.values
    elif method == 'edgeR':
        colors = [map_color(r[1]['logFC'], r[1]['PValue']) for r in data.iterrows()]
        sizes = [12 if r[1]['PValue'] < 0.05 else 6 for r in data.iterrows()]
        data['logp'] = data['PValue'].apply(lambda x: -np.log10(x))
        data['gene'] = data.index.values
    elif method == 'DESeq2':
        colors = [map_color(r[1]['log2FoldChange'], r[1]['pvalue']) for r in data.iterrows()]
        sizes = [12 if r[1]['pvalue'] < 0.05 else 6 for r in data.iterrows()]
        data['logp'] = data['pvalue'].apply(lambda x: -np.log10(x))
        data['gene'] = data.index.values
    elif method == 'wilcoxon':
        colors = [map_color(r[1]['logfoldchanges'], r[1]['pvals']) for r in data.iterrows()]
        sizes = [12 if r[1]['pvals'] < 0.05 else 6 for r in data.iterrows()]
        data['logp'] = data['pvals'].apply(lambda x: -np.log10(x))
        data['gene'] = data.index.values


    # generate data source
    if method == 'limma':
        data_source = ColumnDataSource(
            data=dict(
                x = data['logFC'],
                y = data['logp'],
                gene =  data['gene'],
                t = data['t'],
                pval = data['P.Value'], 
                adjpval = data['adj.P.Val'],
                avgExpr = data['AvgExpr'], 
                B = data['B'],
                colors = colors, 
                sizes = sizes,
            )
        )
        tools = [
        ("Gene", "@gene"),
        ("t-test", "@t"),
        ("P-Value", "@pval"),
        ("Adj. P-Value", "@adjpval"),
        ("log2 Fold Change", "@x"),
        ("log-odds", "@B"),
        ("Avg. Expression", "@avgExpr")
    ]
    if method == 'edgeR':
        data_source = ColumnDataSource(
            data=dict(
                x = data['logFC'],
                y = data['logp'],
                logCPM = data['logCPM'],
                gene =  data['gene'],
                pval = data['PValue'], 
                FDR = data['FDR'],
                colors = colors, 
                sizes = sizes,
            )
        )
        tools = [
        ("Gene", "@gene"),
        ("P-Value", "@pval"),
        ("log2 Fold Change", "@x"),
        ("-log10(p)", "@y"),
        ("logCPM", "@logCPM")
    ]
    if method == 'DESeq2':
        data_source = ColumnDataSource(
            data=dict(
                x = data['log2FoldChange'],
                y = data['logp'],
                gene =  data['gene'],
                pval = data['pvalue'], 
                adjpval = data['padj'], 
                stat = data['stat'],
                lfcSE = data['lfcSE'],
                baseMean = data['baseMean'],
                colors = colors, 
                sizes = sizes,
            )
        )
        tools = [
        ("Gene", "@gene"),
        ("P-value", "@pval"),
        ("log2 Fold Change", "@x"),
        ("-log10(p)", "@y"),
        ("adj. P-value", "@adjpval"),
        ("stat", "@stat"),
        ("base Mean", "@baseMean"),
        ("lfcSE", "@lfcSE")
        ]
    if method == 'wilcoxon':
        data_source = ColumnDataSource(
            data=dict(
                x = data['logfoldchanges'],
                y = data['logp'],
                gene =  data['gene'],
                pval = data['pvals'], 
                adjpval = data['pvals_adj'], 
                scores = data['scores'],
                colors = colors, 
                sizes = sizes,
            )
        )
        tools = [
        ("Gene", "@gene"),
        ("P-value", "@pval"),
        ("log2 Fold Change", "@x"),
        ("-log10(p)", "@y"),
        ("adj. P-value", "@adjpval"),
        ("scores", "@scores")
        ]
    

    # create hover tooltip
    
    plot = figure(
        plot_width=700,
        plot_height=500,
        tooltips=tools
    )
    plot.circle(
        'x', 'y', 
        size='sizes',
        alpha=0.7, 
        line_alpha=0,
        line_width=0.01, 
        source=data_source,
        fill_color='colors', 
        name=None,
    )

    plot.xaxis.axis_label = 'log2(Fold Change)'
    plot.yaxis.axis_label = '-log10(P-value)'
    plot.title.text = f"Differential Gene Expression of {title}"
    plot.title.align = 'center'
    plot.title.text_font_size = '14px'
    print("made_plot")
    return json_item(plot, id_plot)



def str_to_int(string, mod):
    string = re.sub(r"\([^()]*\)", "", string).strip()
    byte_string = bytearray(string, "utf8")
    return int(hashlib.sha256(byte_string).hexdigest(), base=16)%mod

def make_single_visialization_plot(plot_df, values_dict,type, option_list,sample_names, caption_text, category_list_dict=None, location='right', category=True, dropdown=False, additional_info=None):

    # init plot 
    if additional_info is not None:
        source = ColumnDataSource(data=dict(x=plot_df["x"], y=plot_df["y"], values=values_dict[option_list[0]], 
                                        names=sample_names, info=additional_info[option_list[0]]))
    else:
        source = ColumnDataSource(data=dict(x=plot_df["x"], y=plot_df["y"], values=values_dict[option_list[0]], 
                                        names=sample_names))
    # node size
    if plot_df.shape[0] > 1000:
        node_size = 2
    else:
        node_size = 6
        
    if location == 'right':
        plot = figure(plot_width=700, plot_height=600)   
    else:
        plot = figure(plot_width=1000, plot_height=1000+20*len(category_list_dict[option_list[0]]))   
    if category == True:
        unique_category_dict = dict()
        for option in option_list:
            unique_category_dict[option] = sorted(list(set(values_dict[option])))
        
        # map category to color
        # color is mapped by its category name 
        # if a color is used by other categories, use another color
        factors_dict = dict()
        colors_dict = dict()
        for key in values_dict.keys():
            unused_color = list(Category20[20])
            factors_dict[key] = category_list_dict[key]
            colors_dict[key] = list()
            for category_name in factors_dict[key]:
                color_for_category = Category20[20][str_to_int(category_name, 20)]
                
                if color_for_category not in unused_color:
                    if len(unused_color) > 0:
                        color_for_category = unused_color[0]                        
                    else:
                        color_for_category = Category20[20][19]
                
                colors_dict[key].append(color_for_category)
                if color_for_category in unused_color:
                    unused_color.remove(color_for_category)
                    
        color_mapper = CategoricalColorMapper(factors=factors_dict[option_list[0]], palette=colors_dict[option_list[0]])
        legend = Legend()
        
        plot.add_layout(legend, location)
        scatter = plot.scatter('x', 'y', size=node_size, source=source, color={'field': 'values', 'transform': color_mapper}, legend_field="values")
        plot.legend.label_width = 30
        plot.legend.click_policy='hide'
        plot.legend.spacing = 1
        if location == 'below':
            location = 'bottom_left'
        plot.legend.location = location
        plot.legend.label_text_font_size = '10pt'
    # else:
    #     color_mapper = LinearColorMapper(palette=cc.CET_D1A, low=min(values_dict[option_list[0]]), high=max(values_dict[option_list[0]]))
    #     color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12)
    #     plot.add_layout(color_bar, 'right')
    #     plot.scatter('x', 'y', size=node_size,  source=source, color={'field': 'values', 'transform': color_mapper})
    
    if additional_info is not None:
            tooltips = [
            ("Sample", "@names"),
            ("Value", "@values"),
            ("p-value", "@info")
        ]
    else:
        tooltips = [
            ("Sample", "@names"),
            ("Value", "@values"),
        ]
    plot.add_tools(HoverTool(tooltips=tooltips))
    # plot.output_backend = "webgl"
    plot_name = ''
    if type == 'umap':
        plot.xaxis.axis_label = "UMAP_1"
        plot.yaxis.axis_label = "UMAP_2"
        plot_name = 'umap-plot'
    if type == 'tsne':
        plot.xaxis.axis_label = "tSNE_1"
        plot.yaxis.axis_label = "tSNE_2"
        plot_name = 'tsne-plot'

    if type == 'pca':
        plot.xaxis.axis_label = "PCA_1"
        plot.yaxis.axis_label = "PCA_2"
        plot_name = 'pca-plot'
    plot.xaxis.axis_label_text_font_size = "12pt"
    
    plot.yaxis.axis_label_text_font_size = "12pt"
    
    plot.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
    plot.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
    plot.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
    plot.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
    plot.xaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels
    plot.yaxis.major_label_text_font_size = '0pt'  # preferred method for removing tick labels
    
    return json_item(plot, plot_name)




def log2_normalize(x, offset=1.):
    return np.log2(x + offset)


@lru_cache()
def bulk_vis(expr_df, meta_df):

    meta_df = pd.read_csv(s3.open(meta_df), header=0, index_col=0, sep='\t')
    expr_df = pd.read_csv(s3.open(expr_df), header=0, index_col=0, sep='\t')

    expr_df.replace([np.inf, -np.inf],np.nan, inplace=True)
    expr_df = expr_df.dropna()

   
    expr_df = expr_df.transpose()
    expr_df = expr_df.dropna(axis=1)


    df_data_norm = log2_normalize(expr_df, offset=1)

    df_data_norm = quantile_normalize(df_data_norm, axis=0)

    #df_data_norm.replace([np.inf], np.nan, inplace=True)
    #df_data_norm = df_data_norm.dropna()
    #var = df_data_norm.var(axis = 0, numeric_only = True)
    
    return json_item(plot, plot_name)

 
@lru_cache()
def read_anndata_raw(path_to_data):
    return anndata.read_h5ad(s3.open(path_to_data))
@lru_cache()
def read_anndata_h5(path_to_data):
    return h5py.File(s3.open(path_to_data))