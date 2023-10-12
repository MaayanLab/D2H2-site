from functools import lru_cache
import pandas as pd
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction
from maayanlab_bioinformatics.dge.logfc import logfc_differential_expression
from scipy.stats.mstats import zscore
from itertools import combinations
import warnings
import numpy as np
import scipy.stats as ss
from statsmodels.stats.multitest import multipletests
import s3fs
import scanpy as sc
import random
from helpers import read_anndata_h5, read_anndata_raw
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

base_url = os.environ.get('BASE_URL', 'd2h2/data')
endpoint = os.environ.get('ENDPOINT', 'https://minio.dev.maayanlab.cloud/')

s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': endpoint})


def qnormalization(data):

    X_quantile_norm = quantile_normalize(data)
    return X_quantile_norm


def CPM(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)

    return data


def logCPM(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)
        data = np.log2(data+1)

    # Return
    return data


def log(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = data.fillna(0)
        data = np.log2(data+1)

    return data


def get_precomputed_dge(sig, species):
    sig = [sig] + ['Unnamed: 0']
    df_fc = pd.read_csv(f'/Users/giacomomarino/D2H2-site/ETL/signatures/v1.1/rna_{species}_fc.csv', index_col=0, usecols=sig)
    df_fc.columns = ['logFC']
    df_pval = pd.read_csv(f'/Users/giacomomarino/D2H2-site/ETL/signatures/v1.1/rna_{species}_pval.csv', index_col=0, usecols=sig)
    df_pval.columns = ['pval']
    dge_df = pd.merge(left=df_fc, right=df_pval, right_index=True, left_index=True)
    dge_df['adj.pval'] = multipletests(dge_df['pval'].values)[1]
    dge_df.sort_values('adj.pval', inplace=True)
    return dge_df

def get_precomputed_dge_options(gse, species):
    df = pd.read_feather(f'/Users/giacomomarino/D2H2-site/ETL/signatures/v1.1/rna_{species}_pval.f', columns=['index'])
    precomputed_sigs = [sig for sig in df['index'].values if sig.split('-')[0] == gse]
    return precomputed_sigs


def get_signatures(classes, dataset, normalization, method, meta_class_column_name, filter_genes):

    expr_df = dataset['rawdata']
    if filter_genes == True:
        expr_df = dataset['rawdata+filter_genes']

    signatures = dict()

    for cls1, cls2 in combinations(classes, 2):
        cls1_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"]
                                                          [meta_class_column_name] == cls1, :].index.tolist()  # control
        cls2_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"]
                                                          [meta_class_column_name] == cls2, :].index.tolist()  # case
        signature_label = " vs. ".join([cls1, cls2])
        if method == "characteristic_direction":
            expr_df.dropna(inplace=True)
            signature = characteristic_direction(expr_df.loc[:, cls1_sample_ids], expr_df.loc[:, cls2_sample_ids])
            signature['P-value'] = ss.norm.sf(np.abs(zscore(signature['CD-coefficient'])))
            signature['LogFC'] = logfc_differential_expression(expr_df.loc[:, cls1_sample_ids], expr_df.loc[:, cls2_sample_ids]).loc[signature.index.values]['LogFC']
        elif method == "DESeq2":
            dds = DeseqDataSet(
            counts=expr_df.T,
            clinical=dataset['dataset_metadata'],
            design_factors="Condition",
            refit_cooks=True,
            n_cpus=2,
            )
            dds.deseq2()
            stat_res = DeseqStats(dds, n_cpus=2)
            stat_res.summary()
            signature = stat_res.results_df.sort_values("padj", ascending=True)

        signatures[signature_label] = signature

    return signatures, signature_label


def normalize(dataset, current_dataset, logCPM_normalization, log_normalization, z_normalization, q_normalization):
    normalization = current_dataset
    if logCPM_normalization == True:
        data = dataset[normalization]
        normalization += '+logCPM'
        dataset[normalization] = logCPM(data)

    if log_normalization == True:
        data = dataset[normalization]
        normalization += '+log'
        dataset[normalization] = log(data)

    if z_normalization == True:
        data = dataset[normalization]
        normalization += '+z_norm'
        dataset[normalization] = data.T.apply(ss.zscore, axis=0).T.dropna()

    if q_normalization == True:
        data = dataset[normalization]
        normalization += '+q_norm'
        dataset[normalization] = qnormalization(data)
    return dataset, normalization


def check_df(df, col):
    if col not in df.columns:
        raise IOError


def compute_dge(rnaseq_data_filename, meta_data_filename, diff_gex_method, control_name, perturb_name, logCPM_normalization, log_normalization, z_normalization, q_normalization):
    meta_class_column_name = 'Condition'

    meta_df = pd.read_csv(s3.open(meta_data_filename),
                          sep="\t", index_col=0, dtype=str)

    if len(set(meta_df['Group'])) > 1:
        meta_df['Combined'] = meta_df['Condition'] + ' ' + meta_df['Group']
        meta_class_column_name = 'Combined'
    meta_df = meta_df[meta_df[meta_class_column_name].isin([control_name, perturb_name])]

    meta_df.index = meta_df.index.map(str)
    expr_df = pd.read_csv(s3.open(rnaseq_data_filename), index_col=0, sep="\t").sort_index()
    expr_df = expr_df.loc[expr_df.sum(axis=1) > 0, :]

    # Match samples between the metadata and the datasets
    try:
        check_df(meta_df, meta_class_column_name)
    except:
        print(f"Error! Column '{meta_class_column_name}' is not in metadata")

    meta_df = meta_df[meta_df.index.isin(expr_df.columns)]
    low_expression_threshold = .3
    logCPM_normalization = True
    log_normalization = False
    z_normalization = True
    q_normalization = False

    classes = list(meta_df[meta_class_column_name].unique())

    classes.remove(control_name)
    classes.insert(0, control_name)
    meta_df['tmp_class'] = pd.Categorical(
        meta_df[meta_class_column_name], classes)
    meta_df = meta_df.sort_values('tmp_class')
    meta_df = meta_df.drop('tmp_class', axis=1)

    expr_df = expr_df.loc[:, meta_df.index]
    expr_df = expr_df.groupby(expr_df.index).sum()

    assert (meta_df.shape[0] == expr_df.shape[1])

    dataset = dict()
    current_dataset = 'rawdata'
    dataset[current_dataset] = expr_df
    filter_genes = True

    # Filter out lowly expressed genes
    mask_low_vals = (expr_df > low_expression_threshold).sum(axis=1) > 2
    expr_df = expr_df.loc[mask_low_vals, :]
    current_dataset += '+filter_genes'
    dataset[current_dataset] = expr_df

    dataset['dataset_metadata'] = meta_df

    dataset, normalization = normalize(
        dataset, current_dataset, logCPM_normalization, log_normalization, z_normalization, q_normalization)

    signatures, signature_label = get_signatures(
        classes, dataset, normalization, diff_gex_method, meta_class_column_name, filter_genes)

    return signatures[signature_label], signature_label


########## SINGLE CELL DGE METHODS ###########
def get_signatures_single(classes, expr_file, method, meta_class_column_name, cluster=True, filter_genes=True, aggregate=False):

    # Getting the same number of samples from each cluster to use for diffrential gene expression.
    f = read_anndata_h5(expr_file)

    clus_numbers = f["var/leiden/codes"][:]
    leiden_data_vals = list(map(lambda x: "Cluster " + str(x), clus_numbers))

    metadata_dict_counts = pd.Series(leiden_data_vals).value_counts()
    genes = np.array(f['obs/gene_symbols'][:].astype(str))
    cells = f['var/column_names'][:].astype(str)
    cluster_list = []
    if cluster == True and aggregate == True:
        num_to_sample = min(min(metadata_dict_counts), 15)
        list_of_adata = []
        for cls in metadata_dict_counts.keys():

            leiden_data_vals = pd.Series(leiden_data_vals)
            cls_leiden_vals = leiden_data_vals[leiden_data_vals == cls]
            idx = list(sorted(random.sample(
                list(cls_leiden_vals.index.values), k=num_to_sample)))
            adata_sample = pd.DataFrame(
                f['raw/X'][:, idx], index=genes, columns=cells[idx])

            cluster_list += [cls]*num_to_sample
            list_of_adata.append(adata_sample)

        full_adata = pd.concat(list_of_adata, axis=1)
        # Need to repeat for the normalized data.
        expr_df = full_adata
        raw_expr_df = full_adata
        cluster_list = pd.Series(cluster_list)

    signatures = dict()

    if cluster == True and aggregate == True:

        for cls1 in classes:
            signature_label = f"{cls1} vs. rest"
            print("Analyzing.. {} using {}".format(signature_label, method))
            cols = pd.Series(full_adata.columns)
            clus_idx = cluster_list[cluster_list == cls1].index.tolist()
            rest_idx = cluster_list[cluster_list != cls1].index.tolist()
            cls1_sample_ids = list(cols[clus_idx].astype(str))  # case
            non_cls1_sample_ids = list(cols[rest_idx].astype(str)) # control

            tmp_raw_expr_df = raw_expr_df

            tmp_raw_expr_df.columns.name = "Sample"

            if method == "characteristic_direction":
                signature = characteristic_direction(expr_df.loc[:, non_cls1_sample_ids], expr_df.loc[:, cls1_sample_ids]).dropna()
                signature['P-value'] = ss.norm.sf(np.abs(zscore(signature['CD-coefficient'])))
                signature['LogFC'] = logfc_differential_expression(expr_df.loc[:, non_cls1_sample_ids], expr_df.loc[:, cls1_sample_ids]).loc[signature.index.values]['LogFC']
            elif method == "DESeq2":
                condition_labels = ['C'] * expr_df.loc[:, non_cls1_sample_ids].shape[1] + ['RS'] *  expr_df.loc[:, cls1_sample_ids].shape[1]
                sample_names = expr_df.loc[:, non_cls1_sample_ids].columns.tolist() + expr_df.loc[:, cls1_sample_ids].columns.tolist()
                metadata = pd.DataFrame({'Sample': sample_names, 'Condition': condition_labels}).set_index("Sample")
                dds = DeseqDataSet(
                    counts=expr_df[non_cls1_sample_ids + cls1_sample_ids].T,
                    clinical=metadata,
                    design_factors="Condition",
                    refit_cooks=True
                )
                dds.deseq2()
                stat_res = DeseqStats(dds)
                stat_res.summary()
                signature = stat_res.results_df.sort_values("padj", ascending=True).dropna()
            elif method == "wilcoxon":
                dataset = read_anndata_raw(expr_file)
                dataset_raw = dataset.raw.to_adata()
                # Passing the gene information as cells x genes for the differential expression method.
                dataset = dataset.T
                dataset.raw = dataset_raw.T
                if 'log1p' in dataset.uns.keys():
                    del dataset.uns['log1p']
                sc.tl.rank_genes_groups(
                    dataset, meta_class_column_name, method='wilcoxon', use_raw=True)
                dedf = sc.get.rank_genes_groups_df(dataset, group=cls1).set_index(
                    'names').sort_values('pvals', ascending=True)
                dedf = dedf.replace([np.inf, -np.inf], np.nan).dropna()
                dedf = dedf.sort_values("logfoldchanges", ascending=False)
                signature = dedf.dropna()
            signatures[signature_label] = signature
    return signatures


def compute_dge_single(expr_file, diff_gex_method, enrichment_groupby, meta_class_column_name, clustergroup, agg):

    meta_class_column_name = "leiden"
    classes = [clustergroup]
    bool_cluster = True

    signatures = get_signatures_single(classes, expr_file, method=diff_gex_method,
                                       meta_class_column_name=meta_class_column_name, cluster=bool_cluster, aggregate=agg)
    return signatures
