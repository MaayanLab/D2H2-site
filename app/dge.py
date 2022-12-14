from functools import lru_cache
import pandas as pd
import rpy2
from rpy2 import robjects
from rpy2.robjects import r, pandas2ri
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction
from maayanlab_bioinformatics.dge.limma_voom import limma_voom_differential_expression
from itertools import combinations
import warnings
import numpy as np
import scipy.stats as ss
import s3fs
import scanpy as sc
import anndata
import random
from helpers import read_anndata_h5, read_anndata_raw


s3 = s3fs.S3FileSystem(anon=True, client_kwargs={'endpoint_url': 'https://minio.dev.maayanlab.cloud/'})

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

def get_signatures(classes, dataset, normalization, method, meta_class_column_name, filter_genes):
    
	
    robjects.r('''limma <- function(rawcount_dataframe, design_dataframe, filter_genes=FALSE, adjust="BH") {
        # Load packages
        suppressMessages(require(limma))
        suppressMessages(require(edgeR))
        # Convert design matrix
        design <- as.matrix(design_dataframe)
        
        # Create DGEList object
        dge <- DGEList(counts=rawcount_dataframe)
        # Filter genes
        if (filter_genes) {
            keep <- filterByExpr(dge, design)
            dge <- dge[keep,]
        }
        # Calculate normalization factors
        dge <- calcNormFactors(dge)
        # Run VOOM
        v <- voom(dge, plot=FALSE)
        # Fit linear model
        fit <- lmFit(v, design)
        # Make contrast matrix
        cont.matrix <- makeContrasts(de=B-A, levels=design)
        # Fit
        fit2 <- contrasts.fit(fit, cont.matrix)
        # Run DE
        fit2 <- eBayes(fit2)
        # Get results
        limma_dataframe <- topTable(fit2, adjust=adjust, number=nrow(rawcount_dataframe))
        
        # Return
        results <- list("limma_dataframe"= limma_dataframe, "rownames"=rownames(limma_dataframe))
        return (results)
    }''')

    robjects.r('''edgeR <- function(rawcount_dataframe, g1, g2) {
    # Load packages
    suppressMessages(require(limma))
    suppressMessages(require(edgeR))
    
    colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
    rownames(colData) <- c(g1,g2)
    colnames(colData) <- c("group")
    colData$group = relevel(as.factor(colData$group), "Control")
    
    y <- DGEList(counts=rawcount_dataframe, group=colData$group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    res <- topTags(et, n=Inf)
    # Return
    res <- as.data.frame(res)
    results <- list("edgeR_dataframe"= res, "rownames"=rownames(res))
    return (results) 
    }
    ''')
    
    robjects.r('''deseq2 <- function(rawcount_dataframe, g1, g2) {
    # Load packages
    suppressMessages(require(DESeq2))
    colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
    rownames(colData) <- c(g1,g2)
    colnames(colData) <- c("group")
    colData$group = relevel(as.factor(colData$group), "Control")
    dds <- DESeqDataSetFromMatrix(countData = rawcount_dataframe, colData = colData, design=~(group))

    dds <- DESeq(dds)
    res <- results(dds)
    
    res[which(is.na(res$padj)),] <- 1
    res <- as.data.frame(res)
    
    results <- list("DESeq_dataframe"= res, "rownames"=rownames(res))
    return(results)
    
    
    
    }
    ''')
    pandas2ri.activate()

    tmp_normalization = normalization.replace("+z_norm+q_norm","").replace("+z_norm","")
    raw_expr_df = dataset['rawdata']
    expr_df = dataset['rawdata']
    if filter_genes == True:
        expr_df = dataset['rawdata+filter_genes']
        
    signatures = dict()

    for cls1, cls2 in combinations(classes, 2):
        cls1_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls1, :].index.tolist() #control
        cls2_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls2,:].index.tolist() #case
        
        signature_label = " vs. ".join([cls1, cls2])
        
        if method == "limma":
            limma = robjects.r['limma']

            design_dataframe = pd.DataFrame([{'index': x, 'A': int(x in cls1_sample_ids), 'B': int(x in cls2_sample_ids)} for x in raw_expr_df.columns]).set_index('index')

            processed_data = {"expression": raw_expr_df, 'design': design_dataframe}
            
            expr_r = pandas2ri.py2rpy(processed_data['expression'])
            design_r = pandas2ri.py2rpy(processed_data['design'])
            limma_results = pandas2ri.rpy2py(limma(expr_r, design_r))

            signature = pd.DataFrame(limma_results[0])
            signature.columns = ['logFC', 'AvgExpr', 't', 'P.Value', 'adj.P.Val', 'B']
            signature.index = limma_results[1]
            signature = signature.sort_values("t", ascending=False)
            
        elif method == "characteristic_direction":
            signature = characteristic_direction(dataset[tmp_normalization].loc[:, cls1_sample_ids], dataset[normalization].loc[:, cls2_sample_ids], calculate_sig=True)
            signature = signature.sort_values("CD-coefficient", ascending=False)
        elif method == "edgeR":
            edgeR = robjects.r['edgeR']
            
            edgeR_results = pandas2ri.conversion.rpy2py(edgeR(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))
            
            signature = pd.DataFrame(edgeR_results[0])
            signature.index = edgeR_results[1]
            signature = signature.sort_values("logFC", ascending=False)
        elif method == "DESeq2":
            # deseq2 receives raw counts
            DESeq2 = robjects.r['deseq2']
            DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))
            
            signature = pd.DataFrame(DESeq2_results[0])
            signature.index = DESeq2_results[1]
            signature = signature.sort_values("log2FoldChange", ascending=False)
                        
            
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



@lru_cache
def compute_dge(rnaseq_data_filename, meta_data_filename, diff_gex_method, control_name, perturb_name, logCPM_normalization, log_normalization, z_normalization, q_normalization):
	meta_class_column_name = 'Condition'

	meta_df = pd.read_csv(s3.open(meta_data_filename), sep="\t", index_col=0, dtype=str)
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
	meta_df['tmp_class'] = pd.Categorical(meta_df[meta_class_column_name], classes)
	meta_df = meta_df.sort_values('tmp_class')
	meta_df = meta_df.drop('tmp_class', axis=1)


	expr_df = expr_df.loc[:,meta_df.index]
	expr_df = expr_df.groupby(expr_df.index).sum()

	assert(meta_df.shape[0]==expr_df.shape[1])

	dataset = dict()
	current_dataset = 'rawdata'
	dataset[current_dataset] = expr_df
	filter_genes = True

	## Filter out lowly expressed genes
	mask_low_vals = (expr_df > low_expression_threshold).sum(axis=1) > 2
	expr_df = expr_df.loc[mask_low_vals, :]
	current_dataset += '+filter_genes'
	dataset[current_dataset] = expr_df 

	dataset['dataset_metadata'] = meta_df



	dataset, normalization = normalize(dataset, current_dataset, logCPM_normalization, log_normalization, z_normalization, q_normalization)

	signatures, signature_label = get_signatures(classes, dataset, normalization, diff_gex_method, meta_class_column_name, filter_genes)

	return signatures[signature_label], signature_label


########## SINGLE CELL DGE METHODS ###########
def get_signatures_single(classes, expr_file, method, meta_class_column_name, cluster=True, filter_genes=True, aggregate = False):
           
    robjects.r('''edgeR <- function(rawcount_dataframe, g1, g2) {
        # Load packages
        suppressMessages(require(limma))
        suppressMessages(require(edgeR))
        colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
        rownames(colData) <- c(g1,g2)
        colnames(colData) <- c("group")
        colData$group = relevel(as.factor(colData$group), "Control")
        y <- DGEList(counts=rawcount_dataframe, group=colData$group)
        y <- calcNormFactors(y)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        et <- exactTest(y)
        res <- topTags(et, n=Inf)
        # Return
        res <- as.data.frame(res)
        results <- list("edgeR_dataframe"= res, "rownames"=rownames(res), "colnames"=colnames(res))
        return (results)
    }
    ''')

    robjects.r('''deseq2 <- function(rawcount_dataframe, g1, g2) {
        # Load packages
        suppressMessages(require(DESeq2))
        colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
        rownames(colData) <- c(g1,g2)
        colnames(colData) <- c("group")
        colData$group = relevel(as.factor(colData$group), "Control")
        dds <- DESeqDataSetFromMatrix(countData = rawcount_dataframe, colData = colData, design=~(group))
        dds <- DESeq(dds)
        res <- results(dds)
        res[which(is.na(res$padj)),] <- 1
        res <- as.data.frame(res)
        print(dim(res))
        results <- list("DESeq_dataframe"= res, "rownames"=rownames(res), "colnames"=colnames(res))
        return(results)
    }
    ''')
    pandas2ri.activate()
    
    # expr_df = dataset.to_df().T
    # raw_expr_df = dataset.raw.to_adata().to_df().T
    # meta_df = dataset.obs
    #Getting the same number of samples from each cluster to use for diffrential gene expression. 
    f = read_anndata_h5(expr_file)
    
    leiden_data = f["var/leiden/categories"][:].astype(str)
    clus_numbers = f["var/leiden/codes"][:]
    leiden_data_vals = list(map(lambda x: "Cluster " + str(x), clus_numbers))
    classes = sorted(leiden_data)
    classes = sorted(classes, key=lambda x: int(x.replace("Cluster ", "")))
	#Stores the number of of cells correlated to each cluster. 
    metadata_dict_counts = pd.Series(leiden_data_vals).value_counts()
    genes = np.array(f['obs/gene_symbols'][:].astype(str))

    if cluster == True and aggregate == True:
        #adata_raw = dataset.raw.to_adata()
        num_to_sample = min(min(metadata_dict_counts),15)
        list_of_adata = []
        for cls in metadata_dict_counts.keys():
            leiden_data_vals = pd.Series(leiden_data_vals)
            cls_leiden_vals = leiden_data_vals[leiden_data_vals == cls]
            idx = list(sorted(random.sample(list(cls_leiden_vals.index.values), k=num_to_sample)))
            adata_sample = pd.DataFrame(f['raw/X'][:, idx], index=genes, columns=[cls]*num_to_sample)
            #adata_sub = adata_raw[adata_raw.obs['leiden'] == cls]
            #adata_sample = sc.pp.subsample(adata_sub, n_obs=num_to_sample, copy=True)
            list_of_adata.append(adata_sample)
        full_adata = pd.concat(list_of_adata, axis=1)
        #Need to repeat for the normalized data. 
        print(full_adata)
        expr_df = full_adata.T
        raw_expr_df = full_adata.T
        meta_df = leiden_data_vals
        # print(meta_df.head())
        print(meta_df.shape)
    
    signatures = dict()
    #GETTING THE TOP GENES RATHER THAN DOING DGE WITH CERTAIN METHODS
    """if method == 'none':
        for cls1 in classes:
            signature_label = f"{cls1} vs. rest"
            print("Analyzing.. {} using {}".format(signature_label, method))
            raw_adata = dataset.raw.to_adata()
            cluster_data = raw_adata[raw_adata.obs[meta_class_column_name] == cls1]
            df = cluster_data.to_df()
            df = df.divide(df.sum(axis=1), axis=0)
            mean_df = df.mean().sort_values(ascending=False)
            yaxis_gene_names = list(mean_df.index)[:100]
            print(yaxis_gene_names)
            print(cluster_data.shape)
            signatures[signature_label] = yaxis_gene_names
        return signatures """
    
    if cluster == True and aggregate == True:
        
            
        for cls1 in classes:
            signature_label = f"{cls1} vs. rest"
            print("Analyzing.. {} using {}".format(signature_label, method))
            cols = pd.Series(full_adata.columns)
            cls1_sample_ids = cols[cols == cls1].index.tolist() #case
            non_cls1_sample_ids = cols[cols != cls1].index.tolist()#control
            sample_ids = non_cls1_sample_ids.copy()
            sample_ids.extend(cls1_sample_ids)
            tmp_raw_expr_df = raw_expr_df
            if method == "limma":
                signature = limma_voom_differential_expression(tmp_raw_expr_df.iloc[:, non_cls1_sample_ids], tmp_raw_expr_df.iloc[:, cls1_sample_ids])
                signature.rename(columns={"AveExpr": "AvgExpr"}, inplace=True)
            elif method == "characteristic_direction":
                signature = characteristic_direction(expr_df.iloc[:, non_cls1_sample_ids], expr_df.iloc[:, cls1_sample_ids], calculate_sig=False)
            elif method == "edgeR":
                edgeR = robjects.r['edgeR']
                edgeR_results = pandas2ri.conversion.rpy2py(edgeR(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(non_cls1_sample_ids), pandas2ri.conversion.py2rpy(cls1_sample_ids)))

                signature = pd.DataFrame(edgeR_results[0])
                if(signature.shape[0] == 4):
                    signature = signature.T
                signature.index = edgeR_results[1]
                signature.columns = edgeR_results[2]
                signature = signature.sort_values("logFC", ascending=False)
            elif method == "DESeq2":
                DESeq2 = robjects.r['deseq2']
                DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(tmp_raw_expr_df), pandas2ri.conversion.py2rpy(non_cls1_sample_ids), pandas2ri.conversion.py2rpy(cls1_sample_ids)))
                signature = pd.DataFrame(DESeq2_results[0])
                if(signature.shape[0] == 6):
                    signature = signature.T
                signature.index = DESeq2_results[1]
                signature.columns = DESeq2_results[2]
                signature = signature.sort_values("log2FoldChange", ascending=False)
            elif method == "wilcoxon":  
                dataset = read_anndata_raw(expr_file) 
                if 'log1p' in dataset.uns.keys():
		            del dataset.uns['log1p']
                sc.tl.rank_genes_groups(dataset, meta_class_column_name, method='wilcoxon', use_raw=True)
                dedf = sc.get.rank_genes_groups_df(dataset, group=cls1).set_index('names').sort_values('pvals', ascending=True)
                dedf = dedf.replace([np.inf, -np.inf], np.nan).dropna()              
                dedf = dedf.sort_values("logfoldchanges", ascending=False)
                signature = dedf
                
            signatures[signature_label] = signature
    return signatures


#compute_dge_single(adata, diff_gex_method, 'Cluster', 'leiden', True)
def compute_dge_single(expr_file, diff_gex_method, enrichment_groupby, meta_class_column_name,clustergroup, agg):
    if diff_gex_method == "characteristic_direction":
        fc_colname = "CD-coefficient"
        sort_genes_by = "CD-coefficient"
        ascending = False
    elif diff_gex_method == "limma":
        fc_colname = "logFC"
        sort_genes_by = "t"
        ascending = False
    elif diff_gex_method == "edgeR":
        fc_colname = "logFC"
        sort_genes_by = "PValue"
        ascending = True
    elif diff_gex_method == "DESeq2":
        fc_colname = "log2FoldChange"
        sort_genes_by = "padj"
        ascending = True
    elif diff_gex_method == "wilcoxon":
        fc_colname = "logfoldchanges"
        sort_genes_by = "scores"
        ascending = False
        
    #f = read_anndata_h5(expr_file)
        

    meta_class_column_name = "leiden"
    classes = [clustergroup]
    bool_cluster=True
        
    signatures = get_signatures_single(classes, expr_file, method=diff_gex_method, meta_class_column_name=meta_class_column_name, cluster=bool_cluster, aggregate = agg)
    
    return signatures