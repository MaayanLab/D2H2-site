# R script to download selected samples
# Copy code and run on a local machine to initiate download

library("rhdf5") # can be installed using Bioconductor
library("GEOquery")
library("rjson")

organism = "mouse"
wd = paste("/Users/giacomomarino/D2H2-site/static/data/", organism, sep="")

setwd(wd)


destination_file = paste(organism,"_matrix_v11.h5", sep="")

studies = c("GSE101207")

samples = h5read(destination_file, "meta/samples/geo_accession")
genes = h5read(destination_file, "meta/genes/genes")


for (i in 1:length(studies)) {
  study = studies[i]
  gds <- getGEO(study)
  print(gds)
  samp = sampleNames(gds)
  if (!dir.exists(study)) {
    dir.create(study)
  }
  
  extracted_expression_file = paste(study, "_Expression.txt", sep="")
  extracted_meta_file = paste(study, "/", study, "_Metadata.txt", sep="")
  # Identify columns to be extracted
  sample_locations = which(samples %in% samp)
  
  # extract gene expression from compressed data
  expression = t(h5read(destination_file, "data/expression", index=list(sample_locations, 1:length(genes))))
  H5close()
  rownames(expression) = genes
  colnames(expression) = samples[sample_locations]
  
  # Print file
  write.table(expression, file=paste(study, "/", extracted_expression_file, sep=""), sep="\t", quote=FALSE, col.names=NA)
  print(paste0("Expression file was created at ", paste(getwd(),"/", study, sep=""), "/", extracted_expression_file))
}










