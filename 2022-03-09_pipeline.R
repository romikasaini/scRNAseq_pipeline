source("https://raw.githubusercontent.com/romikasaini/scRNAseq_pipeline/main/funcs.R")

# Pipeline ####
#path="./Samples/FH_6463_8/filtered_feature_bc_matrix"
#path="../Immunocore/data/FHRB_743_6/filtered_feature_bc_matrix/"
path = "./data/pmbc/filtered_feature_bc_matrix/"

# This command takes path to sc-RNA data and produces result list containing metrics (cell cycle percentage, ribosome share (range and percentage), mean mitochondrial percentage and erythrocyte percentage creates density plots for ribosomal, hemoglobin and mitochondrial genes), final seurat object, mean and standard deviation values of CPM and raw counts per cluster of in CPM values and top20 genes for each cluster.
# parameters: path to sc-RNA data
# Arguments
# plot: Plots the UMAP in Rstudio. Boolean, default is True.
# save: Save the UMAP in the working directory. Boolean, default is False.
# precise: Boolean, default is False. select precise = F for 1 step annotation and = T for 2-step annotation
# res: Resolution parameter to set the granularity of the clustering, default is 0.8.
# top20: Boolean, default is true. True produces a dataframe with top 20 genes in each cluster based on log2foldchange.

y = SeuratSctype(
  path,
  save = T,
  plot = T,
  top20 = T,
  precise = T,
  res = 0.8,
  N_dims = 10,
  UMAP_dims = 15
)
