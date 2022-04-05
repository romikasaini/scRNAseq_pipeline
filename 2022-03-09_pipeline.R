source("https://raw.githubusercontent.com/romikasaini/scRNAseq_pipeline/main/funcs.R")

# Pipeline ####
path="./Samples/FH_6463_8/filtered_feature_bc_matrix"

y=SeuratSctype(path, save = F, plot=T, top20 = T, precise = T, res=0.8)
