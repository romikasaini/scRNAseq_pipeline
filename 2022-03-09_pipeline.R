source("./funcs.R")

# Pipeline ####
path="./Samples/FH_6385_2/filtered_feature_bc_matrix"
path="../New_research_topics/single_cell_sequencing/multiple_myeloma_2021_originalfiles/Batch1-8_090120-120121/count_210215_A00464_0283_AHYVWLDRXX/FH_6385_2/filtered_feature_bc_matrix/"

#Seurat pipeline and annotation using scType

SeuratSctype <- function(path){
  sample=ReadScData(path = path)
  sample=SampleQC(sample)
  z=GetSampleMetrics(sample, plot=T) #change here to save
  sample=Cluster(sample)
  x=FindLog2FC(sample, as.df = T)
  es=clustScore(Log2FCdata = x[[2]])
  sample=ClustUMAP(sample, es, plot = T, save=T) #scType UMAP
  CPM=CalcCPM(sample, clustname = "customclassif")
  Raw=CalcRawCount(sample, clustname =  "customclassif")
}

SeuratSctype(path)


#Seurat pipeline and annotation using multimodal reference mapping
SeuratMultimodal <- function(path){
  sample=ReadScData(path = path)
  sample=SampleQC(sample)
  z=GetSampleMetrics(sample, plot=T)
  sample=Cluster(sample)
  x=FindLog2FC(sample, as.df = T)
  es=clustScore(Log2FCdata = x[[2]])
  sample <- Multimodal_UMAP(sample, plot = T, save=T) #multimodal reference mapping UMAP
  CPM=CalcCPM(sample, clustname = "predicted.celltype")
  Raw=CalcRawCount(sample, clustname =  "predicted.celltype")
}
SeuratMultimodal(path)
