source("./funcs.R")

# Pipeline ####
path="./Samples/FH_6385_2/filtered_feature_bc_matrix"

#Seurat pipeline and annotation using scType

SeuratSctype <- function(path){
  sample=ReadScData(path = path)
  sample=SampleQC(sample)
  z=GetSampleMetrics(sample, plot=T) #change here to save
  sample=Cluster(sample, res=0.8)
  x=FindLog2FC(sample, as.df = T)
  es=clustScore(Log2FCdata = x[[2]])
  sample=ClustUMAP(sample, es, plot = T, save=T) #scType UMAP
  CPM=CalcCPM(sample, clustname = "customclassif")
  Raw=CalcRawCount(sample, clustname =  "customclassif")
  return(list(Metrics=z, sample=sample, GE=list(CPM=CPM, Raw=Raw)))
}

y=SeuratSctype(path)


#Seurat pipeline and annotation using multimodal reference mapping
SeuratMultimodal <- function(path){
  sample=ReadScData(path = path)
  sample=SampleQC(sample)
  z=GetSampleMetrics(sample, plot=T)
  sample=Cluster(sample, res=0.8)
  x=FindLog2FC(sample, as.df = T)
  es=clustScore(Log2FCdata = x[[2]])
  sample <- Multimodal_UMAP(sample, plot = T, save=T) #multimodal reference mapping UMAP
  CPM=CalcCPM(sample, clustname = "predicted.celltype")
  Raw=CalcRawCount(sample, clustname =  "predicted.celltype")
  return(list(Metrics=z, sample=sample, GE=list(CPM=CPM, Raw=Raw)))
}
x=SeuratMultimodal(path)
