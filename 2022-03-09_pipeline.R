source("./funcs.R")

# Pipeline ####
path="./Samples/FH_6463_8/filtered_feature_bc_matrix"

#Seurat pipeline and annotation using scType

SeuratSctype <- function(path){
  sample=ReadScData(path = path)
  sample=SampleQC(sample)
  z=GetSampleMetrics(sample, plot=T) #change here to save
  sample=Cluster(sample, res=0.8)
  x=FindLog2FC(sample, as.df = T)
  x[[1]] %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
  es=clustScore(Log2FCdata = x[[2]])
  sample=ClustUMAP(sample, es, plot = T, save=F, precise=T) # select precise = F for 1 step annotation and = T for 2-step annotation
  CPM=CalcCPM(sample, clustname = "customclassif")
  Raw=CalcRawCount(sample, clustname =  "customclassif")
  return(list(Metrics=z, sample=sample, GE=list(CPM=CPM, Raw=Raw), Log2FCperClust=top20))
}

y=SeuratSctype(path)
