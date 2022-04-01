source("./funcs.R")

# Pipeline ####
path="./Samples/FH_4263_2/filtered_feature_bc_matrix"
sample=ReadScData(path = path)
sample=SampleQC(sample)
z=GetSampleMetrics(sample, plot=T, save = F)
sample=Cluster(sample)
x=FindLog2FC(sample, as.df = T)

es=clustScore(Log2FCdata = x[[2]])
sample=ClustUMAP(sample, es, plot = T)
CPM=CalcCPM(seu1, clustname = "customclassif")
Raw=CalcRawCount(seu1, clustname =  "customclassif")
