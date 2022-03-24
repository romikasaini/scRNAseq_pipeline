source("./funcs.R")

# Pipeline ####
path="./Samples/FH_6463_8/filtered_feature_bc_matrix"
sample=ReadScData(path = path)
seu1=SampleQC(sample)
#Use both plots for graphical QC representation
z=GetSampleMetrics(sample, plot=T, save = T, plotname = "beforeQC")
c=GetSampleMetrics(seu1, plot=T, save=T, plotname = "afterQC")

## INSERT CLUSTERING FUNCTIONS HERE ##

x=CalcCPM(seu1, clustname = "customclassif")
y=CalcRawCount(seu1, clustname =  "customclassif")


#scaling to UMAP

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE) # less than a minute
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:15) #Define here the number of wanted PC's in dim= ! 
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)


#check the plot
#DimPlot(seurat_obj_list[[i]], reduction = "umap", label = T)

#Finding differentially expressed features (cluster biomarkers) should we skip this if we are using sctype?
#Differential gene expression? should we skip this if we are using sctype?

# Cell type assignment from scType
# Quick start -------------------------------------------------------------

# load libraries and functions
####install.packages("HGNChelper")
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# assign cell types
scRNAseqData = readRDS(gzcon(url('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/exampleData.RDS'))); #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)

# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

for (i in 1:length(seurat_obj_list)){
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seurat_obj_list[[i]][["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj_list[[i]]@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj_list[[i]]@meta.data[seurat_obj_list[[i]]@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj_list[[i]]@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  # Overlay the identified cell types on UMAP plot
  seurat_obj_list[[i]]@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat_obj_list[[i]]@meta.data$customclassif[seurat_obj_list[[i]]@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
}

#check the plot
plot_UMAP=paste("data/pmbc/UMAP.tiff",sep="")
tiff(plot_UMAP,height=1500,width=2000,res=300)
DimPlot(seurat_obj_list[[1]], reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') 
dev.off()








# ##Normalization to HK genes
# 
# hk=read.table("hk.txt")[,1]
# y=FetchData(seurat_obj, vars = hk, slot="counts")
# 
# f=x[hk,]
# x=x[!rownames(x) %in% hk,]
# 
# res=data.frame(gene=NULL, hk=NULL, val=NULL)
# for (i in 1:length(hk)){
#   res=rbind(res, data.frame(gene=rownames(x), hk=rownames(f)[i],val=x[,1]/f$Mean[i]))
# }
