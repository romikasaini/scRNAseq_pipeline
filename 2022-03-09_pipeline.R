#Read all the libraries
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(Seurat)
  library(stringr)
})



set.seed(100101) #for UMAP

sample <- c("pmbc", "other_dataset")


for (file in sample) {
  seurat_data <-Read10X(data.dir = paste0("data", "/" , file, "/filtered_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, project = file, min.cells = 3, min.features = 250)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet( seurat_obj, pattern = "^MT-")
  seurat_obj$log10GenesPerUMI <- log10( seurat_obj$nFeature_RNA) / log10( seurat_obj$nCount_RNA)
  seurat_obj= PercentageFeatureSet(seurat_obj, "^RP[SL]", col.name = "percent_ribo")
  seurat_obj= PercentageFeatureSet(seurat_obj, "^HB[^(P)]", col.name = "percent_hb")
  #Preprocessing/filtering
  seurat_obj=subset( seurat_obj, subset = log10GenesPerUMI > 0.8 & nFeature_RNA > 200)
  seurat_obj=subset( seurat_obj, subset = percent.mt < ifelse(median(seurat_obj$percent.mt)*3 < 20,median(seurat_obj$percent.mt)*3, 20))
  seurat_obj=subset(seurat_obj, subset =nFeature_RNA <sort(seurat_obj$nFeature_RNA)[round(length(seurat_obj$nFeature_RNA)*0.95)])
  counts <- GetAssayData(object = seurat_obj, slot = "counts")
  # Output a logical vector for every gene on whether the more than zero counts per cell
  nonzero <- counts > 0
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  keep_genes <- Matrix::rowSums(nonzero) >= as.numeric(summary(Matrix::rowSums(nonzero))[2])
  # Only keeping those genes expressed in more than 10 cells
  filtered_counts <- counts[keep_genes, ]
  # Reassign to filtered Seurat object
  seurat_obj <- CreateSeuratObject(filtered_counts, meta.data = seurat_obj@meta.data, project = "pmbc")
  assign(file, seurat_obj)                                      
}

seurat_obj_list <- list(pmbc = pmbc, other_dataset = other_dataset) #edit here to get all the seurat objects in environment


## Normalization
for (i in 1:length(seurat_obj_list)){
  seurat_obj_list[[i]] <- NormalizeData(seurat_obj_list[[i]], verbose = TRUE, scale.factor = 1000000) #To obtain CPM values scale factor 10^6
}
## Identify highly variable features: features that exhibit high cell-to-cell 
## variation in the data set
for (i in 1:length(seurat_obj_list)){
  seurat_obj_list[[i]] <- FindVariableFeatures(seurat_obj_list[[i]], selection.method = "vst", nfeatures = 2000)
}

#cell cycle scoring
load("data/cycle.rda")
for (i in 1:length(seurat_obj_list)){
  seurat_obj_list[[i]] <-  CellCycleScoring(seurat_obj_list[[i]], g2m.features=g2m_genes, s.features=s_genes)
}

#scaling to UMAP
for (i in 1:length(seurat_obj_list)){
  seurat_obj_list[[i]] <- ScaleData(seurat_obj_list[[i]])
  seurat_obj_list[[i]] <- RunPCA(seurat_obj_list[[i]], npcs = 30, verbose = FALSE) # less than a minute
  seurat_obj_list[[i]] <- RunUMAP(seurat_obj_list[[i]], reduction = "pca", dims = 1:15) #Define here the number of wanted PC's in dim= ! 
  seurat_obj_list[[i]] <- FindNeighbors(seurat_obj_list[[i]], dims = 1:10)
  seurat_obj_list[[i]] <- FindClusters(seurat_obj_list[[i]], resolution = 0.5)
  }

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

#Raw counts and CPM data
for (i in 1:length(seurat_obj_list)){
  ## Calculation MEAN and SD ##
  DefaultAssay(seurat_obj_list[[i]]) <- "RNA"
  gene_list=rownames(seurat_obj_list[[i]])#get detected gene names 
  gene_list=gene_list[1:200]#small test set
  clusts=unique(as.character(seurat_obj_list[[i]]@meta.data[["customclassif"]])) #get cluster names
  #Raw counts 
  counts_genes=list()
  for (j in 1:length(clusts)){
    x=FetchData(seurat_obj_list[[i]],slot="counts",vars=gene_list, cells = colnames(seurat_obj_list[[i]])[seurat_obj_list[[i]]@meta.data[["customclassif"]]==clusts[j]]) #Fetch gene expression for cells
    #included in each cluster
    x=data.frame(Mean=(apply(x, 2, mean)),SD=(apply(x, 2, sd))) 
    counts_genes[[j]]=x
    names(counts_genes)[j]=clusts[j]
    assign(paste0(levels(seurat_obj_list[[i]]$orig.ident), "_count_genes"), counts_genes) 
  }
  #Normalized counts (CPM) 
  CPM_genes=list() #create a list for each cluster
  for (j in 1:length(clusts)){
    x=FetchData(seurat_obj_list[[i]],slot="data",vars=gene_list, cells = colnames(seurat_obj_list[[i]])[seurat_obj_list[[i]]@meta.data[["customclassif"]]==clusts[j]]) #Fetch gene expression for cells
    #included in each cluster
    x=exp(x)-1 #invert log1p normalization and show zeroes correctly. 
    x=data.frame(Mean=(apply(x, 2, mean)),SD=(apply(x, 2, sd))) 
    CPM_genes[[j]]=x
    names(CPM_genes)[j]=clusts[j]
    assign(paste0(levels(seurat_obj_list[[i]]$orig.ident), "_CPM_genes"), CPM_genes)
  }
}




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







