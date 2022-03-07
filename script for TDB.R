## QC ##
OT.combined=Read10X(data.dir = paste0("./Samples/",ls[i],"/filtered_feature_bc_matrix/"))
names(OT.combined)=ls[i]
OT.combined <- CreateSeuratObject(counts = OT.combined, project = ls[i], min.cells = 3, min.features = 250)
OT.combined[["percent.mt"]] <- PercentageFeatureSet( OT.combined, pattern = "^MT-")
OT.combined$log10GenesPerUMI <- log10( OT.combined$nFeature_RNA) / log10( OT.combined$nCount_RNA)
OT.combined= PercentageFeatureSet(OT.combined, "^RP[SL]", col.name = "percent_ribo")
OT.combined= PercentageFeatureSet(OT.combined, "^HB[^(P)]", col.name = "percent_hb")
OT.combined=subset(OT.combined, subset =
                          log10GenesPerUMI > 0.8 & nFeature_RNA >200 &
                          percent.mt<median(OT.combined$percent.mt)*3)
# pbmc=OT.combined
# pbmc=as.SingleCellExperiment(pbmc)
# pbmc <- scDblFinder(pbmc)
# OT.combined=as.Seurat(pbmc, project = ls[i])

OT.combined=subset(OT.combined, subset =nFeature_RNA <sort(OT.combined$nFeature_RNA)[round(length(OT.combined$nFeature_RNA)*0.95)])
counts <- GetAssayData(object = OT.combined, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= as.numeric(summary(Matrix::rowSums(nonzero))[2])
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
OT.combined <- CreateSeuratObject(filtered_counts, meta.data = OT.combined@meta.data, project = ls[i])
OT.combined<- NormalizeData(OT.combined, verbose = TRUE, scale.factor = 1000000) #To obtain CPM values scale factor 10^6
OT.combined<- CellCycleScoring(OT.combined, g2m.features=g2m_genes, s.features=s_genes)

## Calculation MEAN and SD ##
DefaultAssay(OT.combined) <- "RNA"
gene_list=rownames(OT.combined)#get detected gene names 
gene_list=gene_list[1:200]#small test set
clusts=unique(as.character(OT.combined@meta.data[["manid08"]])) #get cluster names

#Raw counts 
counts_genes=list()
for (i in 1:counts_genes){
  x=FetchData(OT.combined,slot="counts",vars=gene_list, cells = colnames(OT.combined)[OT.combined@meta.data[["manid08"]]==clusts[i]]) #Fetch gene expression for cells
                                                                                                                                      #included in each cluster
  x=data.frame(Mean=(apply(x, 2, mean)),SD=(apply(x, 2, sd))) 
  counts_genes[[i]]=x
  names(counts_genes)[i]=clusts[i]
}

# ##Normalization to HK genes
# 
# hk=read.table("hk.txt")[,1]
# y=FetchData(OT.combined, vars = hk, slot="counts")
# 
# f=x[hk,]
# x=x[!rownames(x) %in% hk,]
# 
# res=data.frame(gene=NULL, hk=NULL, val=NULL)
# for (i in 1:length(hk)){
#   res=rbind(res, data.frame(gene=rownames(x), hk=rownames(f)[i],val=x[,1]/f$Mean[i]))
# }

#Normalized counts (CPM) 
CPM_genes=list() #create a list for each cluster
for (i in 1:length(CPM_genes)){
  x=FetchData(OT.combined,slot="data",vars=gene_list, cells = colnames(OT.combined)[OT.combined@meta.data[["manid08"]]==clusts[i]]) #Fetch gene expression for cells
                                                                                                                                    #included in each cluster
  x=exp(x)-1 #invert log1p normalization and show zeroes correctly. 
  x=data.frame(Mean=(apply(x, 2, mean)),SD=(apply(x, 2, sd))) 
  CPM_genes[[i]]=x
  names(CPM_genes)[i]=clusts[i]
}



