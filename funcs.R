#Read all the libraries
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(Seurat)
  library(stringr)
  library(graphics)
  library(pastecs)
})
load("cycle.rda")
set.seed(100101) #for UMAP
# LOAD FUNCTIONS #
ReadScData <- function(path){
  seurat_data <-Read10X(data.dir = path)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 250)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet( seurat_obj, pattern = "^MT-")
  seurat_obj$log10GenesPerUMI <- log10( seurat_obj$nFeature_RNA) / log10( seurat_obj$nCount_RNA)
  seurat_obj= PercentageFeatureSet(seurat_obj, "^RP[SL]", col.name = "percent_ribo")
  seurat_obj= PercentageFeatureSet(seurat_obj, "^HB[^(P)]", col.name = "percent_hb")
  seurat_obj <-  CellCycleScoring(seurat_obj, g2m.features=g2m_genes, s.features=s_genes)
  return(seurat_obj)
}
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}
CalcPercent=function(vector, min, max){
  x=round((length(which(vector>min & vector<max))/length(vector))*100, 1)
  return(x)
}
CalcTPshare=function(vector){
  d=density(vector)
  ts_y<-ts(d$y)
  tp<-turnpoints(ts_y)
  minima=d$x[tp$tppos][seq(1:length(d$y[tp$tppos]))[!seq(1:length(d$y[tp$tppos])) %in% localMaxima(d$y[tp$tppos])]]
  range=NULL
  res=NULL
  minima=c(0, minima)
  for (i in 1:length(minima)){
    lim=minima[i+1]
    if (is.na(minima[i+1])){
      lim=max(vector)
    }
    range=c(range,paste(round(minima[i],2), round(lim,2), sep="-", collapse = "-"))
    res=c(res,CalcPercent(vector, minima[i], lim))
  }
  return(data.frame(RB_genes_share=range, Pct=res))
}
PlotTP=function(vector){
  d=density(vector)
  ts_y<-ts(d$y)
  tp<-turnpoints(ts_y)
  plot(d, main=deparse(substitute(vector)))
  points(d$x[tp$tppos],d$y[tp$tppos],col="red")
  minima=d$x[tp$tppos][seq(1:length(d$y[tp$tppos]))[!seq(1:length(d$y[tp$tppos])) %in% localMaxima(d$y[tp$tppos])]]
  for (i in 1:length(minima)){
    abline(v=minima[i])
  }
}
GetSampleMetrics <-function(seurat_obj, plot=F, save=F, plotname="stats"){
  cycle=round(table(seurat_obj@meta.data[["Phase"]])/length(seurat_obj@meta.data[["Phase"]])*100, 2)
  rb=seurat_obj@meta.data$percent_ribo
  mt=seurat_obj@meta.data$percent.mt
  hb=seurat_obj@meta.data$percent_hb
  
  if (plot){
    if (save){
      tiff(paste0(plotname, ".tiff"), width = 960, height = 960, res=200)
      par(mfrow=c(3,1))
      PlotTP(rb)
      plot(density(mt), main="mt")
      abline(v=ifelse(median(mt)*3 < 20,median(seurat_obj$percent.mt)*3, 20))
      plot(density(hb), main="hb")
      dev.off()
    }
    par(mfrow=c(3,1))
    PlotTP(rb)
    plot(density(mt), main="mt")
    abline(v=ifelse(median(mt)*3 < 20,median(seurat_obj$percent.mt)*3, 20))
    plot(density(hb), main="hb")
  }
  rb=CalcTPshare(rb)
  mt=round(mean(mt), 2)
  hb=as.numeric(round(table(hb>5)/length(hb)*100, 2)[2])
  z=list(CellCycle.pct=cycle, RibosomeShare=rb, Mean.mt.pct=mt, Erythroid.pct=hb)
  return(z)
}
SampleQC=function(seurat_obj){
  
  #Preprocessing/filtering
  seurat_obj=subset( seurat_obj, subset = log10GenesPerUMI > 0.8 & nFeature_RNA > 200)
  seurat_obj=subset( seurat_obj, subset = percent.mt < ifelse(median(seurat_obj$percent.mt)*3 < 20,median(seurat_obj$percent.mt)*3, 20))
  #remove 5% top outliers
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
  ## Normalization
  seurat_obj <- NormalizeData(seurat_obj, verbose = TRUE, scale.factor = 1000000) #To obtain CPM values scale factor 10^6
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
}
CalcRawCount <- function(seurat_obj, clustname="customclassif") {
  DefaultAssay(seurat_obj) <- "RNA"
  gene_list=rownames(seurat_obj)#get detected gene names 
  gene_list=gene_list[1:200]#small test set
  clusts=unique(as.character(seurat_obj@meta.data[[clustname]])) #get cluster names
  clusts=clusts[1:3]
  #Raw counts 
  counts_genes=data.frame(Mean=NULL, SD=NULL, vars=NULL, clust=NULL)
  for (j in 1:length(clusts)){
    x=FetchData(seurat_obj,slot="counts",vars=gene_list, cells = colnames(seurat_obj)[seurat_obj@meta.data[[clustname]]==clusts[j]]) #Fetch gene expression for cells
    #included in each cluster
    x=data.frame(Mean=(apply(x, 2, mean)),SD=(apply(x, 2, sd)), vars=apply(x, 2, function(y) paste0(y, collapse=";"))) 
    x$clust=clusts[j]
    counts_genes=rbind(counts_genes,x)
  }
  return(counts_genes)
}
CalcCPM<- function(seurat_obj, clustname="customclassif"){
  DefaultAssay(seurat_obj) <- "RNA"
  gene_list=rownames(seurat_obj)#get detected gene names 
  gene_list=gene_list[1:200]#small test set
  clusts=unique(as.character(seurat_obj@meta.data[[clustname]])) #get cluster names
  clusts=clusts[1:3]
  #Normalized counts (CPM) 
  CPM_genes=data.frame(Mean=NULL, SD=NULL, vars=NULL, clust=NULL) #create a list for each cluster
  for (j in 1:length(clusts)){
    x=FetchData(seurat_obj,slot="data",vars=gene_list, cells = colnames(seurat_obj)[seurat_obj@meta.data[[clustname]]==clusts[j]]) #Fetch gene expression for cells
    #included in each cluster
    x=exp(x)-1 #invert log1p normalization and show zeroes correctly. 
    x=data.frame(Mean=(apply(x, 2, mean)),SD=(apply(x, 2, sd)), vars=apply(x, 2, function(y) paste0(y, collapse=";"))) 
    x$clust=clusts[j]
    CPM_genes=rbind(CPM_genes, x)
  }
  return(CPM_genes)
}