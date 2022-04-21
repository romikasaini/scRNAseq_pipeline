#Read all the libraries
suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(Seurat)
  library(stringr)
  library(graphics)
  library(pastecs)
  library(HGNChelper)
  library(SeuratData)
  library(patchwork)
})

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
load(url("https://raw.githubusercontent.com/romikasaini/scRNAseq_pipeline/main/data/cycle.rda"))
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/romikasaini/scRNAseq_pipeline/main/ScTypeDB_full.xlsx", "Immune system")
InstallData("bmcite") #for multimodal reference mapping
set.seed(100101) #for UMAP

# LOAD FUNCTIONS #
ReadScData <- function(path){
  seurat_data <-Read10X(data.dir = path)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 250)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet( seurat_obj, pattern = "^MT-")
  seurat_obj$log10GenesPerUMI <- log10( seurat_obj$nFeature_RNA) / log10( seurat_obj$nCount_RNA)
  seurat_obj= PercentageFeatureSet(seurat_obj, "^RP[SL]", col.name = "percent_ribo")
  seurat_obj= PercentageFeatureSet(seurat_obj, "^HB[^(P)]", col.name = "percent_hb")
  #seurat_obj <-  CellCycleScoring(seurat_obj, g2m.features=g2m_genes, s.features=s_genes)
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
      png(paste0(plotname, ".png"), width = 960, height = 960, res=200)
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
Cluster=function(seurat_obj, UMAP_dims=15, N_dims=10, res=0.5){
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj<- RunPCA(seurat_obj, npcs = 30, verbose = FALSE) # less than a minute
  seurat_obj<- FindNeighbors(seurat_obj, dims = 1:N_dims)
  seurat_obj<- FindClusters(seurat_obj, resolution = res)
  seurat_obj<- RunUMAP(seurat_obj, reduction = "pca", dims = 1:UMAP_dims) #Define here the number of wanted PC's in dim= ! 
}
FindLog2FC=function(seurat_obj, as.df=F, min.pct=0.25, thres=0.25, idents){
  Idents(seurat_obj)=idents
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min.pct, logfc.threshold = thres)
  if (as.df){
    xx=reshape2::dcast(markers, markers$gene~markers$cluster, value.var = "avg_log2FC")
    rownames(xx)=xx[,1]
    xx[,1]=NULL
    return(list(markers, xx))
  }
  return(markers)
}


clustScore=function(Log2FCdata, gs=gs_list$gs_positive, gs2 = gs_list$gs_negative, gene_names_to_uppercase = !0, ...){
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(Log2FCdata) = toupper(rownames(Log2FCdata));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(Log2FCdata) %in% as.character(gs[[d_]]); rownames(Log2FCdata)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(Log2FCdata) %in% as.character(gs2[[d_]]); rownames(Log2FCdata)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  #if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  Z=Log2FCdata
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = (2**Z[cell_markers_genes_score[jj,"gene_"], ]) * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z, na.rm = T) / sqrt(length(gs_z[!is.na(gs_z)]))); sum_t2 = sum(gz_2, na.rm = T) / sqrt(length(gz_2[!is.na(gz_2)]));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
}
choosCl=function(df){
  res=NULL
  for (i in 1:ncol(df)){
    x=rownames(df)[df[,i]==max(df[,i], na.rm = T) & !is.na(df[,i])]
    if (length(x)>1){
      if (length(grep("B cells", x))==length(x)){
        if (max(df[,i], na.rm = T)<10){
          x="Malignant plasma cells"
        }
        else{
          x="B cells"
        }
        
      }
      if (length(grep("CD8", x))==length(x)){
        x="CD8+ T cells"
      }
      if (length(grep("CD4", x))==length(x)){
        x="CD4+ T cells"
        
      }
      res=c(res, x[1])
    }
    else{
      x=sort(decreasing = T, df[,i])
      res=c(res, names(x)[1])
    }
    
  }
  res[grep("NKT", res)]="NK cells"
  res[grep("CD8", res)]="CD8+ T cells"
  res[grep("CD4", res)]="CD4+ T cells"
  return(res)
}
#Sc-type function
ClustUMAP=function(seurat_obj, matrix, plot=F, save=F, precise=F){
  cl=choosCl(matrix)
  names(cl)=levels(seurat_obj$seurat_clusters)
  seurat_obj@meta.data$customclassif = ""
  for(j in unique(seurat_obj$seurat_clusters)){
    seurat_obj@meta.data$customclassif[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl[as.character(j)])
  }
  if (precise){
    seurat_obj=PreciseCluster(seurat_obj)
  }
  if (plot){
    if (save){
      png("clustUMAP.png", width = 800, height = 600, res = 100)
      plot(DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')) + NoLegend()
      dev.off()
    }
    plot(DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')) + NoLegend()
  }
  
  return(seurat_obj)
}

#multimodal reference mapping
###########
PrepRef=function(){
  bm <- LoadData(ds = "bmcite")
  bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                reduction.key = "wnnUMAP_", return.model = TRUE)
  bm <- ScaleData(bm, assay = 'RNA')
  bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')
  #Computing a cached neighbor index
  bm <- FindNeighbors(object = bm, reduction = "spca", dims = 1:50, graph.name = "spca.annoy.neighbors", k.param = 50, cache.index = TRUE, return.neighbor = TRUE, l2.norm = TRUE )
  return(bm)
}
bm=PrepRef()

Multimodal_UMAP <- function(seurat_obj, ref=bm, plot=F, save=F){
  
  #seuraj_obj is clustered seurat object
  anchors<- FindTransferAnchors(
    reference = bm,
    query = seurat_obj,
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
  Z <- MapQuery(
    anchorset = anchors, 
    query = seurat_obj,
    reference = bm, 
    refdata = list(
      celltype = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  
  if (plot){
    if (save){
      tiff("Multimodal_UMAP.tiff", width = 800, height = 600, res = 100)
      plot(DimPlot(Z, reduction = "ref.umap", group.by =  "predicted.celltype", label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5) + NoLegend())
      dev.off()
    }
    plot(DimPlot(Z, reduction = "ref.umap", group.by =  "predicted.celltype", label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5) + NoLegend())
  }
  return(Z)
}

PreciseCluster=function(seurat_obj){
  obj=seurat_obj
  Idents(obj)="customclassif"
  PC <- subset(obj, idents=grep("T cells|CD4|CD8|NK|DC|Dendritic|B cells",levels(obj), value = T))
  PC=Multimodal_UMAP(PC, ref=bm)
  m=FetchData(PC, "predicted.celltype")
  obj$customclassif[rownames(m)]=m[,1]
  return(obj)
}

#########
CalcRawCount <- function(seurat_obj, clustname="customclassif") {
  DefaultAssay(seurat_obj) <- "RNA"
  gene_list=rownames(seurat_obj)#get detected gene names 
  #gene_list=gene_list[1:200]#small test set
  clusts=unique(as.character(seurat_obj@meta.data[[clustname]])) #get cluster names
  #clusts=clusts[1:3]
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
  #gene_list=gene_list[1:200]#small test set
  clusts=unique(as.character(seurat_obj@meta.data[[clustname]])) #get cluster names
  #clusts=clusts[1:3]
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


SeuratSctype <- function(path, save=T, plot=T, top20=T, precise=T, res=0.8){
  sample=ReadScData(path = path)
  sample=SampleQC(sample)
  z=GetSampleMetrics(sample, plot=plot, save=save) #change here to save
  sample=Cluster(sample, res=res)
  x=FindLog2FC(sample, as.df = T, idents="seurat_clusters")
 
  
  es=clustScore(Log2FCdata = x[[2]])
  sample=ClustUMAP(sample, es, plot = plot, save=save, precise=precise) # select precise = F for 1 step annotation and = T for 2-step annotation
  x=FindLog2FC(sample, as.df = T, idents="customclassif")
  CPM=CalcCPM(sample, clustname = "customclassif")
  Raw=CalcRawCount(sample, clustname =  "customclassif")
  if (top20){
    x[[1]] %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top20
    return(list(Metrics=z, sample=sample, GE=list(CPM=CPM, Raw=Raw), Log2FCperClust=top20))
  }
  return(list(Metrics=z, sample=sample, GE=list(CPM=CPM, Raw=Raw), Log2FCperClust=x[[1]]))
}
                                                                           
                                                                           
                                                                           
