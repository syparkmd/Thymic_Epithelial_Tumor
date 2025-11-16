# data process
library(tidyverse)
library(Seurat)


# 1. basic process of single cell data ---------------------------------------------------------------------------------

if(!any(c("oldmouse","merged") %in% ls())){
  if(T){
    oldmouse <- read_rds("data/singlecell/10x.Rds") #n=1791
  }else{
    CalculateMito <- function(object, plot = T){
      object$percent.mito <- grep(pattern = "^[Mm][Tt]-", rownames(object@assays$RNA@counts), value = TRUE) %>%
      {Matrix::colSums(object@assays$RNA@counts[., ])/Matrix::colSums(object@assays$RNA@counts)}
      if(plot == T) {
        par(mfrow = c(1,2),mar = c(3,2,2,2))
        hist(object$percent.mito, breaks = 100)
        hist(object$nCount_RNA, breaks = 100)
      }
      object
    }
    oldmouse <- Read10X("data/singlecell/10x/outs/filtered_gene_bc_matrices/mm10/") 
    #original cell coujnt 1911
    oldmouse <- Read10X("data/singlecell/10x/outs/filtered_gene_bc_matrices/mm10/") %>%
      CreateSeuratObject() %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors() %>%
      FindClusters()  %>%
      RunTSNE() %>%
      CalculateMito()
    oldmouse$remove <- oldmouse$nCount_RNA > 10^4.3 | oldmouse$nCount_RNA < 10^3.3 | oldmouse$percent.mito > 0.08
    table(oldmouse$remove) #120
    # pdf("figures/not_for_publication_or_supplimentary/10x_doublet_removal.pdf", width= 16, height = 5)
    FeaturePlot(oldmouse, c("nCount_RNA", "percent.mito"), ncol = 3)
    FeaturePlot(oldmouse, features = c("Krt8", "Cd3d"), blend = TRUE)
    FeaturePlot(oldmouse, features = c("nCount_RNA", "remove", "Psmb11"), ncol = 3)
    # dev.off()
    # re-scale after doublet removal
    oldmouse <- subset(oldmouse, subset = remove == FALSE) %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors() %>%
      FindClusters() %>%
      RunTSNE()
    oldmouse@meta.data$remove <- NULL
    # pdf("figures/not_for_publication_or_supplimentary/oldmouse_cluster.pdf", width = 10, height = 10)
    # oldmouse %>% DimPlot() %>% LabelClusters(id = "ident")
    oldmouse %>% DimPlot() %>% LabelClusters(id = "ident")
    oldmouse %>% FeaturePlot(c("Psmb11", "Aire", "Fezf2", "Sh2d6", "Col1a1", "Cd3e", "Pdpn"), ncol = 3)
    oldmouse <- RenameIdents(oldmouse, 
                             `0` = "jTEC (11wo)", 
                             `1` = "Fibroblast (11wo)",
                             `2` = "Endothelial (11wo)", 
                             `3` = "Fibroblast (11wo)",
                             `4` = "mTEChi (11wo)",
                             `5` = "jTEC_outer (11wo)",
                             `6` = "Thymocytes (11wo)",
                             `7` = "mTEClo (11wo)",
                             `8` = "Thymocytes (11wo)", 
                             `9` = "Thymocytes (11wo)", 
                             `10` = "cTEC (11wo)", 
                             `11` = "RBC", 
                             `12` = "Fibroblast (11wo)",
                             `13` = "Tuft (11wo)")
    oldmouse %>% DimPlot() %>% LabelClusters(id = "ident")
    # dev.off()
    # oldmouse %>% write_rds("single_cell/data/10x.Rds", compress = "gz")
    rm(oldmouse,excluded, CalculateMito)
  }
}

# Miragaia data 
# 164 cells SMART seq
if(!any(c("miragia","merged")%in%ls())) {
  if(T){
    miragaia <- read_rds("data/singlecell/miragia.Rds")
    dim(miragaia) #n=164
  }else{
    miragaia0 <- read_tsv("data/singlecell/Miragaia/gene_expression.txt")
    bm <- read_tsv("biomart/biomart_mg.txt", skip = 1, col_names = c("X1","gene"))
    meta <- read_tsv("data/singlecell/Miragaia/cluster_info.txt")
    miragaia2 <- miragaia0 %>% 
      left_join(bm) %>%
      select(-X1) %>%
      group_by(gene) %>%
      summarize_all(funs(sum)) %>%
      as.data.frame() %>%
      na.omit() %>%
      {rownames(.) <- .$gene; .[,-1]}
    miragaia3 <- miragaia2 %>% CreateSeuratObject() %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors() %>%
      RunTSNE()
    DimPlot(miragaia3)
    Idents(miragaia3) <- meta$CClus_final
    DimPlot(miragaia3)
    miragaia <- miragaia3
    miragaia %>% Idents %>% levels
    miragaia <- RenameIdents(miragaia,
                             mTEChi = "mTEChi (2-4wo)",
                             mTEClo = "mTEClo (2-4wo)",
                             jTEC = "jTEC (2-4wo)")
    DimPlot(miragaia) %>% LabelClusters(id = "ident")
    miragaia %>% FeaturePlot(c("Psmb11", "Aire", "Fezf2", "Col1a1", "Cd3e"))
    # write_rds(miragaia, "single_cell/data/miragia.Rds")
    rm(miragaia0, miragaia2,miragaia3,meta,bm)
  }
}

# Kernfeld  (developing thymus)
if(!any(c("Kernfeld","merged") %in% ls())){
  if(T){
    Kernfeld <- read_rds("data/singlecell/kernfeld.Rds")
    dim(Kernfeld) #n = 12861
  }else{
    tmp_fx <- function(x, condition){
      x %>%
        read_tsv(col_types = cols(.default = "i", GENE = "c")) %>%
        group_by(GENE) %>%
        summarize_all(funs(sum)) %>%
        as.data.frame() %>%
        column_to_rownames("GENE") %>%
        CreateSeuratObject() %>%
        {Idents(.) <- condition;.}
    }
    
    Kernfeld <- merge(
      x = tmp_fx("data/singlecell/Kernfeld_et_al/E12_5_wholeThy_venus_1.tsv", "E12.5_PAX9sort"),
      y = list(
        tmp_fx("data/singlecell/Kernfeld_et_al/E13_5_wholeThy_1.tsv", "E13.5_rep1"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E13_5_wholeThy_2.tsv", "E13.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E14_5_wholeThy_1.tsv", "E14.5_rep1"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E14_5_wholeThy_2.tsv", "E14.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E15_5_wholeThy_1.tsv", "E15.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E15_5_wholeThy_2.tsv", "E15.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E16_5_wholeThy_1.tsv", "E16.5_rep1"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E16_5_wholeThy_2.tsv", "E16.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E16_5_wholeThy_3.tsv", "E16.5_rep3"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E17_5_wholeThy_1.tsv", "E17.5_rep1"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E17_5_wholeThy_2.tsv", "E17.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E18_5_wholeThy_1.tsv", "E18.5_rep1"),
        tmp_fx("data/singlecell/Kernfeld_et_al/E18_5_wholeThy_2.tsv", "E18.5_rep2"),
        tmp_fx("data/singlecell/Kernfeld_et_al/P0_wholeThy_1.tsv", "P0_rep1"),
        tmp_fx("data/singlecell/Kernfeld_et_al/P0_wholeThy_2.tsv", "P0_rep2")))
    Kernfeld$batch <- Idents(Kernfeld)
    Kernfeld$date <- stringr::str_replace(Idents(Kernfeld), "_.*","")
    Kernfeld <- Kernfeld %>% 
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA %>%
      FindNeighbors() %>%
      FindClusters(resolution = 3) %>%
      RunTSNE()
    Kernfeld <- RenameIdents(
      Kernfeld,
      `0` = "Thymocytes", `1` = "cTEC1", `2` = "Thymocytes", `3` = "Thymocytes",
      `4` = "Thymocytes", `5` = "cTEC2", `6` = "Thymocytes", `7` = "Thymocytes",
      `8` = "Thymocytes", `9` = "Fibroblast", `10` = "Thymocytes", `11` = "Thymocytes",
      `12` = "Thymocytes", `13` = "Thymocytes", `14` = "Thymocytes", `15` = "cTEC1",
      `16` = "Thymocytes", `17` = "Thymocytes", `18` = "Thymocytes", `19` = "Thymocytes",
      `20` = "Thymocytes", `21` = "cTEC1", `22` = "Thymocytes", `23` = "Thymocytes",
      `24` = "Thymocytes", `25` = "mTEC1", `26` = "Thymocytes", `27` = "mono_DC",
      `28` = "Thymocytes", `29` = "Thymocytes", `30` = "Thymocytes", `31` = "mTEC2",
      `32` = "Fibroblast", `33` = "Thymocytes", `34` = "Endothelia", `35` = "cTEC1")
    
    Kernfeld %>% Idents %>% levels
    Kernfeld <- RenameIdents(Kernfeld,
                             Thymocytes = "Thymocytes (E12.5-E18.5)",
                             mono_DC = "mono_DC (E12.5-E18.5)",
                             cTEC1 = "cTEC_early (E12.5-E18.5)",
                             cTEC2 = "cTEC_mature (E12.5-E18.5)",
                             Fibroblast = "Fibroblast (E12.5-E18.5)",
                             mTEC1 = "mTEC_early (E12.5-E18.5)",
                             mTEC2 = "mTEC_mature (E12.5-E18.5)",
                             Endothelia = "Endothelial (E12.5-E18.5)")
    
    write_rds(Kernfeld, "single_cell/data/kernfeld.Rds", compress = "gz")
    rm(tmp_fx)
    # Kernfeld <- read_rds("single_cell/data/kernfeld.Rds")
    DimPlot(Kernfeld) %>% LabelClusters(id = "ident")
    Kernfeld %>% FeaturePlot(c("Tbx1", "Cdh5"), ncol = 2)
  }
}


#  Ibarra (midgestational embryo)
if(!any(c("Ibarra_final","merged")%in%ls())) {
  if(T){
    Ibarra_final <- read_rds("data/singlecell/ibarra_foregut_and_pharyngeal_mesoderm.Rds")
    dim(Ibarra_final) #n=894
  } else {
    bm <- read_tsv("biomart/biomart_mg.txt", skip = 1, col_names = c("gene", "gene_symbol"))
    Ibarra0 <- read_tsv("data/singlecell/Ibarra/rawCounts.tsv") %>%
      left_join(bm, by = "gene") %>%
      select(-gene) %>%
      group_by(gene_symbol) %>%
      summarize_all(funs(sum)) %>%
      na.omit() %>%
      as.data.frame() %>%
      column_to_rownames("gene_symbol")
    dim(Ibarra0)
    Ibarra0[1:3,1:3]
    "Acta2" %in% rownames(Ibarra0)
    Ibarra <- Ibarra0 %>% CreateSeuratObject() %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA %>%
      FindNeighbors() %>%
      FindClusters()  %>%
      RunTSNE()
    Ibarra_meta <- read_tsv("data/singlecell/Ibarra/cellAnnotation.tsv")
    Ibarra_meta %>% head()
    table(is.na(Ibarra_meta$cellType))
    
    tmp <- Ibarra_meta$cellType
    names(tmp) <- Ibarra_meta$cell
    Idents(Ibarra) <- tmp[names(Idents(Ibarra))]
    Ibarra %>% DimPlot %>% LabelClusters(id = "ident")
    Ibarra %>% FeaturePlot(c("Tbx1", "Foxn1", "Pax9", "Pax8"), ncol = 2)
    # Ibarra %>% write_rds("data/singlecell/ibarra.Rds", compress = "gz")
    rm(bm, Ibarra_meta, tmp)
    Ibarra_foregut <- Ibarra %>%
      subset(idents = "foregut") %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA %>%
      FindNeighbors() %>%
      FindClusters()  %>%
      RunTSNE()
    Ibarra_foregut %>% DimPlot %>% LabelClusters(id = "ident")
    Ibarra_foregut %>% FeaturePlot(c("Tbx1", "Hnf1a", "Pax9", "Pax1"), ncol = 2)
    Ibarra_foregut %>% FeaturePlot(c(
      "Gsc", "Trh", "Otx2",  # early foregut endoderm
      "Ttr", "Tbx3", # hepatic progenitor
      "Pax9", "Tbx1", "Ttf1", "Irx3", "Irx5", "Irx1" # thyroid/lung/thymus
    ))
    Ibarra_foregut %>% VlnPlot(features = c("Tbx1", "Foxn1", "Pax9", "Pax1"), ncol = 2)
    Ibarra_foregut <- Ibarra_foregut %>% RenameIdents(
      `2` = "Hepatic_progenitor",
      `0` = "Lung_thyroid_thymus_progenitor", 
      `1` = "early_endoderm")
    Ibarra_foregut %>% RidgePlot(c("Tbx1", "Gsc","Tbx3"))
    # Ibarra_foregut %>% write_rds("data/singlecell/ibarra_foregut.Rds", compress = "gz")
    Ibarra_final <- merge(
      subset(Ibarra, idents = c("pharyngealMesoderm")),
      subset(Ibarra_foregut, idents = c("early_endoderm", 
                                        "Lung_thyroid_thymus_progenitor")))
    Ibarra_final <- Ibarra_final %>% FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors() %>%
      RunTSNE() 
    
    Ibarra_final %>% Idents %>% levels
    Ibarra_final <- RenameIdents(
      Ibarra_final,
      early_endoderm = "Early endoderm (E8.5)",
      Lung_thyroid_thymus_progenitor = "Upper foregut progenitors (E8.5)",
      pharyngealMesoderm = "Pharyngeal Mesoderm (E8.5)")
    
    # write_rds(Ibarra_final,"data/singlecell/ibarra_foregut_and_pharyngeal_mesoderm.Rds",compress = "gz")
    
    # Ibarra <- read_rds("single_cell/data/ibarra.Rds")
    # Ibarra_final <- read_rds("single_cell/data/ibarra_foregut_and_pharyngeal_mesoderm.Rds")
    Ibarra %>% DimPlot %>% LabelClusters(id = "ident")
    Ibarra_final %>% FeaturePlot(features = "Tbx1")
    cowplot::plot_grid(Ibarra_final %>% FeaturePlot(features = "Tbx1"), Ibarra_final %>% DimPlot)
  }
}

#No. of cells in  public dataset
if(F){
  ncol(miragaia) #n=164
  ncol(Kernfeld) #n = 12861
  ncol(Ibarra_final) #894
  
  ncol(miragaia) + ncol(Kernfeld) + ncol(Ibarra_final) #13919
  
}




# 2. bbknn -------------------------------------------------------------------------------------------------------------
if(!"merged"%in%ls()){
  if(T){
    merged <- read_rds("data/singlecell/merged.Rds")
  }else{
    oldmouse <- read_rds("data/singlecell/10x.Rds")
    miragaia <- read_rds("data/singlecell/miragia.Rds")
    Kernfeld <- read_rds("data/singlecell/kernfeld.Rds")
    Ibarra_final <- read_rds("data/singlecell/ibarra_foregut_and_pharyngeal_mesoderm.Rds")
    
    # marking batch (age of mouse)
    oldmouse$age <- "11wo (this study)"
    miragaia$age <- "2-4wo"
    Kernfeld$age <- "E12.5-E18.5"
    Ibarra_final$age <- "E8.5"
    
    # merging, get pca, basic tsne/umap to compare with bbknn result
    merged <- merge(x = oldmouse, y = list(miragaia, Kernfeld, Ibarra_final)) %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      FindNeighbors() %>%
      RunTSNE() %>%
      RunUMAP(dims = 1:30)
    
    rm(oldmouse, miragaia, Kernfeld, Ibarra_final)
    # write_rds(merged, "single_cell/data/merged.Rds")
    Idents(merged)
    # pdf("figures/not_for_publication_or_supplimentary/merged_single_cell_projections.pdf", 16,11)
    DimPlot(merged, reduction = "tsne", group.by = "ident") %>% LabelClusters(id = "ident")
    DimPlot(merged, reduction = "tsne", group.by = "age") %>% LabelClusters(id = "age")
    # DimPlot(merged, reduction = "umap", group.by = "ident") %>% LabelClusters(id = "ident")
    # DimPlot(merged, reduction = "umap", group.by = "age") %>% LabelClusters(id = "age")
    # dev.off()
    
    
    # RunBBKNN (ugly version)-----------------------------------------------------
    RunBBKNN <- function(object,   
                         batch.key = "batch", # slot of object that store batch information
                         dims.use = 1:50,
                         neighbors_within_batch = 6,
                         trim = 15,
                         reduction.use = "pca",
                         python.path = "/home/users/kjyi/anaconda3/bin/python") {
      write_tsv(as.data.frame(object@reductions[[reduction.use]]@cell.embeddings[,dims.use]),"pca.tmp.tsv")
      write_tsv(data.frame(object[[batch.key]]), "batch.tmp.tsv")
      write_lines(paste0("
                         import numpy as np
                         import pandas as pd
                         import scanpy.api as sc
                         import anndata
                         import bbknn
                         import os
                         from scipy import sparse
                         X=pd.read_table('pca.tmp.tsv')
                         obs=pd.read_table('batch.tmp.tsv')
                         adata=anndata.AnnData(X=X.as_matrix(), obs=obs)
                         sc.tl.pca(adata)
                         adata.obsm.X_pca = X.as_matrix()[:,0:(adata.obsm.X_pca[1].size)]
                         bbknn.bbknn(adata, batch_key = '",batch.key,"', neighbors_within_batch = ",neighbors_within_batch,", trim = ",trim,")
                         sc.tl.umap(adata)
                         sc.tl.louvain(adata, resolution=2)
                         np.savetxt('umap.tmp.tsv',adata.obsm.X_umap,delimiter='\t')
                         adata.obs.to_csv('louvain.tmp.csv')
                         "),"bbknn.tmp.py")
      pythonexitstatus <- system(paste(python.path, "bbknn.tmp.py"))
      if(pythonexitstatus != 0) {stop()}
      umap <- read_tsv("umap.tmp.tsv",col_names=c("BBKNN_1","BBKNN_2")) %>%
        as.data.frame %>% as.matrix
      louvain <- read_csv("louvain.tmp.csv")$louvain
      rownames(umap) <- rownames(object@reductions[[reduction.use]]@cell.embeddings)
      object[["bbknn"]] <- CreateDimReducObject(embeddings = umap, assay = "bbknn", key = "BBKNN_")
      object$louvain <- read_csv("louvain.tmp.csv")$louvain
      system("rm -f louvain.tmp.csv umap.tmp.tsv pca.tmp.tsv batch.tmp.tsv bbknn.tmp.py")
      object
    }
    
    # current best parameter ----
    merged <- RunBBKNN(merged,batch.key = "age", dims.use = 1:50, neighbors_within_batch = 6, trim = 15)
    
    # ordering group
    Idents(merged) <- factor(Idents(merged), levels = c("cTEC_early (E12.5-E18.5)",
                                                        "cTEC_mature (E12.5-E18.5)",
                                                        "cTEC (11wo)",
                                                        "Early endoderm (E8.5)",
                                                        "Upper foregut progenitors (E8.5)",
                                                        "mTEC_early (E12.5-E18.5)",
                                                        "mTEC_mature (E12.5-E18.5)",
                                                        "jTEC (2-4wo)",
                                                        "jTEC_outer (11wo)",
                                                        "jTEC (11wo)",
                                                        "mTEClo (11wo)",
                                                        "mTEClo (2-4wo)",
                                                        "mTEChi (11wo)",
                                                        "mTEChi (2-4wo)",
                                                        "Tuft (11wo)",
                                                        "Pharyngeal Mesoderm (E8.5)",
                                                        "Fibroblast (E12.5-E18.5)",
                                                        "Fibroblast (11wo)",
                                                        "Endothelial (E12.5-E18.5)",
                                                        "Endothelial (11wo)",
                                                        "mono_DC (E12.5-E18.5)",
                                                        "Thymocytes (E12.5-E18.5)",
                                                        "Thymocytes (11wo)",
                                                        "RBC"))
    rm(RunBBKNN)
    
    # write_rds(merged, "data/singlecell/merged.Rds") # save point
  }
}

merged <- read_rds("data/singlecell/merged.Rds") # save point
table(merged$simplified,merged$age)
