#public scRNAseq tumor cell only
#convert nnemb plotting to ggplot

library(tidyverse)
library(Seurat)
packageVersion("Seurat") #4.1.0
library(SeuratObject)
library(ggsci)
library(ggrastr)
#library(future)
#plan("multisession", workers = 1) 
#plan("sequential")
options(future.globals.maxSize = 5 * 1024^3) 

out_dir ="Seurat_data"
fig_out_dir="pdf_from_R"

#path
#simple integrated object path
sim_obj_path="Seurat_data/sobj_simple_merged.rds"
gene_cl_path="gene_clusters_in_heatmap.tsv"

#color
histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
names(histo_pal) = c("A","AB","MNT","B1","B2","B3","NE","TC")

group_pal = c("#FDE2A9", "#FDBB2B","#3D2F14")
names(group_pal) <- c("g1", "g2", "g3")

cl_pal = c("#f47d32","#407cab","#00A087")
names(cl_pal) <- c("c1","c2","c3")


if(file.exists(paste0(out_dir,'/ff_sobj_tumor_only.rds'))==F){
  #load simple integrated object
  sobj <- readRDS(sim_obj_path)
  dim(sobj) #27631 113948
  
  #filter only thymus tissues exclude peripheral blood
  f_sobj <- subset(sobj, subset = site == "Thymus")
  dim(f_sobj) #27631 83323
  
  #add percent.mt to object
  f_sobj[["percent.mt"]] <- PercentageFeatureSet(f_sobj, pattern = "^MT-")
  
  #check percent.mt
  if(F){
    
    tmp_dt <- f_sobj@meta.data %>% as.data.frame() %>% rownames_to_column("row_name") %>%
      as_tibble()
    ggplot(tmp_dt, aes(x=sample, y=percent.mt))+
      geom_boxplot()
  }
  
  #filter out cells with percent.mt >=10%
  f_sobj <- subset(f_sobj, percent.mt < 10)
  dim(f_sobj) #27631 81542
  
  
  #select only tumor cells
  meta_dt <- f_sobj@meta.data %>% as.data.frame() %>% rownames_to_column("row_name") %>%
    as_tibble()
  meta_dt %>% select(cluster_L1, cluster_L2) %>% unique() %>%
    arrange(cluster_L1) %>% print(n=50) 
  meta_dt %>% filter(cluster_L1 == "TEC") %>%
    group_by(sample) %>% count()
  
  #sample      n
  #<chr>   <int>
  #  1 MG03_TE   116
  #2 MG21_TE    76
  #3 MG22_TE    79
  #4 MG23_TE    18
  #5 T01        45
  #6 T02      3736
  #7 T03       614
  #8 T04       227
  #9 T05       422
  
  #manually add the WHO classification
  meta_dt <- meta_dt %>% mutate(WHO_class = 
                                  case_when(grepl("MG03", sample) ~ "AB",
                                            grepl("MG21", sample) ~ "B1",
                                            grepl("MG22", sample) ~ "B2",
                                            grepl("MG23", sample) ~ "AB",
                                            sample == "T01" ~ "AB",
                                            sample == "T02" ~ "TC",
                                            sample == "T03" ~ "A",
                                            sample == "T04" ~ "AB",
                                            sample == "T05" ~ "MNT",
                                            sample == "T06" ~ "B3"))
  f_sobj@meta.data <- meta_dt %>% as.data.frame() %>% column_to_rownames("row_name")
  
  f_sobj@meta.data$nCount_RNA %>% summary()
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #265    1850    3781    6016    7362   86136 
  
  #subset only tumor cells
  ff_sobj <- subset(f_sobj, cluster_L1 == "TEC" & nCount_RNA >= 1000)
  ncol(ff_sobj) #5287
  
  #routine processing with subset
  ff_sobj <- RunPCA(ff_sobj, verbose = F, npcs = 50)
  ElbowPlot(ff_sobj, ndims = 50, reduction = "pca")
  n_pcs = 20
  ff_sobj <- RunUMAP(ff_sobj, reduction = "pca", dims = 1:n_pcs, verbose = FALSE)
  ff_sobj <- FindNeighbors(ff_sobj, reduction = "pca", dims = 1:n_pcs)
  ff_sobj <- FindClusters(ff_sobj, resolution = 0.1)
  ff_sobj <- FindClusters(ff_sobj, resolution = 0.5)
  ff_sobj <- FindClusters(ff_sobj, resolution = 1)

  
  #500 cutoff already applied       
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #529    6810   11825   13394   17444   73191 
  
  tmp_meta_dt <- ff_sobj@meta.data %>% as.data.frame() %>%
    rownames_to_column("row_name") %>% as_tibble()
  tmp_meta_dt$nCount_RNA %>% summary()
  ggplot(tmp_meta_dt, aes(x=log10(nCount_RNA)))+
    geom_histogram()
  ggplot(tmp_meta_dt, aes(x=percent.mt))+
    geom_histogram()
  
  
  #save
  ff_sobj %>% saveRDS(paste0(out_dir,'/ff_sobj_tumor_only_n1000.rds'))
}else{
  ff_sobj <- readRDS(paste0(out_dir,'/ff_sobj_tumor_only_n1000.rds'))
}

#stats and manual clustering to three
if(T){
  tmp_meta_dt <- ff_sobj@meta.data %>% as.data.frame() %>%
    rownames_to_column("row_name") %>% as_tibble()
  
  tmp_meta_dt %>% group_by(sample) %>% count() 
  #sample      n
  #<chr>   <int>
  #  1 MG03_TE   116
  #2 MG21_TE    60
  #3 MG22_TE    69
  #4 MG23_TE    18
  #5 T01        41
  #6 T02      3731
  #7 T03       608
  #8 T04       223
  #9 T05       421
  tmp_meta_dt %>% group_by(sample) %>% count()  %>% pull(n) %>% summary()
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #18.0    60.0   116.0   587.4   421.0  3731.0 
  
  tmp_meta_dt %>% select(sample, WHO_class) %>% unique() %>%
    group_by(WHO_class) %>% count()
  #WHO_class     n
  #<chr>     <int>
  #  1 A             1
  #2 AB            4
  #3 B1            1
  #4 B2            1
  #5 MNT           1
  #6 TC            1
  
  #draw umap plot
  DimPlot(ff_sobj, group.by="RNA_snn_res.0.1")+
    scale_color_igv()
  #add umap info to meta table
  umap_dt <- ff_sobj@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% rownames_to_column("row_name") %>%
    as_tibble()
  tmp_meta_dt <- left_join(tmp_meta_dt, umap_dt)
  g <- ggplot(tmp_meta_dt, aes(x=UMAP_1, y=UMAP_2))+
    geom_point_rast(aes(color = sample), size=0.5, alpha=0.7,
                    dpi=300)+
    scale_color_igv()+
    guides(color=guide_legend(override.aes = list(size=2, alpha=1)))+
    theme_cowplot()
  pdf(paste0(fig_out_dir,'/scRNAseq_TEC_umap_raster.pdf'),
      width=8, height=6)
  print(g)
  dev.off()
  
  
  #cluster1: 3
  #clsuter2: 2, 6
  #cluster3: 0,1,4,5
  tmp_meta_dt <- tmp_meta_dt %>% 
    mutate(large_group = case_when(RNA_snn_res.0.1 == 3 ~ "c1",
                                   RNA_snn_res.0.1 %in% c(2,6) ~ "c2",
                                   TRUE ~ "c3"))
  
  
  
    
  #proportion bat plot
  tmp_dt1 <- tmp_meta_dt %>% group_by(large_group, sample) %>% count() %>% ungroup()
  tmp_dt2 <- tmp_meta_dt %>% group_by(large_group) %>% 
    summarise(n_total = n()) %>% ungroup()
  tmp_dt <- left_join(tmp_dt1, tmp_dt2) %>%
    mutate(prop = n/n_total)
  ggplot(tmp_dt, aes(x=large_group, y=prop))+
    geom_bar(stat="identity", aes(fill=sample))+
    scale_fill_igv()+
    theme_cowplot()
  
  #group1: MG21, T04 (2 samples)
  #group2: MG03, MG22, MG23, T01, T03, T05 (6samples)
  #group3: T02 (1 sample)
  
  ff_sobj@meta.data <- tmp_meta_dt %>% as.data.frame() %>%
    column_to_rownames("row_name")
  
  tmp_meta_dt
}

#make pseudobulk of the samples and save
if(F){
  ct_mx <- ff_sobj@assays$RNA@counts %>% as.matrix()
  tmp_meta_dt <- ff_sobj@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
    as_tibble()
  sample_list <- tmp_meta_dt$sample %>% unique() %>% sort()
  length(sample_list) #9
  out_list=list()
  for (i in 1:length(sample_list)){
    t_sample = sample_list[i]
    print(t_sample)
    target_cells <- tmp_meta_dt %>% filter(sample == t_sample) %>% pull(cell_id)
    out_list[[i]] <- tibble(gene_name = rownames(ct_mx),
           {{t_sample}}:=ct_mx[,target_cells] %>%rowSums())
    
  }
  lapply(out_list, nrow) # all same row numbers
  pb_dt <- reduce(out_list, full_join, by="gene_name")
  pb_dt %>% write_tsv(paste0(out_dir,'/tumor_pseudobulk_count_table.tsv'))
  
}


#expression comparison with bulkRNAseq
#use cluster1,2,3 genes in Fig. 3a
if(F){
  #load gene cluster information from bulk RNAseq
  gene_cl_dt <- read_tsv(gene_cl_path)
  #filter only cluster1,2,3
  f_cl_dt <- gene_cl_dt %>% filter(cluster %in% c("c1", "c2","c3"))
  
  #load pseudobulk count table
  pb_dt <- read_tsv(paste0(out_dir,'/tumor_pseudobulk_count_table.tsv'))
  #filter only protein coding genes
  bm_dt <- read_tsv('biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt') %>% 
    dplyr::rename(gene_name = `Gene name`, gene_type = `Gene type`) %>% 
    dplyr::select(gene_name, gene_type) %>% unique()
  pb_dt <- left_join(pb_dt, bm_dt)
  f_pb_dt <- pb_dt %>% filter(gene_type == "protein_coding")
  #make it to cpm 
  pb_ct_mx <- f_pb_dt %>% select(-gene_type) %>% 
    as.data.frame() %>% column_to_rownames("gene_name") %>%
    as.matrix()
  dim(pb_ct_mx) #17527 9
  library_sizes <- colSums(pb_ct_mx)
  pb_cpm_mx <- t(t(pb_ct_mx) / library_sizes) * 1e6
  
  #select only intersect genes
  itx_genes <- intersect(f_cl_dt$gene_name, rownames(pb_cpm_mx))
  length(itx_genes) #1905
  f_pb_cpm_mx <- pb_cpm_mx[itx_genes,]
  dim(f_pb_cpm_mx) #1905 9
  #remove all zero genes
  f_pb_cpm_mx <- f_pb_cpm_mx[rowSums(f_pb_cpm_mx) > 0,]
  dim(f_pb_cpm_mx) #1876 9
  #log10
  lf_pb_cpm_mx <- log10(f_pb_cpm_mx + 0.01)
  #draw heatmap 
  ComplexHeatmap::Heatmap(lf_pb_cpm_mx,
                          clustering_distance_columns = 'pearson',
                          clustering_method_columns =  "average")
  #row scaling
  hm_mx <- t(scale(t(lf_pb_cpm_mx))) # row scaling
  #draw heatmap 
  #make annotation table
  tmp_meta_dt <- ff_sobj@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>%
    as_tibble()

  top_annot_dt <- tmp_meta_dt %>% select(sample, WHO_class ) %>% unique() 
  #cluster1: MG21, T04 (2 samples)
  #cluster2: MG03, MG22, MG23, T01, T03, T05 (6samples)
  #cluster3: T02 (1 sample)
  
  top_annot_dt <- top_annot_dt %>%
    mutate(group = case_when(sample %in% c("MG21_TE", "T04") ~ "g1",
                             sample  == "T02" ~ "g3",
                             TRUE ~ "g2"))
  
  top_annot_df <- top_annot_dt %>% as.data.frame() %>% column_to_rownames("sample") 
  top_annot_df <- top_annot_df[colnames(hm_mx),]
  top_anno = ComplexHeatmap::HeatmapAnnotation(df = top_annot_df,
                                               col=list(WHO_class = histo_pal,
                                                        group = group_pal))
  #make row annotation table
  row_annot_df <- left_join(tibble(gene_name = rownames(hm_mx)), f_cl_dt) %>%
    as.data.frame() %>% column_to_rownames("gene_name")
  row_annot_df <- row_annot_df[rownames(hm_mx),,drop=F]
  left_anno <- ComplexHeatmap::rowAnnotation(df = row_annot_df, 
                                             col=list(cluster = cl_pal))
  
  hm <- ComplexHeatmap::Heatmap(hm_mx,
                          clustering_distance_columns = 'pearson',
                          clustering_method_columns =  "average",
                          clustering_distance_rows = "pearson",
                          clustering_method_rows = "complete",
                          show_row_names = F,
                          top_annotation = top_anno,
                          left_annotation = left_anno,
                          use_raster = T)
  
  #Heatmap looks good. I can also try with DEGs but I didn't
  pdf(paste0(fig_out_dir,'/scRNAseq_heatmap_raster.pdf'),
      width=6, height=6)
  print(hm)
  dev.off()

}


#coembedding with mice single cell data
if(F){
  #load pseudobulk count table
  pb_dt <- read_tsv(paste0(out_dir,'/tumor_pseudobulk_count_table.tsv'))
  #filter only protein coding genes
  bm_dt <- read_tsv('biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt') %>% 
    dplyr::rename(gene_name = `Gene name`, gene_type = `Gene type`) %>% 
    dplyr::select(gene_name, gene_type) %>% unique()
  pb_dt <- left_join(pb_dt, bm_dt)
  f_pb_dt <- pb_dt %>% filter(gene_type == "protein_coding")
  #make it to cpm 
  pb_ct_mx <- f_pb_dt %>% select(-gene_type) %>% 
    as.data.frame() %>% column_to_rownames("gene_name") %>%
    as.matrix()
  dim(pb_ct_mx) #17527 9
  library_sizes <- colSums(pb_ct_mx)
  pb_cpm_mx <- t(t(pb_ct_mx) / library_sizes) * 1e6
  
  # load gene.use for coembedding
  gene.use <- read_rds("data/gene.use.Rds")
  
  #load merged single cells
  merged <- read_rds("data/singlecell/merged.Rds")
  
  # prep pseudobulk data --------------------------------------------------------------------------------
  pb_cpm_mggene <- pb_cpm_mx %>% 
    as.data.frame() %>%
    rownames_to_column("hgname") %>% 
    dplyr::left_join(gene.use, by = "hgname") %>%
    na.omit() %>%
    dplyr::select(-hgname) %>% 
    group_by(mgname) %>% 
    summarize_all(sum) %>%
    column_to_rownames("mgname") # %>%
  # > thymoma_tpm_mggene[1:5,1:7]
  #               SNU_03_C SNU_04_C SNU_05_C SNU_06_C SNU_09_C SNU_10_C TCGA-4X-A9FD
  # 0610040J01Rik     0.37     0.93     1.04     1.05     0.22     0.68         0.64
  # 1520401A03Rik     0.04     5.73     1.27     4.14     3.03     3.73         9.24
  pb_mggene_scaled <- scale(log10(t(pb_cpm_mggene) + .1)) %>% {.[is.nan(.)]=0;.}
  
  #load mice single cell data
  load("data/data.scaled.for_comparison.RData")
  
  T0 <- t(thymus_cluster_scaled[thymus_cluster_label %in% c("progenitor", "cTEC", "mTEC", "Tuft", "jTEC"),])
  Tx <- t(pb_mggene_scaled)
  
  #use only intersect genes
  itx_genes <- intersect(rownames(T0), rownames(Tx))
  T0 <- T0[itx_genes,]
  Tx <- Tx[itx_genes,]
  
  if(T){
    # gene weight --------------------------------------------------------------------------------------
    mg_markers <- merged[gene.use$mgname,] %>% subset(simplified %in% c("cTEC", "progenitor","mTEC","jTEC","Tuft")) %>%
      FindAllMarkers(test.use = "wilcox",logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)
    #FindAllMarkers function. Seurat2 log_e_FC -> Seurat3 log2FC
    
    score1 <- mg_markers %>% 
      group_by(cluster) %>%
      mutate(logq = log10(p_val_adj), combined = -logq*avg_log2FC) %>% #avg_logFC -> avg_log2FC
      ungroup() %>%
      dplyr::select(gene, combined) %>%
      group_by(gene) %>% summarize(combined = min(combined)) %>%
      as.data.frame %>% column_to_rownames("gene") %>%
      {.[rownames(Tx),]} %>%
      {.[is.na(.)] <- min(.,na.rm = T);.} %>%
      {.[.==Inf] <- NA;.} %>%
      {.[is.na(.)] <- max(.,na.rm = T);.} %>%
      {.^0.1} %>%
      {(. - min(.))/(max(.)-min(.))}
    hist(score1)
    table(score1>0.2) # using only 1768 -> 1245 genes
    score1 = score1/sum(score1)
    
    # non-vectorized, two column only function
    wmean=function(x,w){sum(x*w)/sum(w)}
    wcov=function(x,y,w){sum(w*(x-wmean(x,w))*(y-wmean(y,w)))/sum(w)}
    wcor=function(x,y,w){wcov(x,y,w)/sqrt(wcov(x,x,w)*wcov(y,y,w))}
    wcor2 <- function(A,B=NULL,w){
      ws=sum(w)
      x=t(A)-colSums(A*w)/ws
      y=t(B)-colSums(B*w)/ws
      r=tcrossprod(sweep(x,2,w,'*')/sqrt(rowSums(sweep(x^2,2,w,'*'))),y/sqrt(rowSums(sweep(y^2,2,w,'*'))))
      r
    }
  }
  # -------------------------------------------------------------------------------------------------
  
  fast_cor <- function(xt,yt=NULL){
    if(is.null(yt)){
      x <- t(xt) - colMeans(xt)
      return(tcrossprod(x / sqrt(rowSums(x ^ 2))))
    } else {
      x <- t(xt) - colMeans(xt)
      y <- t(yt) - colMeans(yt)
      return(tcrossprod(x / sqrt(rowSums(x ^ 2)),y / sqrt(rowSums(y ^ 2))))  
    }
  }
  
  # calculate nn using 2453 non-stromal genes
  # HiClimR::fastCor(A, nSplit = 1,optBLAS=T)
  nn <- fast_cor(T0, Tx) %>% apply(2, function(x){names(sort(x,decreasing = T)[1:20])})
  
  if(T){
    nnw <- wcor2(T0, Tx, score1) %>% apply(2, function(x){names(sort(x,decreasing = T)[1:20])})
  }
  
  background_plot <- function(){
    #my_color_palette <- c("#F78981", "#CE425A", "#9D0721", "#0BE2A1", "#20A27B", "#00A1FF", "#0B7DC0", "#AB07FF", "#624B92",
    #                      "#5100FF", "#002EFC", "#1F30BF", "#282C4D", "#1C2362", "#E38900", "#8E766B", "#715757", "#926650",
    #                      "#BCBCC2", "#84848C", "#74E74C", "#6FA75A", "#102607", "#F766BF")
    # par(mar = c(2,2,2,2),pty="s")
    #levels(Idents(merged))
    #manually curated color palette (25.08.24)
    my_color_palette <- c("#f08381","#c24961", "#90282b", "#57bf99", "#3b9f7f",
                          "#3499d2", "#297db4", "#7252a2", "#635291",  "#4452a4", 
                          "#4e4e3b", "#2f305e", "#2b5aa9", "#314da1", "#E38900",
                          "#8E766B", "#715757", "#926650","#BCBCC2", "#84848C", "#74E74C", "#6FA75A", "#102607", "#F766BF")
                    
    #show_col(my_color_palette)
    plot(merged@reductions$bbknn@cell.embeddings, 
         pch = 20,
         cex=1,
         col = my_color_palette[as.numeric(Idents(merged))],
         # bty = 'l',
         bty='n',
         xaxt='n',yaxt='n',xlab="",ylab="",
         xlim=c(2.5,11.5),ylim=c(-9,-0.57),
         asp = 1)
    
    # mtext(side =1 ,"UMAP (batch-balanced)",cex=0.8)
  }
  

  
  pdf(paste0(fig_out_dir,"/scRNAseq_nn_embedding.pdf"),
      height = 4, width=4,pointsize = 12*0.7)
  #par(mar=rep(0,4),oma=c(0,0,0,0),pty='s')
  
  #option1 - non-weighted correlation
  if(F){
    background_plot()
    rect(1.5,-10,12.5,1,col="#FFFFFF60")
    box()
    emb=matrix(0,ncol=2,nrow=ncol(nn))
    for(i in 1:ncol(nn)){
      position <- nn[,i] %>%
        merged@reductions$bbknn@cell.embeddings[.,] %>%
        apply(2,median)
      emb[i,]=position
      position %>%
        {points(.[1],.[2], pch = mypch[i], col = "black", bg=mycol[i], cex = 2)}}
  }
  
  #option2 weighted pearson
  if(T){
    emb=matrix(0,ncol=2,nrow=ncol(nnw))
    for(i in 1:ncol(nnw)){
      position <- nnw[,i] %>%
        merged@reductions$bbknn@cell.embeddings[.,] %>%
        apply(2,median)
      emb[i,]=position
    }
  }
  
  rownames(emb) = colnames(nnw)
  colnames(emb)=c("x","y")
  
  #color setting
  #group1: MG21, T04 (2 samples)
  #group2: MG03, MG22, MG23, T01, T03, T05 (6samples)
  #group3: T02 (1 sample)
  color_order <- emb %>% as.data.frame() %>% rownames_to_column("row_name") %>%
    as_tibble() %>% mutate(color= 
                             case_when(row_name %in% c("MG21_TE","T04") ~ "g1",
                                       row_name == "T02" ~ "g3",
                                       TRUE ~ "g2")) %>%
    pull(color)
  mycol = group_pal[color_order]
  
  
  background_plot()
  rect(1.5,-10,12.5,1,col="#FFFFFF50")
  box()
  # points(emb, pch = mypch, col = "black", bg=mycol, cex = 1.8)
  
  
  points(emb, pch =25, col = "black", bg=mycol, cex = 1.8)
  text(emb, rownames(emb))
  legend("bottomright", legend = c("g1", "g2", "g3"),
         pch = c(25,25,25), col = "black", pt.bg=group_pal,
         pt.cex=1.7,bty="n",horiz = F,cex=0.9)
  dev.off()
  
  
  
}