
# nearest neighbor embedding
# param: w
# k=20
# 

library(tidyverse)
library(Seurat)
library(doMC)
registerDoMC(1) #parallel for bootstrapping

out_dir = "coembedding"
meta_dt=read_tsv('thymoma_meta_table.250826.tsv')

gtf2i_pal = c("#3C5488FF", "#E64B35FF", "#00A087FF")
names(gtf2i_pal) = c('m','w','c')
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")


if(F){
  mydend <- read_rds("sequenza.hclust.list.Rds")[[3]]
  plot(mydend)
  rect.hclust(mydend, k = 3, border = 2:7)
  table(cutree(mydend,k=3)) # 2 = wt, 3 = amp, 1 = del
  meta_dt$cn_profile <- NA
  meta_dt$cn_profile[meta_dt$id %in% names(cutree(mydend,k=3))[cutree(mydend,k=3) == 1]] <- "del"
  meta_dt$cn_profile[meta_dt$id %in% names(cutree(mydend,k=3))[cutree(mydend,k=3) == 3]] <- "amp"
  meta_dt$cn_profile[meta_dt$id %in% names(cutree(mydend,k=3))[cutree(mydend,k=3) == 2]] <- "eup"
  boxplot(meta_dt$final_cellularity~meta_dt$cn_profile)
}

merged <- read_rds("merged.Rds")
gene.use <- read_rds("gene.use.Rds")
load("data.scaled.for_comparison.RData")
# merged@reductions$bbknn@assay.used = "RNA"

T0 <- t(thymus_cluster_scaled[thymus_cluster_label %in% c("progenitor", "cTEC", "mTEC", "Tuft", "jTEC"),])
Tx <- t(thymoma_mggene_scaled)

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
colnames(Tx)[1]
meta_dt$GTF2I_status2[meta_dt$id==colnames(Tx)[1]]
thymoma_label <- meta_dt$GTF2I_status2[match(colnames(Tx),meta_dt$id)]
thymoma_histol <- meta_dt$histologic_type[match(colnames(Tx),meta_dt$id)]


mycol <- c(w = "#4dd816", m = "#1c54ac", c = "#f64a1d")[thymoma_label]
mycol <- gtf2i_pal[thymoma_label]

# myhis <- thymoma_histol
mycol2 <- c(A  = "#B60A27",
            AB = "#F4756B",
            B1 = "#C0DDEB",
            B2 = "#6F94C4",
            B3 = "#356199",
            TC = "#936650",
            NE = "#484848",
            "MN-T" = "black")[thymoma_histol]
mycol2 <- histo_pal[thymoma_histol]
mypch <- c(w = 21, m = 23, c = 25)[thymoma_label]

cn.profile <- meta_dt$cn_profile[match(colnames(Tx),meta_dt$id)]
cn.profile[is.na(cn.profile)] <- "NA"
mypch.cn <- c("amp" = 24, "eup" = 23, "del" = 25, "NA" = 21)[cn.profile]
mycol.cn <- c("amp" = "red", "eup" = "#00A0E8", "del" = "#F8B62C", "NA" = "grey")[cn.profile]
cn.pal = c("amp" = "red", "eup" = "#00A0E8", "del" = "#F8B62C", "NA" = "grey")

background_plot <- function(){
  my_color_palette <- c("#F78981", "#CE425A", "#9D0721", "#0BE2A1", "#20A27B", "#00A1FF", "#0B7DC0", "#AB07FF", "#624B92",
                        "#5100FF", "#002EFC", "#1F30BF", "#282C4D", "#1C2362", "#E38900", "#8E766B", "#715757", "#926650",
                        "#BCBCC2", "#84848C", "#74E74C", "#6FA75A", "#102607", "#F766BF")
  # par(mar = c(2,2,2,2),pty="s")
  
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

# draw all points colored by thymoma group
if(F){
  cairo_pdf(paste0(out_dir,"/nearestneighborembedding_orig.pdf"),height = 53.461/25.4,width=53.461/25.4,pointsize = 12*0.7)
  par(mar=rep(0,4),oma=c(0,0,0,0),pty='s')
  
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
  
  background_plot()
  rect(1.5,-10,12.5,1,col="#FFFFFF50")
  box()
  # points(emb, pch = mypch, col = "black", bg=mycol, cex = 1.8)
  points(emb, pch = mypch, col = "black", bg=mycol, cex = 1.8)
  legend("bottomright", legend = c(expression(GTF2I^mut), expression(GTF2I^WT), "TC"),
         pch = c(21,23,25), col = "black", pt.bg=gtf2i_pal,pt.cex=1.7,bty="n",horiz = F,cex=0.9)
  dev.off()
  
}


#draw points colored by histologic types
if(F){
  cairo_pdf(paste0(out_dir,"/supplimentary.nearestneighborembedding.histol.pdf"),height = 10/2.54,width=10/2.54,pointsize = 12*0.7)
  background_plot()
  box()
  rect(1.5,-10,12.5,1,col="#FFFFFF50")
  for(i in 1:ncol(nnw)){nnw[,i] %>%
      merged@reductions$bbknn@cell.embeddings[.,] %>%
      apply(2,median) %>%
      {points(.[1],.[2], pch = mypch[i], col = "black", bg=mycol2[i], cex = 2)}}
  legend("bottomright",  legend = names(histo_pal),
         pch = 23, pt.bg = histo_pal,ncol = 2,pt.cex=1.2, bty = "n")
  dev.off()
}


#try to boot strapping and make confidence interval for coembedding
# Initialize a data frame to store confidence scores
confidence_scores <- data.frame(
  id = colnames(Tx),
  median_x = NA,
  median_y = NA,
  convex_hull_area = NA,
  convex_hull95_area = NA,
  sd_distance_from_median = NA
)
rownames(confidence_scores) <- colnames(Tx)


all_boot_coords <- list()
chull95_coords <- list()
data_dir = paste0(out_dir,'/data')
system(paste0("mkdir -p ", data_dir))
for(i in 1:ncol(nnw)){ # Use nnw (weighted nearest neighbors) consistently
  print(paste0("-----", i, "-----"))
  t_id <- colnames(nnw)[i]
  # Bootstrap process (boot strap genes since gene-gene cor is most important in this result)
  if(file.exists(paste0(data_dir,'/', t_id, ".boot100.rds"))==F){
    x <- foreach(1:100, .combine = rbind) %dopar% { # 100 bootstraps
      boot_genes <- sample(nrow(T0), nrow(T0), replace = TRUE) # Sample genes with replacement
      
      # Calculate weighted correlation for the current thymoma sample against bootstrapped T0
      cor_values <- wcor2(T0[boot_genes,], Tx[boot_genes, i, drop=FALSE], score1[boot_genes])[,1]
      # Identify top 20 nearest neighbors based on this bootstrapped correlation
      top_neighbors <- names(sort(cor_values, decreasing = TRUE))[1:20]
      # Get UMAP coordinates for these neighbors and calculate their median
      merged@reductions$bbknn@cell.embeddings[top_neighbors, ] %>% apply(2, median)
    }
    x %>% saveRDS(paste0(data_dir,'/', t_id, ".boot100.rds"))
  }else{
    x <- readRDS(paste0(data_dir,'/', t_id, ".boot100.rds"))
  }
  
  # Store bootstrap coordinates for this sample
  all_boot_coords[[colnames(nnw)[i]]] <- x
  
  # Calculate confidence scores
  # 1. Median projection (which is already calculated as 'emb' above, but good to re-confirm)
  median_proj <- apply(x, 2, median)
  confidence_scores[colnames(nnw)[i], "median_x"] <- median_proj[1]
  confidence_scores[colnames(nnw)[i], "median_y"] <- median_proj[2]
  
  # 2. Convex Hull Area
  # It's good practice to ensure enough unique points for convex hull calculation (at least 3)
  if (nrow(unique(x)) >= 3) {
    ch_indices <- chull(x)
    # Calculate area of polygon using shoelace formula
    # x and y coordinates of the vertices of the polygon
    ch_coords <- x[ch_indices,]
    area <- 0.5 * abs(sum(ch_coords[-nrow(ch_coords), 1] * ch_coords[-1, 2]) + ch_coords[nrow(ch_coords), 1] * ch_coords[1, 2] -
                        sum(ch_coords[-nrow(ch_coords), 2] * ch_coords[-1, 1]) - ch_coords[nrow(ch_coords), 2] * ch_coords[1, 1])
    confidence_scores[colnames(nnw)[i], "convex_hull_area"] <- area
    
    #save coords for 95% convex hull area (95% confidence area)
    if(T){
      myrank <- ((x[,1] - median(x[,1])) ^ 2 + (x[,2] - median(x[,2])) ^ 2 ) ^ 0.5 %>% rank
      ch95_indices <- chull(x[which(myrank <= 95),]) # Using top 90% as in original script
      ch95 <- x[ch95_indices,]
      chull95_coords[[colnames(nnw)[i]]] <- ch95
      area95 <- 0.5 * abs(sum(ch95[-nrow(ch95), 1] * ch95[-1, 2]) + ch95[nrow(ch95), 1] * ch95[1, 2] -
                          sum(ch95[-nrow(ch95), 2] * ch95[-1, 1]) - ch95[nrow(ch95), 2] * ch95[1, 1])
      confidence_scores[colnames(nnw)[i], "convex_hull95_area"] <- area95
    }
  } else {
    confidence_scores[colnames(nnw)[i], "convex_hull_area"] <- NA # Not enough unique points to form a polygon
    confidence_scores[colnames(nnw)[i], "convex_hull95_area"] <- NA
  }
  
  # 3. Standard Deviation of Distances from Median
  distances <- apply(x, 1, function(row_coords) {
    sqrt(sum((row_coords - median_proj)^2))
  })
  confidence_scores[colnames(nnw)[i], "sd_distance_from_median"] <- sd(distances)
}

#save
print("Saving results")
confidence_scores %>% saveRDS(paste0(out_dir,'/conf_score_df.rds'))
all_boot_coords %>% saveRDS(paste0(out_dir,'/all_boot_coords_list.rds'))
chull95_coords %>% saveRDS(paste0(out_dir,'/chull95_boot_coords_list.rds'))
