library(tidyverse)
library(Seurat)
library(gelnet)
library(plot3D)

out_dir="coembedding"

meta_dt <- read_tsv(paste0('thymoma_meta_table.250826.tsv'))
dim(meta_dt)

load(file="data.scaled.for_comparison.RData")
# thymus_cluster_scaled
# thymoma_mggene_scaled
# gene.use.matched
# thymus_cluster_label

# ** OCLC construction ** ----
# construct one-class models to minimize following objective function:
# -1/n * sum(s[i] - log(1+exp(s[i]))) + R(w)
# where s[i] = w'x[i]
# and R(w) = lambda1*sum_j(d[j]|w[j]|) + (lambda2/2)*(w-m)'P(w-m)
set.seed(42) # random seed may not require
OCLC_progenitor <- thymus_cluster_scaled[thymus_cluster_label == "progenitor", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_cTEC <- thymus_cluster_scaled[thymus_cluster_label == "cTEC", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_mTEC <- thymus_cluster_scaled[thymus_cluster_label == "mTEC", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_Tuft <- thymus_cluster_scaled[thymus_cluster_label == "Tuft", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_jTEC <- thymus_cluster_scaled[thymus_cluster_label == "jTEC", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_mTEClo <- thymus_cluster_scaled[thymus_cluster_label == "mTEClo", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_mTEChi <- thymus_cluster_scaled[thymus_cluster_label == "mTEChi", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_jTECi <- thymus_cluster_scaled[thymus_cluster_label == "jTEC_inner", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)
OCLC_mTEClo_old <- thymus_cluster_scaled[thymus_cluster_label == "mTEClo_old", ] %>% gelnet(X = ., y = NULL, l1 = 0, l2 = 0.01)

thymus_cluster_label[1:10]

# small OCLC umap
merged <- read_rds("merged.Rds")
# thymus_cluster_cpm <- GetAssayData(merged_cpm, slot = "data")[gene.use$mgname,] %>%
#   Matrix::t() %>%
#   as.data.frame() %>%
#   {rbind(
#     cbind(y = "progenitor",.[merged_cpm$simplified == "progenitor",]),
#     cbind(y = "cTEC",      .[merged_cpm$simplified == "cTEC",]),
#     cbind(y = "mTEC",      .[merged_cpm$simplified == "mTEC",]),
#     cbind(y = "Tuft",      .[merged_cpm$simplified == "Tuft",]),
#     cbind(y = "jTEC",      .[merged_cpm$simplified2 == "jTEC",]),
#     cbind(y = "mTEClo",    .[merged_cpm$simplified2 == "mTEClo",]),
#     cbind(y = "mTEChi",    .[merged_cpm$simplified2 == "mTEChi",]),
#     cbind(y = "jTEC_inner",.[Idents(merged_cpm) == "jTEC (11wo)",]),
#     cbind(y = "mTEClo_old",.[Idents(merged_cpm) == "mTEClo (11wo)",])
#   )}
my_color_palette <- c("#F78981", "#CE425A", "#9D0721", "#0BE2A1", "#20A27B", "#00A1FF", "#0B7DC0", "#AB07FF", "#624B92",
                      "#5100FF", "#002EFC", "#1F30BF", "#282C4D", "#1C2362", "#E38900", "#8E766B", "#715757", "#926650",
                      "#BCBCC2", "#84848C", "#74E74C", "#6FA75A", "#102607", "#F766BF")
my_color_palette <- c("cTEC"="#c24057","progenitor"="#2a9b77","mTEC"="#203f97","jTEC"="grey","Tuft"="#da8510","etc"="grey","stroma"="grey")
table(merged$simplified)

cairo_pdf(paste0(out_dir,"/oclc.model.used.cells.pdf"),width=40/25.4,height = 45/25.4,pointsize = 12*0.7,onefile=T)
par(mar=rep(1.5,4))
plot(merged@reductions$bbknn@cell.embeddings, 
     pch = 20,
     cex=1,
     col = my_color_palette[merged$simplified],
     # bty = 'l',
     bty='n',
     xaxt='n',yaxt='n',xlab="",ylab="",
     xlim=c(2.5,11.5),ylim=c(-9,-0.57),
     asp = 1)
dev.off()

# ** predict thymoma cases ** ----
# visualization of four indices on tetrahydron space
# scale scores into 0-1 and scale again for the sum of four scores to be 1
# transformation of coord by multiplying a transformation matrix
# visualization with boxplot/heatmap/scatterplots

X2 <- thymoma_mggene_scaled %>% as.matrix
X0 <- thymus_cluster_scaled %>% as.matrix

thymus_scores <- data.frame(cTEC       = X0 %*% OCLC_cTEC$w,
                            mTEC       = X0 %*% OCLC_mTEC$w,
                            Progenitor = X0 %*% OCLC_progenitor$w,
                            Tuft       = X0 %*% OCLC_Tuft$w,
                            jTEC       = X0 %*% OCLC_jTEC$w,
                            mTEClo     = X0 %*% OCLC_mTEClo$w,
                            mTEChi     = X0 %*% OCLC_mTEChi$w,
                            jTECi      = X0 %*% OCLC_jTECi$w,
                            mTEClo_old = X0 %*% OCLC_mTEClo_old$w)

thymoma_scores <- data.frame(cTEC       = X2 %*% OCLC_cTEC$w,
                             mTEC       = X2 %*% OCLC_mTEC$w,
                             Progenitor = X2 %*% OCLC_progenitor$w,
                             Tuft       = X2 %*% OCLC_Tuft$w,
                             jTEC       = X2 %*% OCLC_jTEC$w,
                             mTEClo     = X2 %*% OCLC_mTEClo$w,
                             mTEChi     = X2 %*% OCLC_mTEChi$w,
                             jTECi      = X2 %*% OCLC_jTECi$w,
                             mTEClo_old = X2 %*% OCLC_mTEClo_old$w)
dim(thymoma_scores)

thymoma_scores[,1:4]
dim(meta_dt)
#save
if(F){
  thymoma_scores %>% saveRDS(paste0(out_dir,"/oclc.Rds"))
}

#compare the values with meta table
if(F){
  meta_dt
  tmp_meta_dt <- left_join(thymoma_scores %>% rownames_to_column("id"),
            meta_dt %>% select(id, cTEC, mTEC, Progenitor,
                               Tuft) %>%
              dplyr::rename(cTEC_prev = cTEC,
                            mTEC_prev = mTEC,
                            Progenitor_prev = Progenitor,
                            Tuft_prev = Tuft))
  tmp_meta_dt %>%
    filter(abs(cTEC -  cTEC_prev) >= 0.001) %>%
    select(cTEC, cTEC_prev)
  #no inconsistent row
  tmp_meta_dt %>%
    filter(abs(mTEC -  mTEC_prev) >= 0.001) 
  #no inconsistent row
  tmp_meta_dt %>%
    filter(abs(Progenitor -  Progenitor_prev) >= 0.001) 
  #no inconsistent row
  tmp_meta_dt %>%
    filter(abs(Tuft -  Tuft_prev) >= 0.001) 
  #no inconsistent row
  
  
  
  
}




#previous plotting
if(F){
  # thymoma_mggene_scaled$y <- score_with_metadata$class[match(rownames(thymoma_mggene_scaled),rownames(score_with_metadata))]
  unique(as.factor(thymus_cluster_label))
  tmp_dt <- thymus_scores[thymus_cluster_label %in% c("cTEC","mTEC","jTEC","Tuft","progenitor"),1:5]
  tmp_label <- thymus_cluster_label[thymus_cluster_label %in% c("cTEC","mTEC","jTEC","Tuft","progenitor")]
  tmp_col <- c(progenitor = "#1F9546", cTEC = "#E41A1C",mTEC = "#0E64C6",Tuft = "#E08307",jTEC = "#5E4FA2")
  tmp_col2 <- c(progenitor = "#1F954640", cTEC = "#E41A1C40",mTEC = "#0E64C640",Tuft = "#E0830740",jTEC = "#5E4FA240")
  
  cairo_pdf("figures/supplimentary.oclc.scatter.pdf",height = 15/2.54,width=15/2.54,pointsize = 12*0.7,onefile=T)
  par(oma=c(0,0,0,0))
  plot(tmp_dt,
       bg=tmp_col2[tmp_label],col="#00000000",pch=21)
  par(new=T,mar=c(0,0,0,0),oma=c(0,0,0,0))
  plot(100,pty='n',bty='n',xaxt='n',yaxt='n',xlab="",ylab="",xlim=c(0,1),ylim=0:1)
  legend("bottom",
         legend=names(tmp_col),
         pt.bg=tmp_col,pch=21,ncol = 5,bty='n')
  mtext("Projection of thymus single cell data to one class models",side=3,line=-1,font=2)
  
  thymoma_label <- meta_dt$GTF2I_status2[match(rownames(thymoma_mggene_scaled),meta_dt$id)]
  gtf2i_pal = c(m="#3C5488FF", w="#E64B35FF", c="#00A087FF")
  
  plot(thymoma_scores[,1:5],
       bg = gtf2i_pal[thymoma_label],
       col="black",
       pch=21,cex=1.5,lwd=0.5)
  par(new=T,mar=c(0,0,0,0),oma=c(0,0,0,0))
  plot(100,pty='n',bty='n',xaxt='n',yaxt='n',xlab="",ylab="",xlim=c(0,1),ylim=0:1)
  legend("bottom",
         legend=c(expression(GTF2I^mut),expression(GTF2I^WT),"Thymic carcinoma"),
         pt.bg=gtf2i_pal,pch=21,ncol = 3,bty='n')
  mtext("Projection of thymoma data to one class models",side=3,line=-1,font=2)
  
  
  
  thymoma_label <- meta_dt$GTF2I_status2[match(rownames(thymoma_mggene_scaled),meta_dt$id)]
  gtf2i_pal = c(m="#3C5488FF", w="#E64B35FF", c="#00A087FF")
  
  
  if(F){
    histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
    names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
    thymoma_histol <- meta_dt$histologic_type[match(rownames(thymoma_mggene_scaled),meta_dt$id)]
    plot(thymoma_scores[,1:5],
         bg = histo_pal[thymoma_histol],
         col="black",
         pch=21,cex=1.5,lwd=0.5)
    par(new=T,mar=c(0,0,0,0),oma=c(0,0,0,0))
    plot(100,pty='n',bty='n',xaxt='n',yaxt='n',xlab="",ylab="",xlim=c(0,1),ylim=0:1)
    legend("bottom",
           legend=names(histo_pal),
           pt.bg=histo_pal,pch=21,ncol = 8,bty='n')
    mtext("Projection of thymoma data to one class models",side=3,line=-1,font=2)
    
  }
  dev.off()
  681/72*2.54
  
  cairo_pdf("figures/oclc.box.pdf",width=67.5/25.4,height = 45/25.4,pointsize = 12*0.7*0.7,onefile=T)
  par(mfrow=c(1,4),mar=c(3,2.7,1,0.1),oma=c(1,0,1,0.5),xpd=NA)
  for(i in c("Progenitor","cTEC","mTEC","Tuft")){
    boxplot(thymoma_scores[[i]]~thymoma_label,col=gtf2i_pal[c(3,1,2)],xlab="",xaxt='n',ylab="")
    mtext(paste(i,"index"),3,0.5,cex=1,font=1)
  }
  par(mfrow=c(1,1),new=T,mar=c(0,0,0,0),oma=c(0,0,0,0))
  plot(100,pty='n',bty='n',xaxt='n',yaxt='n',xlab="",ylab="",xlim=c(0,1),ylim=0:1)
  legend("bottom",horiz = T,bty='n',
         legend=c(expression(GTF2I^mut),expression(GTF2I^WT),"Thymic carcinoma"),
         pt.bg=gtf2i_pal,pch=22,pt.cex=2,cex=1.2)
  dev.off()
  
  
  
  cairo_pdf("figures/supplimentary.oclc.box.extended.pdf",height = 17/2.54,width=30/2.54,pointsize = 12*0.7,onefile=T)
  # pdf("figures/fig6d.pdf",7,6)
  par(mfrow=c(2,9),mar=c(0.5,2.7,1,0))
  for(i in c("Progenitor","cTEC","mTEC","Tuft","mTEClo","mTEChi","jTECi","mTEClo_old")){
    boxplot(thymoma_scores[[i]]~thymoma_label,col=gtf2i_pal[c(3,1,2)],xlab="",xaxt='n',ylab="")
    # mtext(paste(i,"index"),1,1,cex=1,font=1)
  }
  plot(100,pty='n',bty='n',xaxt='n',yaxt='n',xlab="",ylab="",xlim=c(0,1),ylim=0:1)
  legend("left",
         legend=c(expression(GTF2I^mut),expression(GTF2I^WT),"Carcinoma"),
         pt.bg=gtf2i_pal,pch=22,pt.cex=2,ncol = 1,bty='n',cex=1.2)
  
  par(mar=c(3,2.7,1,0))
  thymoma_histol_f <- factor(thymoma_histol, levels = c("A","AB","B1","B2","B3","TC","MN-T","NE"))
  for(i in c("Progenitor","cTEC","mTEC","Tuft","mTEClo","mTEChi","jTECi","mTEClo_old")){
    boxplot(thymoma_scores[[i]]~thymoma_histol_f,col=c("#B60A27","#F4756B","#C0DDEB","#6F94C4","#356199","#936650"),xlab="",xaxt='n',ylab="")
    mtext(paste(i,"index"),1,1,cex=1,font=1)
  }
  # par(mar=c(0,0,0,0))
  plot(100,pty='n',bty='n',xaxt='n',yaxt='n',xlab="",ylab="",xlim=c(0,1),ylim=0:1)
  legend("left",
         legend=c("A","AB","B1","B2","B3","TC","MN-T","NE"),
         pt.bg=c("#B60A27","#F4756B","#C0DDEB","#6F94C4","#356199","#936650"),
         pch=22,pt.cex=2,ncol = 1,bty='n',cex=1.2)
  dev.off()
  
  
  
  
  # rm(tmp_dt,tmp_label,tmp_col,tmp_col2,gtf2i_pal)
  
  
  Dim4Dto3D <- function(A){
    as.matrix(A) %*% matrix(c(0,0,0,1,0,0,1/2,sqrt(3)/2,0,1/2,sqrt(3)/6,sqrt(6)/3),byrow=T,ncol=3)
  }
  range01col <- function(X, ...){
    n = ncol(X);  m = nrow(X)
    Min <- matrix(rep(apply(X,2,min, ...),m),byrow=T,ncol=n)
    Max <- matrix(rep(apply(X,2,max, ...),m),byrow=T,ncol=n)
    (X - Min) / (Max - Min)
  }
  rowsumto1 <- function(X, ...){
    X/matrix(rep(rowSums(X),ncol(X)),ncol = ncol(X))
  }
  
  
  
  {
    cairo_pdf("figures/supplimentary.oclc.4dscatterplot.pdf",height = 8/2.54*0.7,width=26/2.54*0.7,pointsize = 12*0.7)
    
    # http://143.248.19.80:8001/graphics/plot_zoom?width=478&height=591&scale=1
    # tmp_dt <- thymoma_scores[thymoma_histol != "NE",c(3,1,2,4)]
    
    ####
    # par(mfrow=c(1,3))
    
    # layout(matrix(c(1,4,5,6,7,
    #                 1,8,9,10,11,
    #                 1,2,2,3,3),ncol=5,byrow=T),width=c())
    
    # layout(cbind(c(1,1),c(2,3)),width = c(3,1.1))
    layout(matrix(c(1,2,3),ncol=3),width = c(3,1,1))
    par(mar = c(0.5,0.5,0.5,0.5))
    line3Dfx <- function(x,y,z,p=pmat,...){
      XY <- trans3D(x,y,z, p = pmat)
      lines(XY$x, XY$y,...)
    }
    text3Dfx <- function(x,y,z,p=pmat,...){
      XY <- trans3D(x,y,z, p = pmat)
      text(XY$x, XY$y,...)
    }
    lollipop3Dfx <- function(x,y,z,p=pmat,pt.col="black",line.col="grey",...){
      XY <- trans3D(x,y,z, p = pmat)
      XYZ <- trans3D(x,y,0, p = pmat)
      segments(x0 = XY$x,y0 = XY$y,x1 = XYZ$x,y1 = XYZ$y,col=line.col,...)
      points(XY$x, XY$y,col=pt.col,...)
    }
    
    par(mar=c(0.5,0.5,0.5,0.5))
    # par(mar=c(0,0,0,0))
    pmat <- thymoma_scores[,c(3,1,2,4)] %>% range01col %>% rowsumto1 %>% Dim4Dto3D() %>%
      {perspbox(.[,1],.[,2],.[,3],
                xlab = "x", ylab = "y", zlab = "z", bty='n',
                # expand = 0.5, d = 2, 
                xlim = c(0,1), ylim = c(0,sqrt(3)/2),zlim=c(0,sqrt(6)/3+0.1),
                phi = 30, theta = 30,
                colvar=NULL)}
    legend("topleft",c("Thymic carcinoma", expression(GTF2I^mut), expression(GTF2I^WT)),
           col = c("#B60927", "#3B4EC1", "#4BAD29"),cex=0.8,bg = "white", pch = 16)
    
    line3Dfx(x = c(0,1,.5,0,.5,.5),
             y = c(0,0,sqrt(3)/2,0,sqrt(3)/6,sqrt(3)/2),
             z = c(0,0,0,0,sqrt(6)/3,0))
    text3Dfx(c(-.04,1.03,.53,.54), c(-.06,0.04,.87,.25),c(0,0,0,0.88),
             labels = c("P", "C", "M", "T"), cex = 0.8)
    thymoma_scores[,c(3,1,2,4)] %>% range01col %>% rowsumto1 %>% Dim4Dto3D() %>%
      {lollipop3Dfx(.[,1],.[,2],.[,3], pt.col='black', line.col="grey",
                    bg=gtf2i_pal[thymoma_label], cex = 1.6, pch = 21)}
    line3Dfx(x = c(1,.5),
             y = c(0,sqrt(3)/6),
             z= c(0,sqrt(6)/3))
    box()
    #ctec mtec tecp tuft
    thymoma_scores[,c(3,1,2,4)] %>% range01col %>% rowsumto1 %>% Dim4Dto3D() %>%
      {plot(.[,2],.[,3],xlim = c(-0.05,sqrt(3)/2+0.05), ylim = c(-0.05,sqrt(6)/3+0.05),
            bg=gtf2i_pal[thymoma_label], asp=1,cex=1.5,
            xaxt='n',yaxt='n',xlab="",ylab="", bty = "n", pch = 21)}
    lines(c(0,sqrt(3)/2,sqrt(3)/6,0), c(0,0,sqrt(6)/3,0))
    text(c(0,0.29,0.85),c(-0.05,0.85,-0.05),c("P+C","T","M"),cex=0.8)
    box()
    thymoma_scores[,c(3,1,2,4)] %>% range01col %>% rowsumto1 %>% Dim4Dto3D() %>%
      {plot(.[,1],.[,2],xlim = c(-0.05,1.05), ylim = c(-0.05,sqrt(3)/2+0.05),
            bg=gtf2i_pal[thymoma_label], asp=1, cex=1.5,
            xaxt='n',yaxt='n',xlab="",ylab="", bty = "n", pch = 21)}
    lines(c(0,1,.5,0), c(0,0,sqrt(3)/2,0))
    rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0)) %>% Dim4Dto3D() %>%
      {text(.[,1]+c(-0.04,0.04,0.04),.[,2]+c(-0.04,-0.04,+0.04),c("P","C","M"),cex=0.8)}
    box()
    dev.off()
  }
  
  
}

