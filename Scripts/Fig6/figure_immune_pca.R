library(tidyverse)
library(GSVA)
source("scripts/biplotfun.R")

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='',thickness=0.1,bordercol="black") {
  scale = (length(lut)-1)/(max-min)
  # dev.new(width=1.75, height=5)
  plot(c(min,max),c(0,1), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,yaxs='i')
  # axis(1, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(y,0,y+1/scale,thickness, col=lut[i], border=NA)
  }
  rect(min,0,max,thickness,bg=bordercol)
  # mtext(expression(log[10](CPM~"x"~100)),side = 1,line=2.5)
}

color.bar2 <- function(lut, min, max=-min, ylim=c(min,max), nticks=11, ticks=seq(min, max, len=nticks), title='', thickness=c(0,2)) {
  scale = (length(lut)-1)/(max-min)
  # dev.new(width=1.75, height=5)
  plot(c(0,10),c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,xaxs='i',ylim=ylim)
  # axis(2, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(thickness[1],y,thickness[2],y+1/scale, col=lut[i], border=NA)
  }
  rect(thickness[1],min,thickness[2],max)
  # mtext(expression(log[10](CPM~"x"~100)),side = 1,line=2.5)
}




histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

meta_dt <- read_tsv('thymoma_meta_table.250826.tsv')
exp_dt <- read_tsv('thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')


markers <- list("B cell" = c("BLK","CD19","FCRL2","MS4A1","TNFRSF17","TCL1A","SPIB","PNOC"),
                "Cytotoxic" = c("PRF1","GZMA","GZMB","NKG7","GZMH","KLRK1","KLRB1","KLRD1","CTSW","GNLY"),
                "Dendritic cell" = c("CCL13","CD209","HSD11B1", "NR4A3","HAVCR2","KMO","DNASE1L3","ANPEP","CXCL16"),
                "Exhausted" = c("LAG3","TIGIT","HAVCR2","CTLA4"),
                "Macrophage" = c("CD68","CD163","MS4A4A", "CD84"),
                "Neutrophil" = c("FPR1","SIGLEC5","CSF3R","FCAR","FCGR3B"),
                "NK cell" = c("KIR2DL3","KIR3DL1","KIR3DL2","XCL1","XCL2","NCR1"),
                "proT" = c("ADGRG1","ATP6AP1L","CDK6","CEP70","ETS2","FXYD2","GUCY1A1","GUCY1B1","GXYLT2","HIVEP3","HOXA9","JCHAIN","MEST","NDN","RAB13","RGPD1","RIMS3","RRAS2","TLR7","TNFSF4"),
                "Doublepolar" = c("ACSBG1","AL357060.1","BCL6","C3orf52","CALN1","CD1C","CD1D","CD8A","CD8B","CIB2","CPLX1","CYP2U1","DNMBP","EFNB2","ELOVL4","HRK","LHFPL2","LYST","MCTP1","MIR646HG","RASD1","RIPK4","RMND5A","RORC","SH2D1A","SLAMF1","SLC7A3","TBC1D19"),
                "Singlepositive" = c("CLDN1","CTSL","DUSP4","EGR3","HTR2B","IER3","IFI44L","IRAK2","NR4A2","P2RY1","PDE4D","RSAD2","SERPINE2","SPRY2","TPRG1"),
                "naive" = c("ADTRP","AK5","ATP10A","C1orf162","EPPK1","GIMAP5","GIMAP8","IL6R","INPP4B","LDLRAP1","NOG","PASK","PCED1B","PLEKHA1","RAB30","RPS4Y1","SHISAL2A","TAGAP","UPP1","VSIG1"),
                "Late thymocytes" = c("CD1A","CD1B","DNTT"))
markers[["Thymopoiesis"]] <- c(markers[["proT"]],markers[["Doublepolar"]])
for(i in 1:length(markers)){
  cat(names(markers)[i])
  cat("::")
  cat(paste(markers[[i]],collapse=", "))
  cat("\n")
}
# markers[["Cytotoxic/exhausted"]] <- c(markers[["Exhausted"]],markers[["Cytotoxic"]])
# markers[["TNFa signaling"]] <- c("ABCA1","AC129492.1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")
# markers[["TNFa signaling"]] <- markers[["TNFa signaling"]][markers[["TNFa signaling"]] %in% rownames(l10_exp_dt2)]
# markers[["Neutrophil/macrophage"]] <- c(markers[["Neutrophil"]],markers[["Macrophage"]])

l10_exp_dt2 <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01))) %>% column_to_rownames("gene")
exp_dt2 <- exp_dt %>% column_to_rownames("gene")
if(file.exists("immune_ssgsea_all.Rds")==F){
  df22 <- gsva(as.matrix(l10_exp_dt2), markers[c("Thymopoiesis","Cytotoxic","Exhausted", "Neutrophil","NK cell", "Macrophage","Dendritic cell", "B cell")], method = "ssgsea")
  write_rds(df22,"immune_ssgsea_all.Rds")
}else{
  df22<-read_rds("immune_ssgsea_all.Rds")
}

#compare with previous score. scores should be same except for DC score.
if(F){
  tmp_dt <- df22 %>% t() %>% as.data.frame() %>%
    rownames_to_column("id") %>% as_tibble()
  colnames(tmp_dt) <- paste0("new_", colnames(tmp_dt))
  tmp_meta_dt <- left_join(meta_dt, tmp_dt %>%
              dplyr::rename(id = new_id))
  #sanity check
  tmp_meta_dt %>% filter(round(Thymopoiesis,3) != 
                           round(new_Thymopoiesis, 3))
  tmp_meta_dt %>% filter(round(`Cytotoxic T cell`,3) !=
                           round(`new_Cytotoxic`,3))
  tmp_meta_dt %>% filter(round(`Exhausted T cell`,3) != 
                           round(new_Exhausted,3))
  tmp_meta_dt %>% filter(round(Neutrophil,3) != 
                           round(new_Neutrophil,3))
  tmp_meta_dt %>% filter(round(`NK cell`,3) != 
                           round(`new_NK cell`,3))
  tmp_meta_dt %>% filter(round(Macrophage, 3) != 
                           round(new_Macrophage,3))
  tmp_meta_dt %>% filter(round(`B cell`, 3) != 
                           round(`new_B cell`,3))
  tmp_meta_dt %>% filter(round(`Dendritic cell`, 3) != 
                           round(`new_Dendritic cell`,3))
  #sanity check done 
}
dim(df22) #8 137

pca2 <- prcomp(t(df22))
s2 <- summary(pca2)

if(F){ #excluding NE carcinoma. deprecated
  l10_exp_dt2 <- l10_exp_dt2[,!colnames(l10_exp_dt2) %in% c("TCGA-5U-AB0D", "SNU_09_C")]
  exp_dt2 <- exp_dt2[,!colnames(exp_dt2) %in% c("TCGA-5U-AB0D", "SNU_09_C")]
  # l10_exp_dt3 <- l10_exp_dt2[,meta_dt$id[meta_dt$final_cellularity <= 0.5]]
  # df <- gsva(as.matrix(l10_exp_dt3), markers[c("Thymopoiesis","Cytotoxic","Exhausted", "Neutrophil","NK cell", "Macrophage","Dendritic cell", "B cell")], method = "ssgsea")
  if(F){
    df2 <- gsva(as.matrix(l10_exp_dt2), markers[c("Thymopoiesis","Cytotoxic","Exhausted", "Neutrophil","NK cell", "Macrophage","Dendritic cell", "B cell")], method = "ssgsea")
    write_rds(df2,"immune_ssgsea.Rds")
  }else{
    df2 <- read_rds("immune_ssgsea.Rds")
  }
  
  pca2 <- prcomp(t(df2))
  s2 <- summary(pca2)
}

IGH_gene_ids  <- rownames(l10_exp_dt2)[grep("^IGH",rownames(l10_exp_dt2))]
IGH_sum <- colSums(exp_dt2[IGH_gene_ids,]) %>% {log10(.+0.01)}




"------------------------------------------------------------------------------"
"                                                                              "
"                                  Immune PCA                                  "
"                                                                              "
"------------------------------------------------------------------------------"
cairo_pdf("pca_immune_score.pdf",
          width=(200)/25.4,height = 200/3/25.4,pointsize = 12*0.7)
# s <- svglite::svgstring(width=(200)/25.4,height = 200/3/25.4,pointsize = 12*0.7)
layout(matrix(c(1,1,2,3,4,5,
                1,1,6,7,8,9),byrow=T,ncol=6), widths=c(1,1,1,1,1,1),heights=c(1,1))

par(mar=c(3.5,3,0.5,0.5),oma=c(0,0,0,0))

biplot.prcomp2(pca2, bg=gtf2i_pal[meta_dt$GTF2I_status2[match(colnames(l10_exp_dt2),meta_dt$id)]],pt.cex=2,
               xlim=c(-0.37,0.25),ylim=c(-0.2,0.25),axis.text.cex = 1.5,arrow.len=0.05,
               xlab="",ylab="",
               pch=21,lwd=0.8,
               textaxisbg="transparent",textaxisborder="transparent")
mtext(paste("PCA 1 (", round(s2$importance[2,1]*100, 1), "%)", sep = ""),side=1,line=2,cex = 0.8)
mtext(paste("PCA 2 (", round(s2$importance[2,2]*100, 1), "%)", sep = ""),side=2,line=2,cex = 0.8)
abline(v=0,h=0,lty=2,col="grey40")

# ----------------------------------------------------------------------------
"1. CD1A"
# par(mfrow=c(2,4))
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(-1,log10(300+0.01),length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(l10_exp_dt2["CD1A",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("CD1A",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max = log10(300+0.01),thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(0.99,10,50,300)+0.01),labels = c(0,10,50,300), las=1)
# ----------------------------------------------------------------------------
"2. TDT, PSMB11, LGALS"
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(-1,log10(950+0.01),length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(l10_exp_dt2["DNTT",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("DNTT (TdT)",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max = log10(950+0.01),thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(0.99,10,100,950)+0.01),labels = c(0,10,100,950), las=1)

# ----------------------------------------------------------------------------
"3. CTLA4"
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(-0.5,log10(20+.01),length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(l10_exp_dt2["CTLA4",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("CTLA4",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max = log10(20+.01),thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(0.99,5,20)+0.01),labels = c(0,5,20), las=1)


# ----------------------------------------------------------------------------
'4. LAG3 TIGIT PD1 TIM3 c("LAG3","TIGIT","HAVCR2","CTLA4")'
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(-0.72,log10(43.01),length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(l10_exp_dt2["LAG3",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("LAG3",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max = log10(43.01),thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(0.99,5,15,43)+0.01),labels = c(0,5,15,43), las=1)
# ----------------------------------------------------------------------------
'5. Acute inflammation IFNG'
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(0,6,length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(exp_dt2["IFNG",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("IFNG",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max = 6,thickness = 0.05,bordercol="#00000000")
axis(side=1,at = seq(0,6,length.out=4),labels = seq(0,6,length.out=4), las=1)
# "Neutrophil" = c("FPR1","SIGLEC5","CSF3R","FCAR","FCGR3B"),
# "NK cell" = c("KIR2DL3","KIR3DL1","KIR3DL2","XCL1","XCL2","NCR1"),
# ----------------------------------------------------------------------------
'6. Acute inflammation2 TNFa'
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(-2,log10(1000),length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(l10_exp_dt2["IL6",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("IL6",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max =log10(1000.01),thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(0.99,10,100,1000)+0.01),labels = c(0,10,100,1000), las=1)

# ----------------------------------------------------------------------------
'7. IGH sum'
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(1,4,length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(IGH_sum),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("IGHs",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=1,max = 4,thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(10,100,1000,10000)+0.01),labels = c(10,expression(10^2),expression(10^3),expression(10^4)), las=1)
# ----------------------------------------------------------------------------
'CD20'
par(mar=c(2,.5,.5,.5))
biplot.prcomp2(pca2, 
               bg=circlize::colorRamp2(seq(0,log10(230.01),length.out=6),c("#E9E9E9","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
               )(l10_exp_dt2["MS4A1",]),
               col="black", arrow.len=0.05,lwd=0.6,pt.cex=1.3,bty="n",xaxt='n',yaxt='n',xlab='',ylab='',
               textaxis=F,xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),pch=21,draw.box=F);
mtext("MS4A1 (CD20)",cex=0.8,font=2, side=3,line=-1,adj=c(0.03))
par(new=T,mar=c(2,3,3,3))
color.bar(lut = colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),
          min=0,max = log10(230.01),thickness = 0.05,bordercol="#00000000")
axis(side=1,at = log10(c(0.99,10,50,230)+0.01),labels = c(0,10,50,230), las=1)
# ----------------------------------------------------------------------------
dev.off()
htmltools::browsable(htmltools::HTML(s()))

# ---------------

#"Scatter: Thymopoiesis score ~ purity"
pdf("thymopoiesis_purit_cor.pdf",
          width=8, height=8)
plot(1-meta_dt$final_cellularity,
     df22["Thymopoiesis",match(meta_dt$id,colnames(df22))],
     xlab="",ylab="", pch=21,
     bg=gtf2i_pal[meta_dt$GTF2I_status2],cex=2,lwd=0.6,yaxt='n',
     bty="L")
axis(2,cex.axis=0.9,padj=0.5)
mtext(text = "Thymopoiesis score", side=2,line=2,cex=0.65)
mtext(text = "1-purity", side=1,line=2,cex=0.65)
dev.off()


#"Linear model coefficient plot"
fit <- lm((1-meta_dt$final_cellularity) ~ 
            df22["NK cell",match(meta_dt$id,colnames(df22))] +
            df22["Neutrophil",match(meta_dt$id,colnames(df22))] +
            df22["Dendritic cell",match(meta_dt$id,colnames(df22))] +
            df22["Cytotoxic",match(meta_dt$id,colnames(df22))] +  
            df22["Exhausted",match(meta_dt$id,colnames(df22))] +
            df22["Macrophage",match(meta_dt$id,colnames(df22))] +
            df22["B cell",match(meta_dt$id,colnames(df22))] +
            df22["Thymopoiesis",match(meta_dt$id,colnames(df22))])
summary(fit)
y = coef(fit)
x = seq_along(y)
ci = confint(fit)

pdf("multiple_LR_for_nontcf.pdf",
          width=8, height=4)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
ss <-summary(fit)
sss<- ifelse(ss$coefficients[,4] < 0.001, "***", ifelse(ss$coefficients[,4] < 0.05,"*",""))

plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n", main = "Estimated coefficients for 1-puirity",main.cex=1)
axis(2, at=x, labels=c(
  "(Intercept)",
  "NK cell score",
  "Neutrophil score",
  "Dendritic cell score",
  "Cytotoxic score",
  "Exhausted score",
  "Macrophage score",
  "B cell score",
  "Thymopoiesis score"
), tick=FALSE,las=1)
axis(4, at=x, labels=sss, tick=FALSE,las=2)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
dev.off()