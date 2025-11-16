library(tidyverse)
library(Seurat)

merged <- read_rds("merged.Rds")

#stats
if(F){
  DimPlot(merged, reduction = "bbknn", group.by="louvain")
  DimPlot(merged, reduction = "bbknn", group.by="orig.ident")
  
  
  tmp_meta_dt <- merged@meta.data %>% as.data.frame() %>%
    rownames_to_column("row_name") %>% as_tibble()
  nrow(tmp_meta_dt) #15710
  #add umap
  umap_dt <- merged@reductions$bbknn@cell.embeddings %>%
    as.data.frame() %>% rownames_to_column("row_name") %>%
    as_tibble()
  tmp_meta_dt <- left_join(tmp_meta_dt, umap_dt)
  
  #count total cells
  tmp_meta_dt$age %>% table()
  
  tmp_meta_dt$orig.ident %>% table()
  tmp_meta_dt$orig.ident %>% unique() %>% sort()
  #[1] cTEC_early (E12.5-E18.5)         cTEC_mature (E12.5-E18.5)        cTEC (11wo)                     
  #[4] Early endoderm (E8.5)            Upper foregut progenitors (E8.5) mTEC_early (E12.5-E18.5)        
  #[7] mTEC_mature (E12.5-E18.5)        jTEC (2-4wo)                     jTEC_outer (11wo)               
  #[10] jTEC (11wo)                      mTEClo (11wo)                    mTEClo (2-4wo)                  
  #[13] mTEChi (11wo)                    mTEChi (2-4wo)                   Tuft (11wo)                     
  #[16] Pharyngeal Mesoderm (E8.5)       Fibroblast (E12.5-E18.5)         Fibroblast (11wo)               
  #[19] Endothelial (E12.5-E18.5)        Endothelial (11wo)               mono_DC (E12.5-E18.5)           
  #[22] Thymocytes (E12.5-E18.5)         Thymocytes (11wo)                RBC      
  target_cell_types <- tmp_meta_dt$orig.ident %>% unique() %>% sort() %>% .[1:15]
  #filter only TECs 15 cell types
  f_meta_dt <- tmp_meta_dt %>% filter(orig.ident %in% target_cell_types)
  nrow(f_meta_dt) #3481
  
  
  

}


# pdf("figures/singlecellumapwholeshot.pdf", width=300/72, height=400/72)
cairo_pdf("figures/umapwhole.pdf",height = 7.2/2.54,width=12/2.54,pointsize = 12*0.7)
layout(matrix(c(1,2,3,3),ncol=2,byrow=T),widths = c(10,10),heights = c(10,3))
my_color_palette <- c("#F78981", "#CE425A", "#9D0721", "#0BE2A1", "#20A27B", 
                      "#00A1FF", "#0B7DC0", "#AB07FF", "#624B92", "#5100FF", 
                      "#002EFC", "#1F30BF", "#282C4D", "#1C2362", "#E38900", 
                      "#8E766B", "#715757", "#926650", "#BCBCC2", "#84848C", 
                      "#74E74C", "#6FA75A", "#102607", "#F766BF")
par(mar = c(0.5,0.5,0.5,0.5))
plot(merged@reductions$bbknn@cell.embeddings, 
     pch = 20,
     cex=0.4,
     col = my_color_palette[as.numeric(Idents(merged))],
     # bty = 'l',
     xaxt='n',yaxt='n',xlab="",ylab="",
     asp = 1)
# mtext(side =1 ,"UMAP (batch-balanced)",cex=0.8)
rect(1.6,-10.6,12.4,0.56)
# 1.6000000  12.4000000 -10.5664822   0.5664822
plot(merged@reductions$bbknn@cell.embeddings, 
     pch = 20,
     cex=1,
     col = my_color_palette[as.numeric(Idents(merged))],
     # bty = 'l',
     xaxt='n',yaxt='n',xlab="",ylab="",
     xlim=c(1.6,12.4),ylim=c(-10.57,-0.57),
     asp = 1) 
par(mar = c(0, 0, 0, 0))
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legendtxt=levels(Idents(merged))
legendtxt[legendtxt=="jTEC_outer (11wo)"] = "jTEC2 (11wo)"
legend("top",legend=legendtxt, col = my_color_palette,pch=19, ncol = 4,bty="n",cex = 0.7)
dev.off()

epi <- !(merged$simplified %in% c("stroma","etc"))

# # batch figure
# cairo_pdf("figures/umapbatch.pdf",height = 3.3/2.54,width=12/2.54,pointsize = 12*0.7)
batchcolor <- c("#BEBEBE50","#228a3bA0")
# par(mfrow=c(1,4),pty='s',mar=c(0.1,0.1,1.5,0.1))
# for(i in unique(merged$age)){
#   plot(merged@reductions$bbknn@cell.embeddings[epi,], 
#        pch = 21,
#        col="#00000000",
#        cex=1,
#        bg = batchcolor[as.numeric(i == merged$age[epi])+1],
#        bty = 'n',
#        xaxt='n',yaxt='n',xlab="",ylab="",
#        xlim=c(2,12),ylim=c(-9,-1),
#        asp = 1)
#   mtext(i,font=2)
#   
#   points(merged@reductions$bbknn@cell.embeddings[epi,][i == merged$age[epi],], 
#        pch = 21,
#        bg="#228a3bC0",
#        col="#00000000",
#        cex=1)
#        # bg = batchcolor[as.numeric(i == merged$age[epi])+1])
#   mtext(i,font=2)
#   
# }
# dev.off()

# marker figures

# oldmouse %>% FeaturePlot(c("Psmb11", "Aire", "Fezf2", "Sh2d6", "Col1a1", "Cd3e", "Pdpn"), ncol = 3)
markers=c(
  Foxn1="Foxn1",
  Psmb11="Psmb11",
  Lgals7="Lgals7",
  Enpep="Enpep (Ly51)",
  Pax1="Pax1",
  Fezf2="Fezf2", 
  Aire="Aire",
  Cd80="Cd80",
  Sh2d6="Sh2d6",
  L1cam="L1cam",
  Ovol3="Ovol3",
  Prox1="Prox1", 
  Pdpn="Pdpn",
  Tnfrsf11a="Tnfrsf11a (Rank)",
  Ccl21a="Ccl21a",
  Ivl="Ivl (Involucrin)",
  Cldn3="Cldn3",
  Tgfb2="Tgfb2",
  Sall1="Sall1",
  Bmp4="Bmp4",
  Tbx1="Tbx1",
  Irx1="Irx1",
  Irs4="Irs4"
  )#%>%tail(3)
# markers=c("Enpep")
expdat <- GetAssayData(merged)
summary(expdat["Psmb11",epi])
seq(c(0,10))
seq(0,10,length=4)

val2col <- function(x,pal=c("#29394c","#3288bd","#66c2a5","#abdda4","#e6f598",
                            "#ffff33","#fdae61","#f76234","#e84053","#ee165d")){
  # circlize::colorRamp2(seq(0,quantile(x,0.999),length=length(pal)),pal)(x)
  circlize::colorRamp2(seq(0,4,length=length(pal)),pal)(x)
}


cairo_pdf("figures/umapbatchandmarkers.pdf",width=77.5/25.4,height = 120.5/25.4,pointsize = 12*0.7*0.8)
# s <- svglite::svgstring(width=77.5/25.4,height = 132.5/25.4,pointsize = 12*0.7*0.8)
par(mfrow=c(7,4),pty='s',mar=c(0.1,0.1,1.5,0.1))
for(i in unique(merged$age)){
  plot(merged@reductions$bbknn@cell.embeddings[epi,], 
       pch = 21,
       col="#00000000",
       cex=1,
       bg = batchcolor[as.numeric(i == merged$age[epi])+1],
       bty = 'n',
       xaxt='n',yaxt='n',xlab="",ylab="",
       xlim=c(2,12),ylim=c(-9,-1),
       asp = 1)
  # mtext(i,font=2)
  title(i,adj=0,cex=1.1,line=-0.5)
  points(merged@reductions$bbknn@cell.embeddings[epi,][i == merged$age[epi],], 
         pch = 21,
         bg="#228a3bC0",
         col="#00000000",
         cex=1)
  # title(i,adj=0,cex=1.1,line=-0.5)
}



par(mar=c(0.1,0.1,1.5,0.1))
for(i in names(markers)){
  plot(merged@reductions$bbknn@cell.embeddings[epi,][order(expdat[i,epi],decreasing = F),],
       pch = 20,
       cex=0.8,
       col = val2col(expdat[i,epi][order(expdat[i,epi],decreasing = F)]),
       bty = 'n',
       xaxt='n',yaxt='n',xlab="",ylab="",
       xlim=c(2,12),ylim=c(-9,-1),
       asp = 1)
  title(markers[i],adj=0,cex=1.1,line=-0.5)
  # title("Title text", adj = 0.5, line = 0)
}
 
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  # dev.new(width=1.75, height=5)
  plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,yaxs='i')
  axis(1, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(y,0,y+1/scale,2, col=lut[i], border=NA)
  }
  rect(min,0,max,2)
  mtext(expression(log[10]*"("*CPM~"x"~100*")"),side = 1,line=2.5)
}
par(mar=c(5.1,0.5,1.5,0.5),pty='m')
color.bar(colorRampPalette(c("#29394c","#3288bd","#66c2a5","#abdda4","#e6f598",
                             "#ffff33","#fdae61","#f76234","#e84053","#ee165d"))(100),0,4)

dev.off()
htmltools::browsable(htmltools::HTML(s()))

