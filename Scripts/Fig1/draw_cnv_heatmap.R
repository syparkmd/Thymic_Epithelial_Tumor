# Copied from KJYI script
#
#                               Sequenza heatmap
#
# ______________________________________________________________________________

# prep ----
library(tidyverse)
meta_dt <- read_tsv('thymoma_meta_table.250826.tsv')
fl <- system("ls sequenza_wes/rerun/Final*/*/*segments.txt",
             intern = T)
fl.info <- data.frame(file = fl,
                      id = str_extract(fl,"(?<=/)[^/]*(?=.v)"),
                      ploidy = str_extract(fl,"(?<=_p)[^_]*"),
                      purity = str_extract(fl,"(?<=_c)[^_]*"),
                      stringsAsFactors = F)
seg <- lapply(fl, read_tsv, col_types = "ciidddddddddd")
for(i in 1:length(seg)){
  seg[[i]]$case = fl.info$id[i]
}
seg.df <- do.call(rbind, seg)
seg.bound <-group_by(seg.df, chromosome) %>% 
  summarize(start = min(start.pos), end = max(end.pos))
seg.bound$chromosome <- factor(seg.bound$chromosome, 
                               levels = c("1","2","3","4","5","6","7","8","9",
                                          "10","11","12","13","14","15","16",
                                          "17","18","19","20","21","22","X","Y"))
seg.bound <- arrange(seg.bound, chromosome)
seg.bound$size = seg.bound$end - seg.bound$start 
seg.bound$new.start = c(1,cumsum(as.numeric(seg.bound$size))[1:(nrow(seg.bound)-1)])
seg.bound$new.end = cumsum(as.numeric(seg.bound$size))
seg.bound$transform_gap = seg.bound$new.start - seg.bound$start

seg.df$new.start <- seg.df$start.pos + structure(seg.bound$transform_gap, names=as.character(seg.bound$chromosome))[seg.df$chromosome]
seg.df$new.end   <- seg.df$end.pos + structure(seg.bound$transform_gap, names=as.character(seg.bound$chromosome))[seg.df$chromosome]

# filter
seg.df <- seg.df[(seg.df$end.pos - seg.df$start.pos) >= 3000000,]

unique(seg.df$case[is.na(seg.df$A)])
# [1] "SNU_03_C"     "SNU_06_C"     "SNU_09_C"     "SNU_10_C"     "SNU_19_C"     "SNU_22_C"     "SNU_26_C"     "TCGA-3G-AB14" "TCGA-3S-A8YW"
# [10] "TCGA-4V-A9QJ" "TCGA-4V-A9QM" "TCGA-4V-A9QW" "TCGA-4X-A9FA" "TCGA-5K-AAAP" "TCGA-5U-AB0E" "TCGA-5V-A9RR" "TCGA-X7-A8D8" "TCGA-X7-A8DF"
# [19] "TCGA-X7-A8M1" "TCGA-XH-A853" "TCGA-XM-A8RC" "TCGA-XM-A8RG" "TCGA-XM-A8RI" "TCGA-XM-AAZ1" "TCGA-XM-AAZ2" "TCGA-XM-AAZ3" "TCGA-XU-A92O"
# [28] "TCGA-XU-A92Q" "TCGA-XU-A92U" "TCGA-XU-A92X" "TCGA-XU-A92Z" "TCGA-XU-A930" "TCGA-YT-A95E" "TCGA-ZB-A961" "TCGA-ZB-A962" "TCGA-ZB-A963"
# [37] "TCGA-ZB-A964" "TCGA-ZB-A96E" "TCGA-ZB-A96F" "TCGA-ZB-A96H" "TCGA-ZB-A96K" "TCGA-ZB-A96Q" "TCGA-ZB-A96V" "TCGA-ZC-AAA7" "TCGA-ZL-A9V6"

# fill A,B in male X,Y chromosome (as CNt == A == B)
seg.df$A[is.na(seg.df$A)] <- seg.df$CNt[is.na(seg.df$A)]
seg.df$B[is.na(seg.df$B)] <- seg.df$CNt[is.na(seg.df$B)]


# case order??
mat.A = matrix(0, nrow = length(unique(seg.df$case)), ncol = 3000)
mat.B = matrix(0, nrow = length(unique(seg.df$case)), ncol = 3000)
rownames(mat.A) = unique(seg.df$case)
rownames(mat.B) = unique(seg.df$case) 

for(i in 1:3000){
  print(i)
  bin.start = (i-1)*1000000+1
  bin.end = (i)*1000000
  ovl = pmax(pmin(seg.df$new.end,bin.end)-pmax(seg.df$new.start,bin.start),0)
  for(j in 1:nrow(mat.A)){
    ovl.c = ovl * (seg.df$case == rownames(mat.A)[j])
    avg.A = mean(seg.df$A * ovl.c / sum(ovl.c), na.rm=T)
    avg.B = mean(seg.df$B * ovl.c / sum(ovl.c), na.rm=T)
    mat.A[rownames(mat.A)[j],i] = avg.A
    mat.B[rownames(mat.B)[j],i] = avg.B
  }
}

mat.A[1:10,1:30]
mat.A.bak = mat.A
mat.B.bak = mat.B
if(F){
  write_rds(mat.A, "data/mat.A.Rds", compress= "gz")
  write_rds(mat.B, "data/mat.B.Rds", compress= "gz")
}else{
  mat.A <- read_rds("data/mat.A.Rds")
  mat.B <- read_rds("data/mat.B.Rds")
}



mat.A[is.nan(mat.A)] = 0
mat.B[is.nan(mat.B)] = 0
mat.A = mat.A[,1:2798] # exclude sex chromosome
mat.B = mat.B[,1:2798]
image(mat.A[,])

mat.A.t = t(mat.A)
mat.B.t = t(mat.B)

dist.A = dist(mat.A)
dist.B = dist(mat.B)
dist.AB = dist.A + dist.B

# cor.dist.A = as.dist(1-cor(mat.A.t))
# cor.dist.B = as.dist(1-cor(mat.B.t))
# cor.dist.AB = cor.dist.A + cor.dist.B

hclust.AB.wardD2 = hclust(dist.AB,method = "ward.D2")
hclust.AB.ward   = hclust(dist.AB,method = "ward.D")

# hclust.corAB.wardD2 = hclust(cor.dist.AB,method = "ward.D2")
# hclust.corAB.ward   = hclust(cor.dist.AB,method = "ward.D")


# dist.mk0.5.A = dist(mat.A, method="minkowski", p=0.5)
# dist.mk1.5.A = dist(mat.A, method="minkowski", p=1.5)
# dist.mk2.0.A = dist(mat.A, method="minkowski", p=2)
# dist.mk3.0.A = dist(mat.A, method="minkowski", p=3)
# dist.mk0.5.B = dist(mat.B, method="minkowski", p=0.5)
# dist.mk1.5.B = dist(mat.B, method="minkowski", p=1.5)
# dist.mk2.0.B = dist(mat.B, method="minkowski", p=2)
# dist.mk3.0.B = dist(mat.B, method="minkowski", p=3)
# 
# 
# hclust.mk0.5.AB.wardD = hclust(dist.mk0.5.A+dist.mk0.5.B, method = "ward.D")
# hclust.mk1.5.AB.wardD = hclust(dist.mk1.5.A+dist.mk1.5.B, method = "ward.D")
# hclust.mk2.0.AB.wardD = hclust(dist.mk2.0.A+dist.mk2.0.B, method = "ward.D")
# hclust.mk3.0.AB.wardD = hclust(dist.mk3.0.A+dist.mk3.0.B, method = "ward.D")
# hclust.mk3.0.AB.wardD = hclust(dist.mk3.0.A+dist.mk3.0.B, method = "ward.D")
# hclust.mk3.0.AB.wardD = hclust(dist.mk3.0.A, method = "ward.D")

plot(hclust.mk2.0.AB.wardD)
rect.hclust(hclust.mk2.0.AB.wardD, k = 3, border = 2:7)

plot(hclust.AB.ward)
rect.hclust(hclust.AB.ward, k = 3, border = 2:7)
plot(hclust.AB.wardD2)
rect.hclust(hclust.AB.wardD2, k = 3, border = 2:7)

library(seriation)
hclust.AB.HC_ward = seriate(dist.AB, method = "HC_ward")
hclust.AB.GW_ward = seriate(dist.AB, method = "GW_ward")
hclust.AB.OLO_ward = seriate(dist.AB, method = "OLO_ward")
hclust.AB.HC_ward = seriate(dist.AB, method = "DendSer", control = list(method = "ward.D2", criterion = "AR"))
plot(hclust.AB.HC_ward[[1]])
plot(hclust.AB.GW_ward[[1]])
plot(hclust.AB.OLO_ward[[1]])

dendextend::click_rotate(as.dendrogram(hclust.AB.ward))
dendextend::rotate(as.dendrogram(hclust.AB.ward),order=3)
mydend = dendextend::rotate_DendSer(as.dendrogram(hclust.AB.ward))
plot(mydend)
rect.hclust(as.hclust(mydend), k = 3, border = 2:7)

# case_order = structure(1:length(unique(seg.df$case)), names = unique(seg.df$case))
case_order = structure(sort(hclust.AB$order), names = hclust.AB$labels[hclust.AB$order])
# hclust.AB$labels[hclust.AB$order]

# plot ----

pdf("figures/sequenza.heatmap.pdf",195/25.4,93.427/25.4)
for(II in 3){
  hclust.obj = list(hclust.AB.wardD2,
                    hclust.AB.ward,
                    # hclust.AB.HC_ward[[1]],
                    as.hclust(mydend)
                    )[[II]]
  hclust.title = c("euclidean.A+B.wardD2",
                   "euclidean.A+B.ward",
                   # "hclust.AB.HC_ward",
                   "euclidean.A+B.ward.DendSer")[II]
  case_order = structure(sort(hclust.obj$order), names = hclust.obj$labels[hclust.obj$order])
  # group
  colors_ordered =c("m"="#334A77","c"="#268B77","w"="#C44834")[structure(meta_dt$GTF2I_status2, names=meta_dt$id)[names(case_order)]]
  
  # histology
  histo_pal <- c("MN-T"="#564E8B", 
                 "A"="#57246C",
                 "AB"="#334382", 
                 "B1"="#911D51",
                 "B2"="#A72124",
                 "B3"="#D01A1B", 
                 "CA-SqCC"="#137977", 
                 "CA-UN"="black",
                 "NE"="#8AB435")
  colors_ordered3 =histo_pal[structure(meta_dt$histologic_type2, names=meta_dt$id)[names(case_order)]]
  
  # purity
  mypal= circlize::colorRamp2(0:8/8,RColorBrewer::brewer.pal(9,"PuBu"))
  colors_ordered2=mypal(structure(meta_dt$final_cellularity, names=meta_dt$id)[names(case_order)])
  par(mar=c(5,5,5,10))
  
  # par(mfrow=c(1,2))
  layout(matrix(c(1:3),ncol=3), widths=c(1,5,1))
  par(mar=c(2,2,2,0))
  plot(as.dendrogram(hclust.obj),horiz=T,yaxs='i',yaxt="n",xaxt='n',leaflab = "none")
  par(mar=c(2,0,2,1))
  
  plot(1, type = 'n', 
       xlim = c(1,max(seg.bound$new.end)),
       ylim = c(0,length(unique(seg.df$case))),xaxs='i',yaxs='i',
       # main = hclust.title, yaxt="n", xaxt="n")
  main = "", yaxt="n", xaxt="n")
  
  rect(xleft = seg.df$new.start,
       xright = seg.df$new.end,
       ybottom = case_order[seg.df$case]-1,
       ytop = case_order[seg.df$case],
       col = ifelse(seg.df$A > 2, "#FF0000FF", 
                    ifelse(seg.df$A == 2, "#FF000080", "#00000000")),
       border = F)
  rect(xleft = seg.df$new.start,
       xright = seg.df$new.end,
       ybottom = case_order[seg.df$case]-1,
       ytop = case_order[seg.df$case],
       col = ifelse(seg.df$B < 1, "#0000FF80", "#00000000"),
       border = F)

  abline(v=seg.bound$new.start[-1], col="black")
  invisible(lapply(1:nrow(seg.bound),function(i){
    axis(side=3, at = (seg.bound$new.end+seg.bound$new.start)[i]*0.5,
         labels = seg.bound$chromosome[i],tick = F)
  }))
  axis(side=3,at = (seg.bound$new.end+seg.bound$new.start)*0.5, labels = seg.bound$chromosome,tick = F)
  axis(side=3,at = c(1,seg.bound$new.start), labels = F,tick = T)
  
  axis(side=4,at = c(0,case_order),labels = F,las=2, tick=T)

  # invisible(lapply(1:length(colors_ordered),function(i){
  #  axis(side=4,at = case_order[i]-0.5,labels = names(case_order)[i],las=2, tick=F,
  #       cex.axis=0.5,col.axis = colors_ordered[i])
  # }))
  par(mar=c(2,0,2,2))
  plot(1, type = 'n', 
       xlim = c(1,3),
       ylim = c(0,length(unique(seg.df$case))),xaxs='i',yaxs='i',
       # main = hclust.title, yaxt="n", xaxt="n")
       main = "", yaxt="n", xaxt="n")
  rect(xleft = 2,
       xright = 2.5,
       ybottom = case_order-1,
       ytop = case_order,
       col = colors_ordered,
       border = F)
  rect(xleft = 2.5,
       xright = 3,
       ybottom = case_order-1,
       ytop = case_order,
       col = colors_ordered3,
       border = F)
  }
dev.off()


if(F){
  save.image(file = "data/sequenza.hclust.heatmap.RData")
  write_rds(list(hclust.AB.wardD2,
                 hclust.AB.ward,
                 as.hclust(mydend)),"data/sequenza.hclust.list.Rds")
  cutree(as.hclust(mydend),k=3)
}







