#Fig. 1A oncogrid
#drawing code from kjyi

library(ComplexHeatmap)
library(tidyverse)
library(circlize)

meta_data_path='thymoma_meta_table.250826.tsv'
oncogrid_out_path="fig1a_oncogrid.pdf"

meta_dt <- read_tsv(meta_data_path)
ready_dt <- meta_dt %>% dplyr::select(id,
                                      histology = histologic_type2,
                                      cohort, 
                                      purity=final_cellularity, 
                                      TMB=TMB,
                                      GTF2I,
                                      GTF2I_TCGA=TCGA_paper_GTF2Imt,
                                      HRAS,
                                      CYLD,
                                      final_ploidy,
                                      #CN1q=`1q`,
                                      #CN7p=`7p`,
                                      #CN7q=`7q`,
                                      #CN6p=`6p`,
                                      #CN6q=`6q`,
                                      #CN14q=`14q`,
                                      #CN16q=`16q`,
                                      #CNXq=Xq,
                                      amp_prop,
                                      del_prop,
                                      unchanged_prop) %>%
  mutate(rescue = (!is.na(GTF2I_TCGA))&GTF2I_TCGA == "w" & (!is.na(GTF2I))) %>%
  arrange(GTF2I, histology %in% c("CA-SqCC","CA-UN","NE"), histology,HRAS,CYLD,rescue)

ready_dt$GTF2I[is.na(ready_dt$GTF2I)] = "NA"
ready_dt$HRAS[is.na(ready_dt$HRAS)] = "NA"
ready_dt$CYLD[is.na(ready_dt$CYLD)] = "NA"


# background and setting
table(ready_dt$GTF2I)
width1=0.45 # width of a cell is 0.82
height1=0.45
gap1=0.3 # row gap btw big class
gap2 = 1.5 # column gap
x=1:nrow(ready_dt)
x[72:137] = x[72:137] + gap2
x[125:137] = x[125:137] + gap2


cairo_pdf(oncogrid_out_path,
          width=700/72,
          height = 520/72,
          pointsize = 10)
# plot(cars)
# dev.off()

# par(mar=c(1,10,1,10))
par(mar=c(1,10,1,1))

plot(1,type='n',xlim=range(x)+c(-0.5,0.5), ylim=c(25,1)+c(+0.5,-0.5),
     yaxt="n",ylab="",
     xaxt="n",xlab="",
     bty="n",
     yaxs='i', xaxs='i')
# axis(4,at = 1:20,las=2)
current=0


# header grouping
current = current + 1
group_pal=c("#364D7C","#CE4935","#25917C")
rect(xleft  = x[1]-width1,
     xright = x[71]+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = group_pal[1],
     border=NA)
rect(xleft  = x[72]-width1,
     xright = x[124]+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = group_pal[2],
     border=NA)
rect(xleft  = x[125]-width1,
     xright = x[137]+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = group_pal[3],
     border=NA)
text(x=x[1]/2+x[71]/2,y=current,"GTF2I-type",col="white",font=2)
text(x=x[72]/2+x[124]/2,y=current,"CN-type",col="white",font=2)
text(x=x[125]/2+x[137]/2,y=current,"TC",col="white",font=2)


# gap
current = current + gap1

# burden
current = current + 4
ready_dt$TMB

rect(xleft  = x-width1,
     xright = x+width1,
     ybottom = current + height1,
     ytop=current + height1 - pmin(ready_dt$TMB,4),
     col = ifelse(ready_dt$TMB > 3,"red","grey50"),
     border=NA)
at = seq(current+height1,current+height1-4,length.out = 5)
val=seq(0,4,length.out = 5)
axis(side = 2,at=at,labels=val,las=2)

axis(side = 2,at=at[3],labels="Mutational burden    \n (/Mbps)    ",las=2)

# gap
current = current + gap1

# histologic type
current = current + 1
histo_pal <- c("MN-T"="#564E8B", 
               "A"="#57246C",
               "AB"="#334382", 
               "B1"="#911D51",
               "B2"="#A72124",
               "B3"="#D01A1B", 
               "CA-SqCC"="#137977", 
               "CA-UN"="black",
               "NE"="#8AB435")
rect(xleft  = x-width1,
     xright = x+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = histo_pal[ready_dt$histology],
     border=NA)
axis(side = 2,at = current,labels = "Hitologic type",tick = F,las=2)

# gap
current = current + gap1

# GTF2I
current = current + 1
snp_pal = c("Missense mutation"="#09417B",
            "Nonsense mutation"="#D01B1B",
            # "Splice site mutation" = "#42A341",
            # "In-frame indel"="#208AA2",
            "Frame-shift indel"="#84578F",
            "NA"="#00000010")
snp_translate = c("nonframeshift insertion" = "Frame-shift indel",
                  "nonsynonymous SNV" = "Missense mutation",
                  "frameshift deletion" = "Frame-shift indel",
                  "stopgain" = "Nonsense mutation",
                  "NA"="NA")
rect(xleft  = x-width1,
     xright = x+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = snp_pal[snp_translate[ready_dt$GTF2I]],
     border=NA)
rect(xleft  = x-0.35,
     xright = x+0.35,
     ybottom = current - 0.15,
     ytop=current + 0.15,
     col=ifelse(ready_dt$rescue,"yellow",NA),
     border=NA)

axis(side = 2,at = current,labels = expression(italic("GTF2I")),tick = F,las=2)
# HRAS
current = current + 1
rect(xleft  = x-width1,
     xright = x+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = snp_pal[snp_translate[ready_dt$HRAS]],
     border=NA)
axis(side = 2,at = current,labels = expression(italic("HRAS")),tick = F,las=2)
# CYLD
current = current + 1
rect(xleft  = x-width1,
     xright = x+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = snp_pal[snp_translate[ready_dt$CYLD]],
     border=NA)
axis(side = 2,at = current,labels = expression(italic("CYLD")),tick = F,las=2)
# gap
current = current + gap1*2


# cohort info # 2025-07-27
cohort_pal = c("SNUH"="red","TCGA_CancerCell"="blue")
current = current + 1
rect(xleft  = x-width1,
     xright = x+width1,
     ybottom = current - height1,
     ytop=current + height1,
     col = cohort_pal[ready_dt$cohort],
     border=NA)
axis(side = 2,at = current,labels = "Cohort",tick = F,las=2)

cn_pal = colorRamp2(c(0,2,5), c('#253494',"gray90",'#f03b20'))
lgd1 = Legend(title="Histologic type",
              labels = c("Micronodular thymoma",
                         "A","AB","B1","B2","B3",
                         "Squamous cell carcinoma",
                         "Undifferentiated carcinoma",
                         "Neuroendocrine carcinoma"),
              legend_gp=gpar(fill=histo_pal),ncol=3)
lgd2 = Legend(title="Mutation type",
              labels = names(snp_pal)[1:3],
              legend_gp=gpar(fill=snp_pal[1:3]))
lgd3 = Legend(col_fun = cn_pal, title = "Copy number",direction = "horizontal",
              legend_width = unit(3.5, "cm"))
lgd4 = Legend(title = "Cohort",direction = "horizontal",
              legend_gp = gpar(fill=cohort_pal), labels=c("TCGA_CancerCell"="TCGA","SNUH"="This Study")[names(cohort_pal)])
pd = packLegend(list = list(lgd1, lgd2, lgd3, lgd4),direction = "horizontal")
draw(pd, x = unit(0.5, "npc"), y = unit(0.5, "cm"), just = c("bottom"))

dev.off()

