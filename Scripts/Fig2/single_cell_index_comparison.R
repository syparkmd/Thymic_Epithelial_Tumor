

library(tidyverse)
#library(ComplexHeatmap)
library(ggsci)
#library(scales)
#library(circlize)
library(cowplot)
library(ggsignif)
library(patchwork)
library(ggpubr)

out_dir = "pdf_from_R"

#data path
meta_path="thymoma_meta_table.250826.tsv"
exp_tpm_path="thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv"
mixcr_path='Clone3_Fraction0.01_summary_tbl.tsv'

#reference path
bm_path='biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt'
cgs_path="Cancer_gene_census_GRCh37_v89.tsv"


# color setting
if(T){
  my_pal = pal_npg("nrc")(10)
  #scales::show_col(my_pal)
  histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
  #scales::show_col(histo_pal)
  names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
  
  stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'white')
  #show_col(stage_pal)
  names(stage_pal) <- c('I','II','III','IVa','IVb', '0')
  
  cohort_pal = pal_jama("default", alpha = 0.8)(7)[c(1,7)]
  names(cohort_pal) = c("SNUH","TCGA_CancerCell")
  #show_col(cohort_pal)
  
  gtf2i_pal = pal_npg("nrc")(10)[c(4,1,3)]
  #show_col(gtf2i_pal)
  names(gtf2i_pal) = c('m','w','c')
  
  #exp_pal = circlize::colorRamp2(c(-3,0,3), c('#253494',"gray90",'#f03b20'))
  
  myred = pal_npg("nrc")(10)[8]
  myblue = pal_npg("nrc")(10)[4]
  
}

#load meta data table
meta_dt <- read_tsv(meta_path)


#filter subgroup only
if(T){
  tcga_dt <- meta_dt %>% filter(cohort == "TCGA_CancerCell")
  snuh_dt <- meta_dt %>% filter(cohort == "SNUH")
}

#single-cell index
if(F){
  #total
  g1 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=Progenitor))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(2.6, 2.8,3.0)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Progenitor index")+
    theme(legend.position = "none")
  
  g2 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=cTEC))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(6, 7, 8)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("cTEC index")+
    theme(legend.position = "none")
  
  g3 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=mTEC))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(9, 10, 11)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("mTEC index")+
    theme(legend.position = "none")
  
  g4 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=Tuft))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(4.4, 4.8, 5.2)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Tuft index")+
    theme(legend.position = "none")
  
  gg <- wrap_plots(g1, g2, g3, g4, nrow=1)
  pdf(paste0(out_dir , "/sc_index_group.pdf"),
      width=8, height=4)
  print(gg)
  dev.off()
  
  
  #Progenitor - TCGA
  g11 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=Progenitor))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(1.6, 1.8, 2.0)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    scale_y_continuous(breaks = seq(-10, 10, 1))+
    labs(x="", y="")+
    theme_cowplot()+ggtitle("Progenitor") + 
    theme(legend.position = "none")
  
  g12<- ggplot(snuh_dt, aes(x=GTF2I_status2, y=Progenitor))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(2.2, 2.4, 2.6)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    scale_y_continuous(breaks = seq(-10, 10, 1))+
    labs(x="", y="")+
    theme_cowplot()+ggtitle("Progenitor") + 
    theme(legend.position = "none")
  
  
  g21 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=cTEC))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(6, 7, 8)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("cTEC")+
    theme(legend.position = "none")
  
  g22 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=cTEC))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(6, 7, 8)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("cTEC")+
    theme(legend.position = "none")
  
  g31 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=mTEC))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(9, 10, 11)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("mTEC")+
    theme(legend.position = "none")
  
  g32 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=mTEC))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(9, 10, 11)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("mTEC")+
    theme(legend.position = "none")
  
  g41 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=Tuft))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(4.4, 4.8, 5.2)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Tuft")+
    theme(legend.position = "none")
  
  g42 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=Tuft))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(4.4, 4.8, 5.2)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Tuft")+
    theme(legend.position = "none")
  
  gg2 <- wrap_plots(g11, g12, 
                    g21, g22,
                    g31, g32,
                    g41, g42,
                    nrow=2, byrow=F)
  pdf(paste0(out_dir , "/sc_index_subgroup.pdf"),
      width=8, height=4)
  print(gg2)
  dev.off()
  
}

#compare age among TET groups
if(F){
  nrow(meta_dt)
  a <- meta_dt %>% filter(GTF2I_status2 == "m") %>% pull(age_at_diagnosis)
  summary(a)
  b <- meta_dt %>% filter(GTF2I_status2 == "w") %>% pull(age_at_diagnosis)
  summary(b)
  wilcox.test(a,b)
  
  ggplot(meta_dt, aes(x=GTF2I_status2, y=age_at_diagnosis))+
    geom_boxplot(aes(fill = GTF2I_status2), outlier.size=-1)+
    geom_signif(comparisons = list(c(1,2), c(2,3), c(1,3)),
                y_position= c(85,90,95))+
    geom_jitter( width=0.2, height=0, alpha=0.7)+
    scale_fill_manual(values = gtf2i_pal)+
    labs(x="", y="Age at Diagnosis")+
    scale_x_discrete(labels = c("m" = "GTF2I-\ntype",
                                "w" = "CN-\ntype",
                                "c" = "Thymic\ncarcinoma"))+
    theme_cowplot()
  
}
