#Tumor cell fraction(TCF), Tumor mutational burden (TMB), Copy number alteration(CNA)


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

#tcf and N. point mutations check
if(F){
  tmp_dt <- meta_dt
  tmp_dt$final_cellularity[is.na(tmp_dt$final_cellularity)] <- 0.1
  ggplot(tmp_dt, aes(x=final_cellularity, y=n_pointmt/bait_size))+
    geom_point(aes(color = histologic_type))+
    scale_color_manual(values = histo_pal)+
    coord_cartesian(ylim=c(0,2.5))+
    theme_cowplot()
  
}

#TCF and Nmut cor
if(F){
  cor.test(meta_dt$TMB, meta_dt$final_cellularity)
  tmp_dt <- meta_dt %>% filter(TMB < 5)
  tmp_dt$final_cellularity[is.na(tmp_dt$final_cellularity)] <- 0.1
  cor.test(tmp_dt$TMB, tmp_dt$final_cellularity)
  #data:  tmp_dt$TMB and tmp_dt$final_cellularity
  #t = 7.5771, df = 134, p-value = 5.199e-12
  #alternative hypothesis: true correlation is not equal to 0
  #95 percent confidence interval:
  #  0.4178580 0.6555643
  #sample estimates:
  #  cor 
  #0.5476683 
  
}


#copy number altered proportion
if(F){
  #total
  ggplot(meta_dt, aes(x=GTF2I_status2, y=1-unchanged_prop))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(1, 1.1,1.2)
    )+
    geom_boxplot(aes(fill=GTF2I_status2))+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()
  a <- meta_dt %>% filter(GTF2I_status2 == "m" &
                            final_selection != "low") %>% pull(unchanged_prop)
  b <- meta_dt %>% filter(GTF2I_status2 == "w" &
                            final_selection != "low") %>% pull(unchanged_prop)
  t.test(1-a, 1-b) #5.1e-8
  wilcox.test(1-a, 1-b)  #5.8e-12
  
  
  #TCGA
  g1 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=1-unchanged_prop))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(1, 1.1,1.2)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    scale_y_continuous(breaks = seq(0, 1, 0.5))+
    labs(x="", y="Copy number altered proportion")+
    theme_cowplot()+ggtitle("TCGA")
  #ggsignif: wilcox
  
  #SNUH
  g2 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=1-unchanged_prop))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.85, 0.9, 0.95)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    scale_y_continuous(breaks = seq(0, 1, 0.5))+
    labs(x="", y="Copy number altered proportion")+
    theme_cowplot()+ggtitle("SNUH")
  
  gg <- wrap_plots(g1, g2, ncol=1)
  pdf(paste0(out_dir,'/cna_prop_subgroup.pdf'), 
      width=6, height=6)
  print(gg)
  dev.off()
}

#mutational burden (TMB = n_pointmt/bait_size)
if(F){
  #total
  ggplot(meta_dt, aes(x=GTF2I_status2, y=TMB))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(10, 12,14)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    #coord_cartesian(ylim = c(0, 3))+
    theme_cowplot()
  
  #wilcox
  # m, w, :0.031
  # w, c : 5.2e-05
  # m, c: 5.1e-05
  
  
  #TCGA
  g1 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=TMB))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    scale_fill_manual(values= gtf2i_pal)+
    scale_y_continuous(breaks = seq(0, 10, 1))+
    coord_cartesian(ylim=c(0,3))+
    labs(x="", y="Mutational burden (/Mbps)")+
    theme_cowplot()+ggtitle("TCGA")
  #ggsignif: wilcox
  #m, w, : 0.0029
  #w,c, : 0.00039
  #m, c: 0.001
  
  #SNUH
  g2 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=TMB))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(2.6, 2.7, 2.8)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    scale_y_continuous(breaks = seq(0, 3,1))+
    labs(x="", y="Mutational burden (/Mbps)")+
    theme_cowplot()+ggtitle("SNUH")
  
  gg <- wrap_plots(g1, g2, ncol=1)
  pdf(paste0(out_dir,'/tmb_subgroup_box.pdf'), 
      width=6, height=6)
  print(gg)
  dev.off()
}

