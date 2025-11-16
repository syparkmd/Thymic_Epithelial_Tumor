

library(tidyverse)
#library(ComplexHeatmap)
library(ggsci)
#library(scales)
#library(circlize)
library(cowplot)
library(ggsignif)
library(patchwork)
library(ggpubr)

#data path
meta_path="thymoma_meta_table.250826.tsv"
exp_tpm_path="IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv"
cibersort_mcpcount_path="cibersort_mcpcount_estimation.tsv"
mixcr_path='Clone3_Fraction0.01_summary_tbl.tsv'
pl_path="Table_S6_plasmacell_count_tsv_version.txt"


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

##immune score group comparison1
if(F){
  ggplot(meta_dt, aes(x=GTF2I_status2, y=Thymopoiesis))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.78, 0.8,0.82)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Thymopoiesis")+
    theme(legend.position = "none")
  
  g11 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=Thymopoiesis))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.78, 0.8,0.82)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Thymopoiesis")+
    theme(legend.position = "none")
  
  g12 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=Thymopoiesis))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.78, 0.8,0.82)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Thymopoiesis")+
    theme(legend.position = "none")
  
  g21 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=`Cytotoxic T cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.8, 0.82,0.84)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Cytotoxic T cell")+
    theme(legend.position = "none")
  
  g22 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=`Cytotoxic T cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.8, 0.82,0.84)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Cytotoxic T cell")+
    theme(legend.position = "none")
  
  g31 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=`Exhausted T cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.82, 0.85,0.88)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Exhausted T cell")+
    theme(legend.position = "none")
  
  g32 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=`Exhausted T cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.65, 0.68,0.71)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Exhausted T cell")+
    theme(legend.position = "none")
  
  g41 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=`B cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.82, 0.87,0.92)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("B cell")+
    theme(legend.position = "none")
  
  g42 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=`B cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.72, 0.77,0.82)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("B cell")+
    theme(legend.position = "none")
  
  gg <- wrap_plots(g11, g12, 
                   g21, g22,
                   g31, g32,
                   g41, g42,
                   nrow=2, byrow=F)

  print(gg)

}



#immune score group comparison2 (DC update)
if(F){
  #total
  g1 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=Neutrophil))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.6, 0.64,0.68)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Neutrophil")+
    theme(legend.position = "none")
  
  g2 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=`NK cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.42, 0.46,0.50)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("NK cell")+
    theme(legend.position = "none")
  
  g3 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=`Dendritic cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.68, 0.72,0.76)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Dendritic cell")+
    theme(legend.position = "none")
  
  g4 <- ggplot(meta_dt, aes(x=GTF2I_status2, y=Macrophage))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.84, 0.86,0.88)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Macrophage")+
    theme(legend.position = "none")
  
  gg <- wrap_plots(g1, g2, g3, g4, nrow=1)

  print(gg)

  
  #individual cohort
  g11 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=Neutrophil))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.6, 0.64,0.68)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Neutrophil")+
    theme(legend.position = "none")
  
  g12 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=Neutrophil))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.46, 0.5,0.54)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Neutrophil")+
    theme(legend.position = "none")
  
  g21 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=`NK cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.44, 0.48,0.52)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("NK cell")+
    theme(legend.position = "none")
  
  g22 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=`NK cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.3, 0.34,0.38)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("NK cell")+
    theme(legend.position = "none")
  
  g31 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=`Dendritic cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.68, 0.72,0.76)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Dendritic cell")+
    theme(legend.position = "none")
  
  g32 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=`Dendritic cell`))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.62, 0.66,0.70)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Dendritic cell")+
    theme(legend.position = "none")
  
  g41 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=Macrophage))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.84, 0.86,0.88)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Macrophage")+
    theme(legend.position = "none")
  
  g42 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=Macrophage))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(0.78, 0.80,0.82)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("Macrophage")+
    theme(legend.position = "none")
  
  gg2 <- wrap_plots(g11, g12, 
                    g21, g22,
                    g31, g32,
                    g41, g42,
                    nrow=2, byrow=F)

  print(gg2)

}

#immune score validation 
if(F){
  im_dt2 <- read_tsv(cibersort_mcpcount_path)
  im_dt2$estimation
  #[1] "B cell naive_CIBERSORT"                     "B cell memory_CIBERSORT"                   
  #[3] "B cell plasma_CIBERSORT"                    "T cell CD8+_CIBERSORT"                     
  #[5] "T cell CD4+ naive_CIBERSORT"                "T cell CD4+ memory resting_CIBERSORT"      
  #[7] "T cell CD4+ memory activated_CIBERSORT"     "T cell follicular helper_CIBERSORT"        
  #[9] "T cell regulatory (Tregs)_CIBERSORT"        "T cell gamma delta_CIBERSORT"              
  #[11] "NK cell resting_CIBERSORT"                  "NK cell activated_CIBERSORT"               
  #[13] "Monocyte_CIBERSORT"                         "Macrophage M0_CIBERSORT"                   
  #[15] "Macrophage M1_CIBERSORT"                    "Macrophage M2_CIBERSORT"                   
  #[17] "Myeloid dendritic cell resting_CIBERSORT"   "Myeloid dendritic cell activated_CIBERSORT"
  #[19] "Mast cell activated_CIBERSORT"              "Mast cell resting_CIBERSORT"               
  #[21] "Eosinophil_CIBERSORT"                       "Neutrophil_CIBERSORT"                      
  #[23] "T cell_MCPCOUNTER"                          "T cell CD8+_MCPCOUNTER"                    
  #[25] "cytotoxicity score_MCPCOUNTER"              "NK cell_MCPCOUNTER"                        
  #[27] "B cell_MCPCOUNTER"                          "Monocyte_MCPCOUNTER"                       
  #[29] "Macrophage/Monocyte_MCPCOUNTER"             "Myeloid dendritic cell_MCPCOUNTER"         
  #[31] "Neutrophil_MCPCOUNTER"                      "Endothelial cell_MCPCOUNTER"               
  #[33] "Cancer associated fibroblast_MCPCOUNTER"   
  
  tmp_dt <- left_join(meta_dt, im_dt2 %>% as.data.frame() %>% column_to_rownames("estimation") %>%
                        t() %>% as.data.frame() %>% rownames_to_column("id")  %>% as_tibble())
  tmp_dt <- tmp_dt %>% mutate(CIBERSORT_B = 
                                `B cell naive_CIBERSORT`+ `B cell memory_CIBERSORT`+
                                `B cell plasma_CIBERSORT`,
                              CIBERSORT_T = 
                                `T cell CD4+ naive_CIBERSORT` + `T cell CD4+ memory resting_CIBERSORT` + 
                                `T cell CD4+ memory activated_CIBERSORT` + `T cell follicular helper_CIBERSORT` +
                                `T cell regulatory (Tregs)_CIBERSORT`+`T cell gamma delta_CIBERSORT`,
                              CIBERSORT_NK = 
                                `NK cell resting_CIBERSORT`+`NK cell activated_CIBERSORT`,
                              CIBERSORT_DC = 
                                `Myeloid dendritic cell resting_CIBERSORT`+`Myeloid dendritic cell activated_CIBERSORT`,
                              CIBERSORT_Mac_Mono=
                                `Monocyte_CIBERSORT`+`Macrophage M0_CIBERSORT`+
                                `Macrophage M1_CIBERSORT`+`Macrophage M2_CIBERSORT`
  )
  
  #correlation
  #Thymopoiesis vs T
  g11 <- ggplot(tmp_dt, aes(x=Thymopoiesis, y=CIBERSORT_T))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  
  g12 <- ggplot(tmp_dt, aes(x=Thymopoiesis, y=`T cell_MCPCOUNTER`))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  #B cell
  g21 <- ggplot(tmp_dt, aes(x=`B cell`, y=CIBERSORT_B))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  g22 <- ggplot(tmp_dt, aes(x=`B cell`, y=`B cell_MCPCOUNTER`))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  #NK cell
  g31 <- ggplot(tmp_dt, aes(x=`NK cell`, y=CIBERSORT_NK))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  g32 <- ggplot(tmp_dt, aes(x=`NK cell`, y=`NK cell_MCPCOUNTER`))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  #Macrophage
  g41 <- ggplot(tmp_dt, aes(x=`Macrophage`, y=CIBERSORT_Mac_Mono))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  g42 <- ggplot(tmp_dt, aes(x=`Macrophage`, y=`Macrophage/Monocyte_MCPCOUNTER`))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  #Dendritic cell
  g51 <- ggplot(tmp_dt, aes(x=`Dendritic cell`, y=CIBERSORT_DC))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  g52 <- ggplot(tmp_dt, aes(x=`Dendritic cell`, y=`Myeloid dendritic cell_MCPCOUNTER`))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  
  #Neurtophil
  g61 <- ggplot(tmp_dt, aes(x=`Neutrophil`, y=Neutrophil_CIBERSORT))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  g62 <- ggplot(tmp_dt, aes(x=`Neutrophil`, y=`Neutrophil_MCPCOUNTER`))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_cowplot()
  
  gg <- wrap_plots(g11, g21, g31, g41, g51, g61,
                   g12, g22, g32, g42, g52, g62, 
                   nrow=2, byrow=T)
  print(gg)

  
  #Boxplot for B cell
  g1 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=`B cell`))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="B cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("Our B cell")
  
  g2 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=CIBERSORT_B))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="B cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("CIBERSORT B cell")
  
  g3 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=`B cell_MCPCOUNTER`))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="B cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("MCPcounter B cell")
  
  wrap_plots(g1, g2, g3, nrow=1)  
  
  
  #Boxplot for T cell
  g1 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=`Thymopoiesis`))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="Immature T cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("Our immature T cell")
  
  g2 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=CIBERSORT_T))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="T cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("CIBERSORT T cell")
  
  g3 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=`T cell_MCPCOUNTER`))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="T cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("MCPcounter T cell")
  
  wrap_plots(g1, g2, g3, nrow=1)  
  
  #Boxplot for NK cell
  g1 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=`NK cell`))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="NK cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("Our NK cell")
  
  g2 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=CIBERSORT_NK))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="NK cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("CIBERSORT NK cell")
  
  g3 <- ggplot(tmp_dt, aes(x=GTF2I_status2, y=`NK cell_MCPCOUNTER`))+
    geom_boxplot(aes(fill = GTF2I_status2))+
    stat_compare_means(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      method = "wilcox")+
    scale_fill_manual(values = gtf2i_pal)+
    scale_x_discrete(label = c("m" = "GTF2I-type",
                               "w" = "CN-type",
                               "c" = "Thymic carcinoma"))+
    labs(x="", y="NK cell score")+
    theme_cowplot()+theme(legend.position = "none")+
    ggtitle("MCPcounter NK cell")
  
  wrap_plots(g1, g2, g3, nrow=1)  
  
  
}


#Plasma cell count and B cell score relationship
if(F){
  pl_dt <- read_tsv(pl_path)
  colnames(pl_dt) <- gsub(" ", "_", colnames(pl_dt))
  pl_dt <- pl_dt %>% dplyr::rename(B_cell_read = Read_sum_of_expanded_B_cell_clones,
                                   B_cell_score = B_cell_signature_enrichment_score,
                                   plasma_infiltration = `Histology_reivew:_\nPlasma_cell_infiltration`,
                                   TLS = `Histology_review:\nTertiary_lymphoid_structure`)
  f_pl_dt <- pl_dt %>% filter(is.na(plasma_infiltration)==F)
  f_pl_dt <- f_pl_dt %>% mutate(plasma_cell_group = 
                                  case_when(plasma_infiltration %in% c(0,1) ~ "Rare",
                                            plasma_infiltration == 2 ~ "Moderate",
                                            plasma_infiltration == 3 ~ "Frequent"))
  f_pl_dt$plasma_cell_group <- factor(f_pl_dt$plasma_cell_group, 
                                      levels = c("Rare", "Moderate", "Frequent"))
  g <- ggplot(f_pl_dt, aes(x=plasma_cell_group, y=B_cell_score))+
    geom_boxplot(aes(fill = plasma_cell_group))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3),
                         c(1,3)),
      y_position = c(0.8, 0.85, 0.9))+
    geom_jitter(height=0, width=0.2)+
    scale_fill_manual(values = 
                        c("Rare" = "#91CDDB",
                          "Moderate" ="#3A8EA2",
                          "Frequent" = "#1F4B56"))+
    labs(x= "Plasma cell frequency in H&E",
         y=" B cell score in RNAseq")+
    theme_cowplot()+ theme(legend.position = "None")
  
  print(g)

  
  #add cytotoxic T cell score  
  f_pl_dt <- left_join(f_pl_dt,
                       meta_dt %>% dplyr::select(id, `Cytotoxic T cell`) %>%
                         dplyr::rename(ID = id))
  g <- ggplot(f_pl_dt, aes(x=B_cell_score, y=`Cytotoxic T cell`))+
    geom_point(aes(color=TLS))+
    geom_text_repel(data=subset(f_pl_dt, TLS==TRUE),
                    aes(label= ID))+
    scale_color_manual(values = c("TRUE" = "#3D2F14",
                                  "FALSE" = "#FDBB2B"))+
    theme_cowplot()
  

  print(g)

  
}


#immune score and Myasthenia gravis
if(F){
  meta_dt$history_myasthenia_gravis %>% table()
  # NO YES 
  #101  35 
  cell_type_list <- c("Thymopoiesis","Cytotoxic T cell","Exhausted T cell","B cell",
                      "Neutrophil" ,"NK cell", "Dendritic cell", "Macrophage"
  ) 
  
  plot_list=list()
  for(i in 1:length(cell_type_list)){
    t_cell_type = cell_type_list[i]
    print(t_cell_type)
    t_dt <- meta_dt %>% dplyr::select(id, history_myasthenia_gravis,
                               all_of(t_cell_type)) %>% 
      dplyr::rename(score = t_cell_type) %>%
      filter(is.na(history_myasthenia_gravis)==F)
    plot_list[[i]] <- ggplot(t_dt, aes(x=history_myasthenia_gravis, y=score))+
      geom_boxplot()+
      geom_signif(comparisons = list(c(1,2)))+
      theme_cowplot()+
      ggtitle(t_cell_type)
  }
  
  wrap_plots(plot_list)
  
  #only B cell score showed significant difference
  a <- meta_dt %>% filter(history_myasthenia_gravis == "YES") %>% pull(`B cell`)
  b <- meta_dt %>% filter(history_myasthenia_gravis == "NO") %>% pull(`B cell`)
  mean(a) #0.447
  mean(b) #0.392
  wilcox.test(a,b) #p = 0.02
  
  tmp_dt <- meta_dt %>% filter(is.na(history_myasthenia_gravis)==F)
  tmp_dt$history_myasthenia_gravis <- factor(tmp_dt$history_myasthenia_gravis,
                                             levels = c("YES","NO"))
  g <- ggplot(tmp_dt, aes(x=history_myasthenia_gravis, y=`B cell`))+
    geom_boxplot(aes(fill= history_myasthenia_gravis))+
    geom_signif(comparisons = list(c(1,2)))+
    scale_fill_manual(values = c("YES" = "#7d5e25",
                                 "NO" = "#fde2a9"))+
    labs(x = "History of myasthenia gravis",
         y= "B cell score")+
    theme_cowplot()+theme(legend.position = "none")
  

  print(g)

  
  #only in cn type
  tmp_dt <- meta_dt %>% filter(GTF2I_status2 == "w")
  a <- tmp_dt %>% filter(history_myasthenia_gravis == "YES") %>% pull(`B cell`)
  b <- tmp_dt %>% filter(history_myasthenia_gravis == "NO") %>% pull(`B cell`)
  mean(a) #0.48
  mean(b) #0.46
  wilcox.test(a,b) #p = 0.3815 not sig
  
  #muliple linear regression
  lm_res <- lm(`B cell` ~ GTF2I_status2 + history_myasthenia_gravis, 
               data = meta_dt)
  
  s_lm_res <- summary(lm_res)
  print(s_lm_res)
  
  #make lm table
  lm_dt <- s_lm_res$coefficients %>% as.data.frame() %>% 
    rownames_to_column("item") %>% as_tibble()
  #add CI
  ci_dt <- confint(lm_res) %>% as.data.frame() %>%
    rownames_to_column("item") %>% as_tibble()
  lm_dt <- left_join(lm_dt, ci_dt) 
  lm_dt <- lm_dt %>% dplyr::rename(pval = `Pr(>|t|)`,
                                   low_ci = `2.5 %`,
                                   high_ci = `97.5 %`)
  
  #item                         Estimate `Std. Error` `t value`     pval  low_ci high_ci
  #<chr>                           <dbl>        <dbl>     <dbl>    <dbl>   <dbl>   <dbl>
  #1 (Intercept)                    0.340        0.0161    21.1   4.30e-44  0.308   0.372 
  #2 GTF2I_status2w                 0.121        0.0260     4.65  7.92e- 6  0.0696  0.173 
  #3 GTF2I_status2c                 0.157        0.0398     3.94  1.33e- 4  0.0780  0.235 
  #4 history_myasthenia_gravisYES   0.0177       0.0286     0.621 5.36e- 1 -0.0388  0.0743
  
  g <- ggplot(lm_dt, aes(x=item, y=Estimate))+
    geom_point()+
    geom_errorbar(aes(ymin = low_ci, ymax = high_ci))+
    geom_hline(yintercept=0, linetype = "dashed")+
    coord_flip()+
    theme_cowplot()
  

  print(g)

  
}


