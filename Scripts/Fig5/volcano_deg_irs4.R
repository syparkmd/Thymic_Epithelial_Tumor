

library(tidyverse)
#library(ComplexHeatmap)
library(ggsci)
#library(scales)
#library(circlize)
library(cowplot)
library(ggsignif)
library(patchwork)
library(ggpubr)
library(ggrepel)

out_dir = "pdf_from_R"

#data path
meta_path="thymoma_meta_table.250826.tsv"
exp_tpm_path="thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv"
mixcr_path='Clone3_Fraction0.01_summary_tbl.tsv'

#reference path
bm_path='biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename_190522.txt'
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

#function
my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


#load meta data table
meta_dt <- read_tsv(meta_path)

# load biomart data
bm_dt <- read_tsv(bm_path)
bm_dt <- bm_dt %>% rename(gene = `Gene name`, gene_type = `Gene type`) %>% select(gene, gene_type) %>% unique()

# load expression data
exp_dt <- read_tsv(exp_tpm_path)
exp_dt_pcg <- left_join(exp_dt, bm_dt) %>% filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), funs(log10(.+0.01)))
nrow(l10_exp_dt)
l10_exp_dt <- left_join(l10_exp_dt, bm_dt) %>% filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
nrow(l10_exp_dt)

# load cancer gene sensus
cgs_dt <- read_tsv(cgs_path)
onc_dt <- cgs_dt %>% rename(gene = `Gene Symbol`, role_in_cancer = `Role in Cancer`) %>% select(gene, role_in_cancer) %>% filter(grepl('oncogene',role_in_cancer)== T)
onc_dt$role_in_cancer <- 'oncogene'


#filter subgroup only
if(T){
  tcga_dt <- meta_dt %>% filter(cohort == "TCGA_CancerCell")
  snuh_dt <- meta_dt %>% filter(cohort == "SNUH")
}

#draw volcano plot
#define WT samples and Mut samples
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']
CA_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'c']

# volcano plot - total sample
volc_dt <- l10_exp_dt
volc_dt$MT_mean <- rowMeans(volc_dt[,MT_ids])
volc_dt$WT_mean <- rowMeans(volc_dt[,WT_ids])
volc_dt$CA_mean <- rowMeans(volc_dt[,CA_ids])
volc_dt$WTCA_mean <- rowMeans(volc_dt[,c(CA_ids,WT_ids)])
volc_dt$WTMT_mean <- rowMeans(volc_dt[,c(MT_ids,WT_ids)])



volc_dt$t.test.pvalue <- apply(volc_dt[,c(MT_ids, WT_ids)],1,function(x) my.t.test.p.value(x[MT_ids], x[WT_ids]))

volc_dt <- left_join(volc_dt, onc_dt)

ggplot(volc_dt, aes(x=WT_mean - MT_mean, y=-log10(t.test.pvalue)))+
  geom_point(aes(color = role_in_cancer), size=3, alpha=0.7)+
  geom_text_repel(data = subset(volc_dt, abs(WT_mean - MT_mean) >=2 | -log10(t.test.pvalue) >=40), aes(label = gene))+
  geom_hline(yintercept = -log10(0.05/nrow(volc_dt)), linetype="longdash")+
  scale_x_continuous(limits = c(-3,3))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  ggtitle('Total samples')


#SNUH cohort only
MT_ids <- snuh_dt$id[snuh_dt$GTF2I_status2 == 'm']
WT_ids <- snuh_dt$id[snuh_dt$GTF2I_status2 == 'w']

volc_dt <- l10_exp_dt[,c('gene',MT_ids, WT_ids)]
volc_dt$MT_mean <- rowMeans(volc_dt[,MT_ids])
volc_dt$WT_mean <- rowMeans(volc_dt[,WT_ids])

volc_dt$t.test.pvalue <- apply(volc_dt[,c(MT_ids, WT_ids)],1,function(x) my.t.test.p.value(x[MT_ids], x[WT_ids]))
volc_dt <- left_join(volc_dt, onc_dt)

ggplot(volc_dt, aes(x=WT_mean - MT_mean, y=-log10(t.test.pvalue)))+
  geom_point(aes(color = role_in_cancer), size=3, alpha=0.7)+
  geom_text_repel(data = subset(volc_dt, abs(WT_mean - MT_mean) >=2 | -log10(t.test.pvalue) >=40), aes(label = gene))+
  geom_hline(yintercept = -log10(0.05/nrow(volc_dt)), linetype="longdash")+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  ggtitle('SNUH samples')


#TCGA cohort only
MT_ids <- tcga_dt$id[tcga_dt$GTF2I_status2 == 'm']
WT_ids <- tcga_dt$id[tcga_dt$GTF2I_status2 == 'w']

volc_dt <- l10_exp_dt[,c('gene',MT_ids, WT_ids)]
volc_dt$MT_mean <- rowMeans(volc_dt[,MT_ids])
volc_dt$WT_mean <- rowMeans(volc_dt[,WT_ids])

volc_dt$t.test.pvalue <- apply(volc_dt[,c(MT_ids, WT_ids)],1,function(x) my.t.test.p.value(x[MT_ids], x[WT_ids]))
volc_dt <- left_join(volc_dt, onc_dt)

ggplot(volc_dt, aes(x=WT_mean - MT_mean, y=-log10(t.test.pvalue)))+
  geom_point(aes(color = role_in_cancer), size=3, alpha=0.7)+
  geom_text_repel(data = subset(volc_dt, abs(WT_mean - MT_mean) >=2 | -log10(t.test.pvalue) >=40), aes(label = gene))+
  geom_hline(yintercept = -log10(0.05/nrow(volc_dt)), linetype="longdash")+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  ggtitle('TCGA samples')


#IRS4 expression comparison
if(F){
  #total
  ggplot(meta_dt, aes(x=GTF2I_status2, y=corrected_IRS4_TPM))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(480, 500,520)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="")+
    ggtitle("IRS4")+
    theme(legend.position = "none")
  
  #TCGA 
  g1 <- ggplot(tcga_dt, aes(x=GTF2I_status2, y=corrected_IRS4_TPM))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(100, 110,120)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="IRS4 expression (TPM)")+
    ggtitle("TCGA")+
    theme(legend.position = "none")
  
  g2 <- ggplot(snuh_dt, aes(x=GTF2I_status2, y=corrected_IRS4_TPM))+
    geom_boxplot(aes(fill=GTF2I_status2))+
    geom_signif(
      comparisons = list(c(1,2), c(2,3), c(1,3)),
      y_position = c(460, 480,500)
    )+
    scale_fill_manual(values= gtf2i_pal)+
    theme_cowplot()+
    labs(x="", y="IRS4 expression (TPM)")+
    ggtitle("SNUH")+
    theme(legend.position = "none")
  gg <- wrap_plots(g1, g2, nrow=1)
  
  pdf(paste0(out_dir , "/irs4_subgroup.pdf"),
      width=8, height=4)
  print(gg)
  dev.off()
}

