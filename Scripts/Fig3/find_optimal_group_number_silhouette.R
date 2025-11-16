

library(factoextra)
library(cluster) # For silhouette, pam, etc.
library(tidyverse)
library(ggsci)
library(circlize)



# sypark's code --------------------------------------------------------------------------------------------------------
meta_dt <- read_tsv('thymoma_meta_table.250826.tsv')
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']
CA_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'c']
bm_dt <- read_tsv('biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename_190522.txt') %>% 
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()
# load expression data
exp_dt <- read_tsv('thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')
# exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
#   filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only

#color setting
my_pal = pal_npg("nrc")(10)
meta_dt$histologic_type %>% unique()
histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
# histo_pal = c("#5C2371", "#37468A", "#5A5193", "#981C54", "#B11D23", "#DE1115", "#157D7B", "#90BD32")

histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
# show_col(histo_pal)
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")


names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'white')
names(stage_pal) <- c('I','II','III','IVa','IVb', '0')
cohort_pal = pal_jama("default", alpha = 0.8)(7)[c(1,7)]
names(cohort_pal) = c("SNUH","TCGA_CancerCell")
gtf2i_pal = pal_npg("nrc")(10)[c(4,1,3)]
names(gtf2i_pal) = c('m','w','c')
exp_pal = colorRamp2(c(-3,0,3), c('#253494',"gray90",'#f03b20'))
exp_pal2 = colorRamp2(c(-1,0,1), c('#253494',"gray90",'#f03b20'))
# purity_pal = colorRamp2(c(0,0.3,1), c("red","#272822","#999488"))
purity_pal = colorRamp2(c(0,1), c("white", pal_jama("default", alph = 1)(7)[1]))
myred = pal_npg("nrc")(10)[8]
myblue = pal_npg("nrc")(10)[4]

# unsupervised clustering
hm_dt <- l10_exp_dt
hm_dt <- hm_dt %>% as.data.frame() %>% column_to_rownames('gene') 
hm_dt$var <- apply(hm_dt, 1, function(x) var(x))
hm_dt <- hm_dt %>% rownames_to_column('gene')
var_genes <- hm_dt %>% arrange(desc(var))%>%.$gene%>%.[1:2500] # top 2500 high variance genes, variance in log10 scale
hm_dt <- hm_dt %>% dplyr::filter(gene %in% var_genes) %>% dplyr::select(-var) %>% column_to_rownames('gene') %>%
  as.matrix()
hm_annot_dt <- meta_dt %>% dplyr::select(id, cohort, GTF2I_status2, histologic_type, Stage, ImmuneScore,Purity=final_cellularity) %>% 
  as.data.frame() %>% column_to_rownames('id')
hm_annot_dt <- hm_annot_dt[colnames(hm_dt),]
cohort_pal2 = cohort_pal
names(cohort_pal2)[2] = "TCGA"

hm_dt <- t(scale(t(hm_dt))) # row scaling
dim(hm_dt)
if(F){
  top_anno <- ComplexHeatmap::HeatmapAnnotation("Group"= hm_annot_dt$GTF2I_status2,
                                                "Histologic type"= hm_annot_dt$histologic_type,
                                                "Purity" = hm_annot_dt$Purity,
                                                # "Cohort" = hm_annot_dt$cohort %>% str_replace("_CancerCell",""),
                                                col = list("Group" = gtf2i_pal,
                                                           "Histologic type" = histo_pal,
                                                           "Purity" = purity_pal),
                                                # gp = gpar(cex=0.75),
                                                annotation_legend_param = list(
                                                  "Histologic type" = list(
                                                    title="Histologic type",
                                                    at = c("A","AB","MN-T","B1","B2","B3","NE","TC"),
                                                    labels = c("Type A thymoma",
                                                               "Type AB thymoma",
                                                               "Micronodular thymoma with lymphoid stroma",
                                                               "Type B1 thymoma",
                                                               "Type B2 thymoma",
                                                               "Type B3 thymoma",
                                                               "Thymic NEC",
                                                               # "Thymic carcinoma (Squamous cell carcinoma,Undifferentiated ca)")
                                                               "Thymic carcinoma")
                                                  ),
                                                  "Group" = list(
                                                    title = "Group",
                                                    at = c("m", "w","c"),
                                                    labels = c("GTF2I-mutant", "GTF2I wild-type","Thymic carcinoma"))))
  
  ComplexHeatmap::Heatmap(hm_dt, clustering_distance_columns = 'pearson',
                          clustering_method_columns =  "average",
                          clustering_method_rows = "ward.D",
                          top_annotation = top_anno, col= exp_pal, show_row_names = F)
}

#find the optimal number of cluster
exp_dt <- t(hm_dt)


# Silhouette method 
pdf("silhouette.pdf",
    width=6, height=6)
fviz_nbclust(exp_dt, kmeans, method = "silhouette", k.max = 10) +
  labs(subtitle = "Silhouette method")
dev.off()
# definitely k=3
