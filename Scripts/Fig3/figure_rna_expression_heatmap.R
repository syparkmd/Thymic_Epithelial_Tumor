#copied from kjyi script
# save g1,2,3 genes.


library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(ggrepel)
# require(GSEABase)
library(seriation)

# sypark's code --------------------------------------------------------------------------------------------------------
meta_dt <- read_tsv(paste0('thymoma_meta_table.250826.tsv'))
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']
CA_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'c']
bm_dt <- read_tsv('biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt') %>% 
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()
# load expression data
exp_dt <- read_tsv('thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')
# exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
#   filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
nrow(l10_exp_dt) #19900

#color setting
if(T){
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
}

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
  top_anno <- HeatmapAnnotation(df = hm_annot_dt)
  Heatmap(hm_dt, clustering_distance_columns = 'pearson',
          clustering_method_columns =  "average",
          clustering_method_rows = "ward.D",
          top_annotation = top_anno, col= exp_pal, show_row_names = F,
          show_column_names = F, row_split = 6)
}



# gene set prep --------------------------------------------------------------------------------------------------------
h_list <- GSEABase::getGmt('h.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
# c2_list <- getGmt('c2.cp.v6.2.symbols.gmt') %>% geneIds()
c5_list <- GSEABase::getGmt('c5.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
# all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
h_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
c("GDF5","GDF1","BMP2","BMP4","LEFTY2","BMP8B","LRG1","SMAD7","GDF6","BMP6","TGFB2") %in% h_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
c("ADAMTS20","COL9A3","MMP3","THSD4","COL11A1") %in% h_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

# names(all_list)[grepl("DEVELOP",names(all_list))&grepl("THYM",names(all_list))]
# 
# names(all_list)[grepl("GO",names(all_list))&grepl("_T_CELL",names(all_list))]

all_list <- list("T cell development" = c("CD1A","DNTT","CD3D","CD3E", "GPAP2","LCK", "CD8A"),
                 "Metabolic process"=c("CYP2C9","CYP3A5","CYP3A4","STAR","UGT2B7","PPARGC1A","AKR1D1","SRD5A1","WNT4","SLC34A1","ASS1"),
                 # "ECM organization" = c("ADAMTS20","COL9A3","MMP3","THSD4","COL11A1"),
                 # "TGFβ signaling"=c("GDF5","GDF1","BMP2","BMP4","LEFTY2","BMP8B","LRG1","SMAD7","GDF6","BMP6","TGFB2"),
                 "Epithelial-mesenchymal transition"=c("COL11A1","COL1A1", "MMP3","CDH11", "CDH6", "WNT5A"),
                 "Development"=c("TBX1","MYH6","WNT2","HMGA2","WNT5A","BMP4","CX3CR1","IRX2","IRX4","SALL1"),
                 "TNFα signaling/inflammation"=c("CCL7","CXCL10","IL6","CCL20","CXCL13"))
selected_names <- names(all_list)

# library(piano)
# gslist2gsc <- function(gslist){
#   lapply(names(gslist),function(gsname){cbind(gslist[[gsname]],gsname)}) %>% do.call(rbind,.) %>% loadGSC}
# gsares <- runGSAhyper(genes=all_names,
#                       gsc = gslist2gsc(all_list))

if(F){
  gsares$pvalues %>% sort %>% names %>% head(100) %>% paste0('"',.,'"') %>% cat(sep=",\n")
  
  gsares$pvalues %>% sort %>% names %>% head(200) %>% .[grepl("EXTRACELLULAR",.)] -> x
  gsares$pvalues %>% sort %>% names %>% head(200) %>% .[grepl("EMBRY",.)] -> x
  structure(lapply(x,function(x){table(all_list[[x]]%in%rownames(hm_dt))}),names=x)
}



length(selected_names)
myninecolor=RColorBrewer::brewer.pal(9,"Set1")[c(1,5,2,4,3,6)]
myninecolor=c("#E41A1C",
              "#FF7F00",
              "#377EB8",
              "#984EA3",
              "#4DAF4A",
              "#008280")
all_names = rownames(hm_dt)
all_color=rep("black",length(all_names))
for(i in 1:length(selected_names)){
  all_color[all_names %in% all_list[[selected_names[i]]]] = myninecolor[i]
}
table(all_color)


o1pgw = hclust(as.dist(1-cor(t(hm_dt)))) # about genes
o2pgw = seriate(as.dist(1-cor(hm_dt)), method = "GW") #columns
summary(o1pgw)
cutree(o1pgw,4) %>% {names(.)[.==4]} %>% {c("CD3E", "CD8A", "CD1A") %in% .}
cutree(o1pgw,4) %>% {names(.)[.==2]} %>% {table(h_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION %in% .)}
cutree(o1pgw,4) %>% {names(.)[.==2]} %>% {.[.%in%h_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION]}
cutree(o1pgw,4) %>% {names(.)[.==1]} %>% {.[. %in% h_list$HALLMARK_TNFA_SIGNALING_VIA_NFKB]}
cutree(o1pgw,4) %>% {names(.)[.==1]} %>% {.[. %in% h_list$HALLMARK_INTERFERON_GAMMA_RESPONSE]}
cutree(o1pgw,4) %>% {names(.)[.==1]} %>% {.[. %in% h_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION]}
cutree(o1pgw,4) %>% {names(.)[.==2]} %>% {.[. %in% h_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION]}
cutree(o1pgw,4) %>% {names(.)[.==3]} %>% {.[. %in% h_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION]}
cutree(o1pgw,4) %>% {names(.)[.==4]} %>% {.[. %in% h_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION]}

cutree(o1pgw,4) %>% {names(.)[.==1]} %>% {.[. %in% c5_list$GO_OXIDATIVE_PHOSPHORYLATION]}
cutree(o1pgw,4) %>% {names(.)[.==2]} %>% {.[. %in% c5_list$GO_OXIDATIVE_PHOSPHORYLATION]}
cutree(o1pgw,4) %>% {names(.)[.==3]} %>% {.[. %in% c5_list$GO_OXIDATIVE_PHOSPHORYLATION]}
cutree(o1pgw,4) %>% {names(.)[.==4]} %>% {.[. %in% c5_list$GO_OXIDATIVE_PHOSPHORYLATION]}
table(c5_list$GO_OXIDATIVE_PHOSPHORYLATION %in% rownames(hm_dt) )
table(h_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION %in% rownames(hm_dt) )


cutree(o1pgw,4) %>% {names(.)[.==1]} %>% {.[. %in% c5_list$GO_LINOLEIC_ACID_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==2]} %>% {.[. %in% c5_list$GO_LINOLEIC_ACID_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==3]} %>% {.[. %in% c5_list$GO_LINOLEIC_ACID_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==4]} %>% {.[. %in% c5_list$GO_LINOLEIC_ACID_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==1]} %>% {.[. %in% c5_list$GO_DRUG_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==2]} %>% {.[. %in% c5_list$GO_DRUG_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==3]} %>% {.[. %in% c5_list$GO_DRUG_METABOLIC_PROCESS]}
cutree(o1pgw,4) %>% {names(.)[.==4]} %>% {.[. %in% c5_list$GO_DRUG_METABOLIC_PROCESS]}


cutree(o1pgw,4) %>% {names(.)[.==3]} %>% cat(sep="\n")


table(c5_list$GO_DRUG_METABOLIC_PROCESS %in% rownames(hm_dt) )
table(h_list$HALLMARK_OXIDATIVE_PHOSPHORYLATION %in% rownames(hm_dt) )


hm_dt[1:3,1:3] 

o1pgw$labels[o1pgw$order]


top_anno <- HeatmapAnnotation("Group"= hm_annot_dt$GTF2I_status2,
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
ha = rowAnnotation(foo = anno_mark(at = which(all_color != "black"),
                                   labels = all_names[all_color != "black"], 
                                   lines_gp = gpar(col=all_color[which(all_color != "black")]),
                                   labels_gp = gpar(col=all_color[which(all_color != "black")],cex=1)))
# lgd = Legend(labels = c("Keratinization", "TCR signaling", "TGFβ signaling"), title = "Gene ontology", legend_gp = gpar(fill = myninecolor))
lgd = Legend(labels = selected_names, title = "Gene ontology", legend_gp = gpar(col = myninecolor,bg="white"),
             labels_gp = gpar(col = myninecolor),
             type = "lines")

lgd2 = Legend(col_fun = exp_pal, title = "Normalized gene expression")

lgd12=packLegend(lgd,lgd2)
#draw(lgd12)
get_order(o1pgw)
get_order(o1pgw)
get_order(o1pgw)
#plot(o1pgw)
table(cutree(o1pgw,4)) %>% cumsum
table(get_order(o1pgw)[c(which(cutree(o1pgw,4)==1),
                     which(cutree(o1pgw,4)==4),
                     which(cutree(o1pgw,4)==2),
                     which(cutree(o1pgw,4)==3))] ==get_order(o1pgw))
class(o1pgw)
h1 <- Heatmap(hm_dt,
              # row_order = get_order(o1pgw)[c(1:270,2166:2500,271:1553,1554:2165)],
              # row_order = get_order(o1pgw)[c(which(cutree(o1pgw,4)==1),
              #                                which(cutree(o1pgw,4)==4),
              #                                which(cutree(o1pgw,4)==2),
              #                                which(cutree(o1pgw,4)==3))],
              width = unit(9,"cm"),height = unit(11,"cm"),
              show_row_dend = F,
              show_column_names = F,
              cluster_columns = as.dendrogram(o2pgw[[1]]),
              cluster_rows = rev(as.dendrogram(o1pgw)),
              top_annotation = top_anno, col= exp_pal, show_row_names = F,
              right_annotation = ha,show_heatmap_legend = F,
              row_split = 4,
              left_annotation = rowAnnotation(width = unit(4.5, "mm"),d=anno_block(
                gp = gpar(fill = c("#FF7F00","#008280","#EE0000","#377EB8")),
                labels = c("cluster1", "cluster2", "cluster3","cluster4"),
                labels_gp = gpar(col = "white", fontsize = 10)
              )),
              row_title_gp = gpar(col="#00000000"),
              row_gap = unit(0.5, "mm"))

draw(h1)
draw(h1, annotation_legend_list = lgd12)

#h1
#cluster1 = cluster1 in fig
#cluster2 = cluster3 in fig
#cluster3 = cluster4 in fig
#cluster4 = cluster2 in fig
#save gene list
if(F){
  row_order_list <- row_order(h1)
  lapply(row_order_list, length)
  #612, 270, 335, 1283
  #save clutser name after the current figure name
  cl1_genes <- rownames(hm_dt)[row_order_list[[1]]]
  cl2_genes <- rownames(hm_dt)[row_order_list[[4]]]
  cl3_genes <- rownames(hm_dt)[row_order_list[[2]]]
  cl4_genes <- rownames(hm_dt)[row_order_list[[3]]]
  
  dt1 <- tibble(cluster = "c1", gene_name = cl1_genes)
  dt2 <- tibble(cluster = "c2", gene_name = cl2_genes)
  dt3 <- tibble(cluster = "c3", gene_name = cl3_genes)
  dt4 <- tibble(cluster = "c4", gene_name = cl4_genes)
  
  gene_name_dt <- bind_rows(dt1, dt2) %>% bind_rows(dt3) %>%
    bind_rows(dt4)
  gene_name_dt %>% write_tsv("~/00_Project/01_thymoma/10_Final_data/02_metadata/gene_clusters_in_heatmap.tsv")
  
}


cairo_pdf("figures/heatmap_rna_expression.1.pdf",height = 15/2.54,width=25/2.54,pointsize = 12*0.7)
draw(h1, annotation_legend_list = lgd12)
dev.off()

