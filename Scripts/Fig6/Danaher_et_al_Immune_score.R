
#Danaher et al. Journal for ImmunoTherapy of Cancer (2017)5:18
#### make markers using only tumor samples
#We refined markers using our 137 samples

library(tidyverse)

tpm_dt <- read_tsv('thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')
exp_dt <- tpm_dt %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% +0.01 %>% log10() %>% as.data.frame() %>% rownames_to_column('gene') %>% as.tibble()

Bcells <- c("BLK","CD19","FCRL2","MS4A1","TNFRSF17","TCL1A","SPIB","PNOC")
#FAM30A is not exist in dataset
length(Bcells)
tmp_dt <- exp_dt %>% filter(gene %in% Bcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
library(PerformanceAnalytics)
chart.Correlation(tmp_dt)

Cytotoxic <- c("PRF1","GZMA","GZMB","NKG7","GZMH","KLRK1","KLRB1","KLRD1","CTSW","GNLY")
length(Cytotoxic)
tmp_dt <- exp_dt %>% filter(gene %in% Cytotoxic) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

#since DC marker genes were few I aded other genes using Davoli et al. Science 2017 paper supple table S4d
#DC <- c("CCL13","CD209","HSD11B1")
DC <- c("CCL13","CD209","HSD11B1", "NR4A3","HAVCR2","KMO","DNASE1L3","ANPEP","CXCL16")
length(DC)
tmp_dt <- exp_dt %>% filter(gene %in% DC) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)
corrplot::corrplot(cor(tmp_dt))
#All genes show some level of positive correlation. All genes were included


Exhausted <- c("LAG3","TIGIT","HAVCR2","CTLA4")
#CD244 and PTGER4 were excluded
#Exahaustd markers were edited by article search
length(Exhausted)
tmp_dt <- exp_dt %>% filter(gene %in% Exhausted) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

Macrophages <- c("CD68","CD163","MS4A4A", "CD84")
length(Macrophages)
tmp_dt <- exp_dt %>% filter(gene %in% Macrophages) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

Mast <- c("TPSB2","TPSAB1","CPA3","MS4A2","HDC")
length(Mast)
tmp_dt <- exp_dt %>% filter(gene %in% Mast) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

Neutrophils <- c("FPR1","SIGLEC5","CSF3R","FCAR","FCGR3B")
#S100A12 and CEACAM3 was excluded
length(Neutrophils)
tmp_dt <- exp_dt %>% filter(gene %in% Neutrophils) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)


NKCD56dim <- c("KIR2DL3","KIR3DL1","KIR3DL2")
#IL21R was excluded
length(NKCD56dim)
tmp_dt <- exp_dt %>% filter(gene %in% NKCD56dim) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

NKcells <- c("XCL1","XCL2","NCR1")
length(NKcells)
tmp_dt <- exp_dt %>% filter(gene %in% NKcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

Tcells <- c("CD6","CD3D","CD3E","SH2D1A","TRAT1","CD3G")
length(Tcells)
tmp_dt <- exp_dt %>% filter(gene %in% Tcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

Earlythymo <- c("IL2RA","KIT")
length(Earlythymo)
tmp_dt <- exp_dt %>% filter(gene %in% Earlythymo) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)

Latethymo <-c("CD1A","CD1B","DNTT")
length(Latethymo)
tmp_dt <- exp_dt %>% filter(gene %in% Latethymo) %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix() %>% t()
ncol(tmp_dt)
chart.Correlation(tmp_dt)


CD45 <- "PTPRC"
Th1cells <- "TBX21"
Treg <- "FOXP3"
CD8Tcells <- c("CD8A","CD8B")



###calculate score in each samples
#in final version scoreo was calculated using GSVA package (other script)
sample_list <- colnames(exp_dt)[2:ncol(exp_dt)]
score_tbl <- read.table(text='', col.names = sample_list)
score_tbl["B_cells",] <- exp_dt %>% filter(gene %in% Bcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Cytotoxic_cells",] <- exp_dt %>% filter(gene %in% Cytotoxic) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["DC",] <- exp_dt %>% filter(gene %in% DC) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Exhausted_CD8",] <- exp_dt %>% filter(gene %in% Exhausted) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Macrophages",] <- exp_dt %>% filter(gene %in% Macrophages) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Mast_cells",] <- exp_dt %>% filter(gene %in% Mast) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Neutrophils",] <- exp_dt %>% filter(gene %in% Neutrophils) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["NK_CD56dim_cells",] <- exp_dt %>% filter(gene %in% NKCD56dim) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["NK_cells",] <- exp_dt %>% filter(gene %in% NKcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["T_cells",] <- exp_dt %>% filter(gene %in% Tcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Early_thymocytes",] <- exp_dt %>% filter(gene %in% Earlythymo) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Late_thymocytes",] <- exp_dt %>% filter(gene %in% Latethymo) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["CD45",] <- exp_dt %>% filter(gene %in% CD45) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Th1_cells",] <- exp_dt %>% filter(gene %in% Th1cells) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["Treg",] <- exp_dt %>% filter(gene %in% Treg) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["CD8_T_cells",] <- exp_dt %>% filter(gene %in% CD8Tcells) %>% as.data.frame() %>% column_to_rownames('gene') %>% colMeans()
score_tbl["CD4_T_cells",] <- as.numeric(score_tbl["T_cells",]) - as.numeric(score_tbl["CD8_T_cells",])
colnames(score_tbl) <- gsub('\\.','-',colnames(score_tbl))
score_tbl <- score_tbl %>% as.data.frame() %>% rownames_to_column('cell_type') %>% as.tibble()
write_tsv(score_tbl, './13_Danaher_immune_cell/thym_137s_Danaher_score_250818.tsv')
