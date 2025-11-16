#using NMF

library(tidyverse)
library(NMF)

#load meta table
meta_dt <- read_tsv("thymoma_meta_table.250826.tsv")
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']
CA_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'c']

#load biormart
bm_dt <- read_tsv('biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt') %>% 
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()

# load expression data
exp_dt <- read_tsv('thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')

# exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
#   filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only

#important thing is using log10(TPM+1) not 0.01 since you cannot use negative value in NMF
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+1)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only

#select 2500 variable protein-coding genes for NMF
hm_dt <- l10_exp_dt
hm_dt <- hm_dt %>% as.data.frame() %>% column_to_rownames('gene') 
hm_dt$var <- apply(hm_dt, 1, function(x) var(x))
hm_dt <- hm_dt %>% rownames_to_column('gene')
var_genes <- hm_dt %>% arrange(desc(var))%>%.$gene%>%.[1:2500] # top 2500 high variance genes, variance in log10 scale

exp_mx <- l10_exp_dt %>% as.data.frame() %>% column_to_rownames("gene") %>% as.matrix()
input_mx <- exp_mx[var_genes,]

estimate_k <- nmf(input_mx,
                  rank = 2:6,
                  method = 'brunet',
                  nrun = 30,
                  .options = list(parallel = 4),
                  seed = 123456)



pdf("NMF_results.pdf",
    width=8, height=8)
plot(estimate_k)
dev.off()


