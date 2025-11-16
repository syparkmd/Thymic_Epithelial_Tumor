
library(tidyverse)
library(gridExtra)
library(fgsea)
require(GSVA)
require(GSEABase)
library(grid)


out_dir = "pdf_from_R"
out_pdf_path="bulkRNA_gsea_250908.pdf"
out_tsv_path="gsea_results.tsv"

#data path
meta_path="thymoma_meta_table.250826.tsv"
exp_tpm_path="thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv"
mixcr_path='Clone3_Fraction0.01_summary_tbl.tsv'

#reference path
bm_path='biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt'
cgs_path="Cancer_gene_census_GRCh37_v89.tsv"
h_list_path='h.all.v6.2.symbols.gmt'
c5bp_list_path='c5.bp.v6.2.symbols.gmt'
kegg_list_path='c2.cp.kegg.v6.2.symbols.gmt'
reactome_list_path='c2.cp.reactome.v6.2.symbols.gmt'


# Helper function for cleaning up pathway names
clean_pathway_name <- function(x) {
  x %>%
    tolower() %>%
    str_remove("[^_]*_") %>%
    str_replace("rna", "RNA") %>%
    str_replace("dna", "DNA") %>%
    str_replace("_i_", "_I_") %>%
    str_replace("calcium", "Ca") %>%
    str_replace("dependent", "dep.") %>%
    str_replace("_iii_", "_III_") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# prep ----------------------------------------------------------------------------------------------------------------
meta_dt <- read_tsv(meta_path)

volc_dt <- read_tsv(exp_tpm_path) %>%
  mutate(across(-gene, ~log10(. + 0.01)))

# Filter metadata to get samples of interest
meta_filtered <- meta_dt %>%
  filter(!is.na(final_cellularity) & final_cellularity >= 0.2)
nrow(meta_filtered) #77
meta_filtered$GTF2I_status2 %>% table()
#c  m  w 
#12 49 16 
# m + w = 65 thymomas
 
# Define the groups of IDs
group_ids <- list(
  m = meta_filtered %>% filter(GTF2I_status2 == 'm') %>% pull(id),
  c = meta_filtered %>% filter(GTF2I_status2 == 'c') %>% pull(id),
  w = meta_filtered %>% filter(GTF2I_status2 == 'w') %>% pull(id),
  C = meta_filtered %>% filter(GTF2I_status2 != 'c') %>% pull(id)
)

# Calculate mean expression for each group, handling empty lists of IDs
meanexp <- map(group_ids, ~ {
  if (length(.x) > 0) {
    # If the group is not empty, calculate row means on the selected columns.
    # We use `select(all_of(.x))` for robustness, as single-column data frames
    # can be coerced to vectors by `[` and cause issues.
    rowMeans(volc_dt %>% dplyr::select(all_of(.x)), na.rm = TRUE)
  } else {
    # If the group is empty, return a vector of NAs with the correct length.
    # This prevents the "length must be the same" error later.
    rep(NA_real_, nrow(volc_dt))
  }
})

stats <- list(
  C = setNames(meanexp$c - meanexp$C, volc_dt$gene),
  mw = setNames(meanexp$w - meanexp$m, volc_dt$gene)
)

# Read GMT files and combine gene sets
h_list <- getGmt(h_list_path) %>% geneIds()
c5bp_list <- getGmt(c5bp_list_path) %>% geneIds()
kegg_list <- getGmt(kegg_list_path) %>% geneIds()
reactome_list <- getGmt(reactome_list_path) %>% geneIds()

all_list <- c(h_list, c5bp_list, kegg_list, reactome_list) %>%
  .[str_replace(names(.), "_.*", "") %in% c("GO", "HALLMARK", "KEGG", "REACTOME")]

length(all_list) #5346

# Run gsea
set.seed(1)
fgseaRes <- fgseaMultilevel(
  pathways = all_list,
  stats = stats$mw,
  minSize = 15,
  maxSize = 500,
  nproc = 10
) 
dim(fgseaRes) #4158 8

#collapse pathway
set.seed(1)
collapsedPathways_all_mw <- collapsePathways(
  fgseaRes[order(pval)][padj < 0.05],
  pathways= all_list,
  stats$mw
)

# Tidy selection of top pathways
topPathwaysUp <- fgseaRes %>%
  filter(pathway %in% collapsedPathways_all_mw$mainPathways & NES > 0 &
           padj < 0.05) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10) %>%
  pull(pathway)

topPathwaysDown <- fgseaRes %>%
  filter(pathway %in% collapsedPathways_all_mw$mainPathways & NES < 0 &
           padj < 0.05) %>%
  arrange(NES) %>%
  slice_head(n = 10) %>%
  pull(pathway)

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

resmat <- fgseaRes %>%
  filter(pathway %in% topPathways) %>%
  mutate(pathway = fct_reorder(pathway, -NES)) %>%
  mutate(pathway_clean = clean_pathway_name(pathway))

# draw -------
path_order <- resmat %>% arrange(NES) %>% pull(pathway_clean)
g <- ggplot(resmat, aes(x = NES, y = pathway_clean)) +
  geom_segment(aes(xend = 0, yend = pathway_clean), color = "grey50", size = 2, lineend = "round") +
  geom_point(aes(fill = -log10(padj)), shape = 21, size = 4) +
  scale_y_discrete(limits = path_order)+
  scale_fill_viridis_c(
    name = expression(-log[10] * "q"),
    direction = -1,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  labs(
    title = "Top enriched gene sets",
    subtitle = "GTF2I-mutant ↔ wild-type",
    x = "Normalized Enrichment Score",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    legend.position = "right",
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.y = element_blank()
  )

pdf(out_pdf_path,
    width=6, height=6)
print(g)
dev.off()

#save table for Supplementary Table
out_dt <- fgseaRes %>% rowwise() %>%
  mutate(group = unlist(strsplit(as.character(pathway), "_"))[1]) %>%
  ungroup() 
out_dt <- out_dt %>% rowwise() %>%
  mutate(genes = paste(leadingEdge, collapse=',')) %>% ungroup()
s_out_dt <- out_dt %>% dplyr::select(pathway, pval, padj, ES, NES, size, genes, group)
nrow(s_out_dt) #4158
fs_out_dt <- s_out_dt %>% filter(padj < 0.05)
nrow(fs_out_dt) #197
#ordering
fs_out_dt <- fs_out_dt %>% arrange(desc(NES))
#scientific notation
fs_out_dt$pval <- formatC(fs_out_dt$pval, format = "E", digits = 2)
fs_out_dt$padj <- formatC(fs_out_dt$padj, format = "E", digits = 2)
fs_out_dt$ES <- round(fs_out_dt$ES, 3)
fs_out_dt$NES <- round(fs_out_dt$NES, 3)
fs_out_dt %>% write_tsv()
