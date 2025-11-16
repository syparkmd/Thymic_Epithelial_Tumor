#run infer cnv using publec scRNAseq

library(infercnv)

#path
count_matrix_path= "Yasumizu_2024/tables_from_h5ad/counts_layer.mtx"
X_genes_path="Yasumizu_2024/tables_from_h5ad/X_matrix_genes.tsv"
X_cells_path="Yasumizu_2024/tables_from_h5ad/X_matrix_cells.tsv"
obs_csv_path="Yasumizu_2024/tables_from_h5ad/obs.csv"

#ref path
gene_range_path="biomart_human_gene_range.txt"

#load count matrix
ct_matrix_sparse <- readMM(count_matrix_path)
X_genes <- read.delim(X_genes_path, header = FALSE, stringsAsFactors = FALSE)[,1]
X_cells <- read.delim(X_cells_path, header = FALSE, stringsAsFactors = FALSE)[,1]
rownames(ct_matrix_sparse) <- X_genes
colnames(ct_matrix_sparse) <- X_cells

#load obs
obs_dt <- data.table::fread(obs_csv_path)
obs_dt <- obs_dt %>% as_tibble()
#filter only thymus samples
f_obs_dt <- obs_dt %>% filter(site == "Thymus")
#make column tn_group 
f_obs_dt <- f_obs_dt %>% mutate(tn_group = 
                      case_when(cluster_L1 == "TEC" ~ "tumor",
                                TRUE ~ "normal"))
f_obs_dt %>% filter(tn_group == "tumor") %>%
  group_by(sample) %>% count()
#filter only tumor cells
tumor_dt <- f_obs_dt %>% filter(tn_group == "tumor") 
#annotate n_total (per sample)
tumor_dt <- left_join(tumor_dt, tumor_dt %>% group_by(sample) %>% summarise(n_total =n()) %>%
  ungroup())
nrow(tumor_dt) #5446
#random sampling per sample max 100
set.seed(123)
r_tumor_dt <- tumor_dt %>% group_by(sample) %>%
  sample_n(min(100, n_total)) %>% ungroup()
nrow(r_tumor_dt) #758

#filter only normal cells
normal_dt <- f_obs_dt %>% filter(tn_group == "normal")
nrow(normal_dt) #77877
#random sampling per sample max 50
set.seed(123)
r_normal_dt <- normal_dt %>% group_by(sample) %>%
  sample_n(50) %>% ungroup()
nrow(r_normal_dt) #600
#bind tumor and normal dt
tn_dt <- bind_rows(r_tumor_dt, r_normal_dt)
#make tn_group2
tn_dt <- tn_dt %>% mutate(tn_group2 = 
                   case_when(tn_group == "tumor" ~ sample,
                             TRUE ~ tn_group)) 

#select only target cells in the exp matrix
target_ids <- tn_dt$V1
s_ct_mx <- ct_matrix_sparse[,target_ids] %>% as.matrix()
dim(s_ct_mx) #27631 1358


#make group df
annot_df <- tn_dt %>%
  select(V1, tn_group2) %>%
  as.data.frame() %>% column_to_rownames("V1")
#check
all(rownames(annot_df) == colnames(s_ct_mx))

#make gene_order_file
gene_range_dt <- read_tsv(gene_range_path)
colnames(gene_range_dt) <- gsub("/", "_", gsub(" ", "_", colnames(gene_range_dt)))
s_gr_dt <- gene_range_dt %>% select(Gene_name, Chromosome_scaffold_name,
                         `Gene_start_(bp)`,`Gene_end_(bp)`) %>%
  unique() %>%
  dplyr::rename(gene_start = `Gene_start_(bp)`,
                gene_end = `Gene_end_(bp)`, 
                chrom = Chromosome_scaffold_name)
#filter genes only in chrom1:22
#there are many genes the location is overlapped in chrX and Y
fs_gr_dt <- s_gr_dt %>% filter(chrom %in% c(1:22)) 
nrow(fs_gr_dt) #75022

#filter only intersect genes and dedup
itx_genes <- intersect(rownames(s_ct_mx), fs_gr_dt$Gene_name)
length(itx_genes) #18817
fs_ct_mx <- s_ct_mx[itx_genes,]
ffs_gr_dt <- fs_gr_dt %>% filter(Gene_name %in% itx_genes) %>% unique()
nrow(ffs_gr_dt) #18843

dup_genes <- ffs_gr_dt %>% 
  group_by(Gene_name) %>% count() %>% filter(n>1) %>%
  pull(Gene_name)
length(dup_genes) #26

ffs_gr_dt %>% filter(Gene_name %in% dup_genes) %>% arrange(Gene_name, chrom, gene_start) %>%
  print(n=60)
#all have same chrom but little different start and end

#make single position and gene_order_df
gene_order_df = ffs_gr_dt %>% group_by(chrom, Gene_name) %>%
  summarise(start = min(gene_start),
            end = max(gene_end)) %>% ungroup() %>%
  as.data.frame() %>% column_to_rownames("Gene_name")
gene_order_df = gene_order_df[itx_genes,]
all(rownames(gene_order_df) == rownames(fs_ct_mx)) #TRUE

#make infercnv object
ifc_obj <- CreateInfercnvObject(
  raw_counts_matrix = s_ct_mx,
  gene_order_file = gene_order_df,
  annotations_file =  annot_df,
  ref_group_names = "normal")

#Run InferCNV
infercnv_results <- run(ifc_obj,
                        cutoff = 1,  # Set the cutoff for CNV detection
                        out_dir = "inferCNV_2",  # Output directory
                        cluster_by_groups = TRUE,  # Cluster by cell types
                        denoise = TRUE,  # Denoise the data
                        plot_steps = TRUE,# Plot the steps
                        HMM = FALSE)  

# Plot the results
plot_cnv(infercnv_results, 
         out_dir = "inferCNV_2",
         output_format = "pdf")

