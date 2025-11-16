# plotting spatial data

hm_dt <- read_tsv("Yasumizu_Visium/heaatmap_table.tsv")
hm_mx <- hm_dt %>% as.data.frame() %>% column_to_rownames("...1") %>%
  as.matrix()



annot_dt <- read_tsv("Yasumizu_Visium/heaatmap_annot.tsv")
annot_df <- annot_dt %>% select(-`...1`)  %>% as.data.frame() %>%
  column_to_rownames("gene_name") 
annot_df <- annot_df[rownames(hm_mx),,drop=F]


left_anno <- ComplexHeatmap::rowAnnotation(df = annot_df,
                                           col= list(cluster = c("c1" ='#f67c32',
                                                                 "c2" ='#407cac')))
hm <- ComplexHeatmap::Heatmap(hm_mx, left_annotation = left_anno,
                        show_row_names = F, 
                        show_column_names = T,
                        clustering_distance_columns = "pearson",
                        clustering_distance_rows = "pearson",
                        clustering_method_columns = "average",
                        use_raster = T)

pdf("Yasumizu_Visium/heaatmap_raster.pdf",
    width=5, height=6)
print(hm)
dev.off()