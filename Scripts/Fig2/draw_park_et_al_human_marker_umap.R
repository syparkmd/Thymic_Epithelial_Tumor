
library(tidyverse)
library(patchwork)
library(ggrastr)


park_etal <- read_csv("extracted_cell_data_inc_irs4.csv")
nrow(park_etal)
park_etal <- park_etal %>% filter(Anno_level_2 == "TEC")
nrow(park_etal)

max_val=6

# Define the color palette used in the original plot for cell type and age group
batchcolor <- c("#BEBEBE50", "#228a3bA0")

# Define marker genes
markers_human <- c('FOXN1', "PSMB11", "LGALS7", "ENPEP", "PAX1", "FEZF2", "AIRE", "CD80", "SH2D6", "L1CAM",
                   "OVOL3", "PROX1", "PDPN", "TNFRSF11A", "CCL21", "IVL", "CLDN3", "TGFB2", "SALL1", "BMP4", "TBX1", "IRX1")

# Create age groups
park_etal <- park_etal %>%
  mutate(
    AgeGroup = case_when(
      Age %in% c("7w", "8w", "9w", "10w", "11w", "12w", "13w", "14w", "16w", "17w") ~ "7w - 17w",
      Age %in% c("3m", "6m", "10m", "15m", "30m") ~ "3m - 30m",
      Age %in% c("13y", "24y", "35y") ~ "13y - 35y",
      TRUE ~ "Other"
    )
  )

# Function to create a consistent UMAP plot with common theme settings
create_umap_plot <- function(data, aes_mapping, title_text, color_scale = NULL, size_val = 1, alpha_val = 1) {
  p <- ggplot(data, aes_mapping) +
    geom_point_rast(size = size_val, shape = 21, stroke = 0,
                    dpi=300) + # Use shape 21 for fill only, no border
    labs(x = NULL, y = NULL) +
    coord_fixed(ratio = 1, xlim = c(6.7, 12), ylim = c(0.7, 8)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0, size = 11, face = "bold", margin = margin(t = 0, r = 0, b = 0, l = 0)),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      legend.position = "none"
    ) +
    ggtitle(title_text)
  
  if (!is.null(color_scale)) {
    p <- p + color_scale
  }
  return(p)
}

# 1. Plots for Anno_level_3 (Cell Types)
unique_anno_levels <- c("cTEC", "mTEC", "Epi_GCM2", "TEC(myo)", "TEC(neuro)")
plots_anno <- lapply(unique_anno_levels, function(level) {
  plot_data <- park_etal %>%
    mutate(color_group = factor(Anno_level_3 == level, levels = c(FALSE, TRUE)))
  
  create_umap_plot(
    data = plot_data,
    aes_mapping = aes(x = umap_1, y = umap_2, fill = color_group),
    title_text = level,
    color_scale = scale_fill_manual(values = batchcolor),
    size_val = 1.2
  )
})

# 2. Plots for AgeGroup
unique_age_groups <- c("7w - 17w", "3m - 30m", "13y - 35y")
plots_age <- lapply(unique_age_groups, function(group) {
  plot_data <- park_etal %>%
    mutate(color_group = factor(AgeGroup == group, levels = c(FALSE, TRUE)))
  
  create_umap_plot(
    data = plot_data,
    aes_mapping = aes(x = umap_1, y = umap_2, fill = color_group),
    title_text = group,
    color_scale = scale_fill_manual(values = batchcolor),
    size_val = 1.2
  )
})


# Plots for Marker Genes
gene_colors <- c("#29394c", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
                 "#ffff33", "#fdae61", "#f76234", "#e84053", "#ee165d")

plots_genes <- lapply(markers_human, function(gene) {
  # Sort the data by gene expression to ensure higher values are plotted on top,
  plot_data <- park_etal %>%
    dplyr::select(umap_1, umap_2, !!sym(gene)) %>%
    drop_na() %>%
    arrange(!!sym(gene)) %>%
    mutate(exp = case_when(!!sym(gene) > max_val ~ max_val,
                           TRUE ~ !!sym(gene)))
  
  
  p <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = exp)) +
    geom_point_rast(shape = 20, alpha = 0.8, size = 0.8, dpi=300) +
    labs(x = NULL, y = NULL) +
    coord_fixed(ratio = 1, xlim = c(6.7, 12), ylim = c(0.7, 8)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0, size = 11, face = "bold", margin = margin(t = 0, r = 0, b = 0, l = 0)),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      legend.position = "none"
    ) +
    scale_color_gradientn(colors = gene_colors, limits = c(0,max_val)) +
    ggtitle(gene)
  
  return(p)
})

# Create a color bar legend manually using `ggplot2`.
legend_data <- data.frame(
  x = c(0, max_val),
  y = c(0, 1)
)

legend_plot <- ggplot(legend_data, aes(x, y)) +
  geom_raster(aes(fill = x)) +
  scale_fill_gradientn(
    colors = gene_colors,
    name = expression(log[10]*"("*CPM~"x"~100*")"),
    limits = c(0, max_val)
  ) +
  theme_void() +
  labs(x = expression(log[10]*"("*CPM~"x"~100*")")) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.title = element_text(size = 9),
    plot.margin = margin(0.5, 0, 0.5, 0, "cm")
  )

#add plot_spacer to plots_genes
plots_age[[4]] <- plot_spacer()
plots_age[[5]] <- plot_spacer()
gg <- wrap_plots(c(plots_anno, plots_age, plots_genes), 
           ncol=5)

# Save the final plot
pdf("park_et_al_marker_gene_expression_ggplot_v2.pdf",
    width=12, height=12)
print(gg)
dev.off()

#save legend
pdf("park_et_al_marker_gene_expression_ggplot_v2.legend.pdf",
     width=3, height=3)
print(legend_plot)
dev.off()


#draw IRS4
if(F){
  irs4_data <- park_etal %>%
    dplyr::select(umap_1, umap_2, IRS4) 
  irs4_data <- irs4_data %>% mutate(exp = case_when(IRS4 > max_val ~ max_val,
                                       TRUE ~ IRS4))
  
  g <- ggplot(irs4_data, aes(x = umap_1, y = umap_2, color = exp)) +
    geom_point_rast(shape = 20, alpha = 0.8, size = 0.8, dpi=300) +
    labs(x = NULL, y = NULL) +
    coord_fixed(ratio = 1, xlim = c(6.7, 12), ylim = c(0.7, 8)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0, size = 11, face = "bold", margin = margin(t = 0, r = 0, b = 0, l = 0)),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      legend.position = "none"
    ) +
    scale_color_gradientn(colors = gene_colors, limits = c(0, max_val)) +
    ggtitle("IRS4")
  
  pdf("park_et_al_irs4_expression_ggplot_v2.pdf",
      width=4, height=4)
  print(g)
  dev.off()
  
}

