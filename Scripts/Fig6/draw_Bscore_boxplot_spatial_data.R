#250929 drawing boxplot for spatial B score


library(tidyverse)
library(ggrastr)



dt <- data.table::fread("Yasumizu_Visium/adata_thy_obs.csv") %>% as_tibble()


sample_order <- dt %>% group_by(sample_id) %>%
  summarise(mean_B = mean(B_score)) %>%
  arrange(desc(mean_B)) %>% pull(sample_id)

g <- ggplot(dt, aes(x=sample_id, y=B_score))+
  geom_boxplot(outlier.size=-1)+
  geom_jitter_rast(width=0.2, height=0, alpha=0.1)+
  scale_x_discrete(limits=sample_order)+
  theme_cowplot()

pdf("Yasumizu_Visium/B_score_box_ggrast.pdf",
    width=8, height=4)
print(g)
dev.off()