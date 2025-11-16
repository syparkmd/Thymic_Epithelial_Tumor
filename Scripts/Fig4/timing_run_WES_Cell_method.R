# annotate mutcn to wes for timing

library(tidyverse)
library(cowplot)

out_dir = "timing_analysis_wes"

#data path
meta_path="thymoma_meta_table.250826.tsv"
pm_dir = "01_annotated_file"
seqz_dir = "Final_selected_solutions"

#ref path
cyto_path="hg19_cytoBand.txt"
bait_path="S04380219_Regions_cut_v37.bed"

#meta table load
meta_dt <- read_tsv(meta_path)



#meta data update for a sample
if(T){
  meta_dt$final_cellularity[meta_dt$id == "TCGA-ZB-A96F"] <- 0.35
  meta_dt$final_ploidy[meta_dt$id == "TCGA-ZB-A96F"]  <- 3.1
  meta_dt$final_selection[meta_dt$id == "TCGA-ZB-A96F"]   <- "sqz3_manual"
}


#check
nrow(meta_dt) #137
meta_dt %>% group_by(GTF2I_status2) %>% count() 
#m                71
#2 w                53
#3 c                13

#load cyto path to check centromere position
cyto_dt <- read_tsv(cyto_path, col_names=c("chrom","start","end","band","stain"))
acen_dt <- cyto_dt %>% filter(stain=="acen")
m_acen_dt <- acen_dt %>% group_by(chrom) %>% summarise(start=min(start), end=max(end)) %>%
  dplyr::rename(chromosome = chrom) %>%
  mutate(mid_acen_pos = (start+end)/2)


#make sh scripts to annotate early/late probability

if(F){
 
  #filter only CN available
  f_meta_dt <- meta_dt %>% filter(final_selection != "low")
  
  nrow(f_meta_dt) #89
  f_meta_dt$GTF2I_status2 %>% table()
  #c  m  w 
  #12 54 23 
  
  #annot seqz seg path
  merged_dt <- tibble()
  for(i in 1:nrow(f_meta_dt)){
    print(i)
    t_id <- f_meta_dt$id[i]
    t_sqz_id <- list.files(seqz_dir, pattern=t_id)
    t_seg_path = paste0(seqz_dir, '/', t_sqz_id, '/',
                        t_sqz_id, '_segments.txt')
    t_out <- tibble(id = t_id,
                    seqz_id = t_sqz_id,
                    seg_path = t_seg_path)
    merged_dt <- bind_rows(merged_dt, t_out)
  }
  #run_annot_seqzCN_ind/snv_2505.sh
  merged_dt <- merged_dt %>% mutate(ind_seg_cmd = 
                                      case_when(grepl("SNU", id) ~
                         paste0("python SNV_seqz_CN_annot.py ",
                                id, ".indel.varscan2_strelka2_union.vcf.readinfo.readc.rasm.thyN13.snuN30.bgiN24.anv.cap.fi.np_fi.man_edit.sorted.add_fi ",
                                seg_path),
                         TRUE ~ paste0("python SNV_seqz_CN_annot.py ",
                                       id, ".indel.varscan2_strelka2_union.vcf.readinfo.readc.rasm.snuN30.bgiN24.thyN112.rasm_re.anv.cap.fi.np_fi.man_edit.sorted.add_fi ",
                                       seg_path)))
  merged_dt$snv_seg_cmd <- gsub("indel", "snv", merged_dt$ind_seg_cmd)
  merged_dt %>% select(ind_seg_cmd) %>%
    write_tsv("16_point_mutation_final_call/01_annotated_file/run_annot_seqzCN_ind_2505.sh",
              col_names = F)
  merged_dt %>% select(snv_seg_cmd) %>%
    write_tsv("16_point_mutation_final_call/01_annotated_file/run_annot_seqzCN_snv_2505.sh",
              col_names = F)
  
  #run_annot_mutCN_ind/snv_2505.sh
  merged_dt <- left_join(merged_dt, f_meta_dt %>% select(id, final_cellularity))
  merged_dt <- merged_dt %>% mutate(ind_mutcn_cmd = 
                         case_when(grepl("SNU", id) ~
                                     paste0("python2.7 SNV_mutCN_CCF_annot.py ", id,
                                            ".indel.varscan2_strelka2_union.vcf.readinfo.readc.rasm.thyN13.snuN30.bgiN24.anv.cap.fi.np_fi.man_edit.sorted.add_fi.seqzcn ",
                                            final_cellularity),
                                   TRUE ~ 
                                     paste0("python2.7 SNV_mutCN_CCF_annot.py ", id,
                                            ".indel.varscan2_strelka2_union.vcf.readinfo.readc.rasm.snuN30.bgiN24.thyN112.rasm_re.anv.cap.fi.np_fi.man_edit.sorted.add_fi.seqzcn ",
                                            final_cellularity)))
  merged_dt$snv_mutcn_cmd <- gsub("indel", "snv", merged_dt$ind_mutcn_cmd)
  merged_dt %>% select(ind_mutcn_cmd) %>%
    write_tsv("16_point_mutation_final_call/01_annotated_file/run_annot_mutCN_ind_2505.sh",
              col_names = F)
  merged_dt %>% select(snv_mutcn_cmd) %>%
    write_tsv("16_point_mutation_final_call/01_annotated_file/run_annot_mutCN_snv_2505.sh",
              col_names = F)
  
  
  #run_annot_timingProba_2505.sh
  merged_dt <- merged_dt %>% mutate(ind_proba_cmd = 
                                      case_when(grepl("SNU", id) ~
                                                  paste0("python snv_timing_probability_1-4.py ", id,
                                                         ".indel.varscan2_strelka2_union.vcf.readinfo.readc.rasm.thyN13.snuN30.bgiN24.anv.cap.fi.np_fi.man_edit.sorted.add_fi.seqzcn.scF ",
                                                         final_cellularity),
                                                TRUE ~ 
                                                  paste0("python2.7 snv_timing_probability_1-4.py ", id,
                                                         ".indel.varscan2_strelka2_union.vcf.readinfo.readc.rasm.snuN30.bgiN24.thyN112.rasm_re.anv.cap.fi.np_fi.man_edit.sorted.add_fi.seqzcn.scF ",
                                                         final_cellularity)))
  merged_dt$snv_proba_cmd <- gsub("indel", "snv", merged_dt$ind_proba_cmd)
  merged_dt %>% select(ind_proba_cmd) %>%
    write_tsv("16_point_mutation_final_call/01_annotated_file/run_annot_timingProba_ind_2505.sh",
              col_names = F)
  merged_dt %>% select(snv_proba_cmd) %>%
    write_tsv("16_point_mutation_final_call/01_annotated_file/run_annot_timingProba_snv_2505.sh",
              col_names = F)
  
}


#run shell scripts: done

#loop sample by sample
#Only samples with CNV analysis available
f_meta_dt <- meta_dt %>% filter(final_selection != "low")
nrow(f_meta_dt) #89
f_meta_dt %>% filter(final_cellularity >= 0.3 & GTF2I_status2 == "w") %>%
   nrow() #12

if(file.exists(paste0(out_dir,'/chrom_level_early_late_muts_v2.tsv'))==F){
  merged_dt <- tibble()
  for(i in 1:nrow(f_meta_dt)){
    print(i)
    t_id = f_meta_dt$id[i]
    print(t_id)
    t_tcf = f_meta_dt$final_cellularity[i]
    print(t_tcf)

    #load pm dt
    t_snv_path = list.files(pm_dir, pattern = paste0(t_id,'\\.snv.*bino1-4P'),
                            full.names = T)
    t_snv_dt <- read_tsv(t_snv_path, col_types = cols(`#CHROM` = 'c',
                                                      majCN = 'n',
                                                      minCN = 'n',
                                                      pEarly ='n',
                                                      `pLate/minor` = 'n',
                                                      pClonal = 'n',
                                                      pSubclonal = 'n'))
    
    t_amp_dt <- t_snv_dt %>% filter(majCN >=2)
    #assume same number of mutation exist in minor allele
    t_amp_dt <- t_amp_dt %>% mutate(p_early_adj = case_when(majCN == minCN  ~ pEarly,
                                                            minCN == 0 ~ pEarly,
                                                            majCN > minCN ~ pEarly*2))
    t_chr_dt <- t_amp_dt %>% group_by(`#CHROM`) %>%
      summarise(n_mut = n(), 
                p_clonal_ct = sum(pClonal),
                p_early_ct = sum(pEarly),
                p_early_adj_ct = sum(p_early_adj),
                p_late_ct = sum(`pLate/minor`),
                p_subclonal_ct = sum(pSubclonal)) %>% ungroup()
    
    t_chr_dt <- t_chr_dt %>% rowwise() %>%
      mutate(p_late_adj_ct = max(p_clonal_ct - p_early_adj_ct, 0)) %>%
      ungroup() %>% mutate(sample_id = t_id)
    merged_dt <- bind_rows(merged_dt, t_chr_dt)
  }
  merged_dt %>% write_tsv(paste0(out_dir,'/chrom_level_early_late_muts_v2.tsv'))
}else{
  merged_dt <- read_tsv(paste0(out_dir,'/chrom_level_early_late_muts_v2.tsv'))
}

if(T){
  #add histologic_type GTF2I_status2
  merged_dt <- left_join(merged_dt, f_meta_dt %>% select(id, histologic_type, GTF2I_status2,
                                                         final_cellularity, final_ploidy) %>%
                           dplyr::rename(sample_id= id))
  merged_dt %>% filter(GTF2I_status2=='w') %>% arrange(desc(n_mut))
  merged_dt %>% filter(GTF2I_status2=='w' & p_clonal_ct >=3 &
                         final_cellularity >= 0.3) %>% nrow()
  
  #only filter in chroms in GTF2I WT and No. of muts >=4 & final_cellularity >= 0.3
  f_merged_dt <- merged_dt %>% filter(GTF2I_status2=='w' & n_mut >=4 &
                                        final_cellularity >= 0.3)
  nrow(f_merged_dt) #10
  #annot p_early_prop
  f_merged_dt <- f_merged_dt %>% rowwise() %>%
    mutate(p_early_prop = min(1, p_early_adj_ct/n_mut)) %>%
    ungroup()
  f_merged_dt <- f_merged_dt %>% mutate(p_late_prop = 1-p_early_prop)
  #annot amp_id
  f_merged_dt <- f_merged_dt %>% mutate(amp_id = paste0(sample_id, '.chr', `#CHROM`))
}


#plotting
if(F){
  x_order <- f_merged_dt %>% arrange(p_early_prop) %>% pull(amp_id)
  plot_dt <- f_merged_dt %>% select(amp_id, p_early_prop, p_late_prop)  %>%
    gather(-amp_id, key="group", value="prop")
  plot_dt$group <- factor(plot_dt$group, levels = c("p_late_prop", "p_early_prop"))
  g <- ggplot(plot_dt, aes(x=amp_id, y=prop))+
    geom_bar(stat="identity", aes(fill=group))+
    scale_x_discrete(limits=x_order)+
    scale_fill_manual(labels = c("p_early_prop" = "Pre-amplification",
                                 "p_late_prop" = "Post-amplification"),
                      values = c("p_early_prop" = "#33ba8a",
                                 "p_late_prop" = "#ea833a"))+
    labs(x="", y="Proportion")+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  
  pdf(paste0(out_dir,'/early_late_prop.pdf'),
      width=9, height=4)
  print(g)
  dev.off()
}


#example plotting.
#early example: TCGA-XU-A930, chr1 and chr2
if(F){
  t_id = "TCGA-XU-A930"
  t_sqz_id = list.files(seqz_dir, t_id)
  t_seg_path=paste0(seqz_dir, '/', t_sqz_id, '/',
                    t_sqz_id, '_segments.txt')
  t_snv_path=list.files(pm_dir, pattern = paste0(t_id,'\\.snv.*bino1-4P'),
                        full.names = T)
  #load seg dt
  t_seg_dt <- read_tsv(t_seg_path)
  #make absolute start and end position
  tmp_dt <- t_seg_dt %>% group_by(chromosome) %>% summarise(chrom_end = max(end.pos)) %>%
    ungroup()
  tmp_dt$chromosome <- factor(tmp_dt$chromosome,
                              levels = c(1:22, "X"))
  tmp_dt <- tmp_dt %>% arrange(chromosome)
  tmp_dt <- tmp_dt %>% mutate(lag_end = lag(chrom_end)) 
  tmp_dt$lag_end[is.na(tmp_dt$lag_end)] <- 0
  tmp_dt$cumsum <- cumsum(tmp_dt$lag_end)
  t_seg_dt <- left_join(t_seg_dt, tmp_dt %>% select(chromosome, cumsum))
  t_seg_dt <- t_seg_dt %>% mutate(abs_start = cumsum + start.pos,
                      abs_end = cumsum+ end.pos)
  #adjust
  t_seg_dt <- t_seg_dt %>% mutate(majCN = case_when(A == B ~ A+0.05,
                                        TRUE ~ A),
                      minCN = case_when(A==B ~ B-0.05,
                                        TRUE ~ B))
  #filter out short seg
  ft_seg_dt <- t_seg_dt %>% filter(end.pos - start.pos >= 3e6)
  #load snv dt
  t_snv_dt <- read_tsv(t_snv_path, col_types = cols(pEarly='n'))
  t_snv_dt <- left_join(t_snv_dt, tmp_dt %>% select(chromosome, cumsum) %>%
              dplyr::rename(`#CHROM` = chromosome))
  t_snv_dt <- t_snv_dt %>% mutate(abs_pos = cumsum + POS)
  t_snv_dt <- t_snv_dt %>% mutate(early_group = case_when(pEarly > 0.5 ~ "early",
                                              TRUE ~ "late"))
  
  #chr1
  fft_seg_dt <- ft_seg_dt %>% filter(chromosome == "1")
  ft_snv_dt <- t_snv_dt %>% filter(`#CHROM` == '1')
  g1 <- ggplot(fft_seg_dt)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                 y = majCN, yend=majCN), color="#e26733",
                 size=4)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = minCN, yend=minCN),
                 color="#1f449f",
                 size=4)+
    geom_point(data=ft_snv_dt, aes(x=POS, y=mutCN, color=pEarly),
               size=3)+
    scale_color_gradientn(colors = c("#b6b7b9", "#ed2c2c"),
                          values = c(0, 1),
                          limits = c(0, 1), 
                          na.value = "#b6b7b9")+
    scale_y_continuous(limits=c(0, 4.2), expand = c(0,0))+
    theme_cowplot()+theme(legend.position= "none")+
    labs(x="", y="Copy number")+
    ggtitle("chr1")
    
  
  #chr2
  fft_seg_dt <- ft_seg_dt %>% filter(chromosome == "2")
  ft_snv_dt <- t_snv_dt %>% filter(`#CHROM` == '2')
  g2 <- ggplot(fft_seg_dt)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = majCN, yend=majCN), color="#e26733",
                 size=4)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = minCN, yend=minCN),
                 color="#1f449f",
                 size=4)+
    geom_point(data=ft_snv_dt, aes(x=POS, y=mutCN, color=pEarly),
               size=3)+
    scale_color_gradientn(colors = c("#b6b7b9", "#ed2c2c"),
                          values = c(0, 1),
                          limits = c(0, 1))+
    scale_y_continuous(limits=c(0, 4.2), expand = c(0,0))+
    theme_cowplot()+theme(legend.position = "none")+
    labs(x="", y="Copy number")+
    ggtitle("chr2")
  
  #chr3
  fft_seg_dt <- ft_seg_dt %>% filter(chromosome == "3")
  ft_snv_dt <- t_snv_dt %>% filter(`#CHROM` == '3')
  g3 <- ggplot(fft_seg_dt)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = majCN, yend=majCN), color="#e26733",
                 size=4)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = minCN, yend=minCN),
                 color="#1f449f",
                 size=4)+
    geom_point(data=ft_snv_dt, aes(x=POS, y=mutCN, color=pEarly),
               size=3)+
    scale_color_gradientn(colors = c("#b6b7b9", "#ed2c2c"),
                          values = c(0, 1),
                          limits = c(0, 1))+
    scale_y_continuous(limits=c(0, 4.2), expand = c(0,0))+
    theme_cowplot()+theme(legend.position = "none")+
    labs(x="", y="Copy number")+
    ggtitle("chr3")
  
  gg <- plot_grid(g1, g2, g3, nrow=1, rel_widths= c(1,1,0.8)) 
  
  pdf(paste0(out_dir,'/TCGA-XU-A930_chr1to3.pdf'),
      width=9, height=4)
  print(gg)
  dev.off()
  
}


#relatively late example:SNU_18_C chr5
if(F){
  t_id = "SNU_18_C"
  t_sqz_id = list.files(seqz_dir, t_id)
  t_seg_path=paste0(seqz_dir, '/', t_sqz_id, '/',
                    t_sqz_id, '_segments.txt')
  t_snv_path=list.files(pm_dir, pattern = paste0(t_id,'\\.snv.*bino1-4P'),
                        full.names = T)
  #load seg dt
  t_seg_dt <- read_tsv(t_seg_path)
  #make absolute start and end position
  tmp_dt <- t_seg_dt %>% group_by(chromosome) %>% summarise(chrom_end = max(end.pos)) %>%
    ungroup()
  tmp_dt$chromosome <- factor(tmp_dt$chromosome,
                              levels = c(1:22, "X"))
  tmp_dt <- tmp_dt %>% arrange(chromosome)
  tmp_dt <- tmp_dt %>% mutate(lag_end = lag(chrom_end)) 
  tmp_dt$lag_end[is.na(tmp_dt$lag_end)] <- 0
  tmp_dt$cumsum <- cumsum(tmp_dt$lag_end)
  t_seg_dt <- left_join(t_seg_dt, tmp_dt %>% select(chromosome, cumsum))
  t_seg_dt <- t_seg_dt %>% mutate(abs_start = cumsum + start.pos,
                                  abs_end = cumsum+ end.pos)
  #adjust
  t_seg_dt <- t_seg_dt %>% mutate(majCN = case_when(A == B ~ A+0.05,
                                                    TRUE ~ A),
                                  minCN = case_when(A==B ~ B-0.05,
                                                    TRUE ~ B))
  #filter out short seg
  ft_seg_dt <- t_seg_dt %>% filter(end.pos - start.pos >= 3e6)
  #load snv dt
  t_snv_dt <- read_tsv(t_snv_path, col_types = cols(`#CHROM` = 'c', pEarly='n'))
  t_snv_dt <- left_join(t_snv_dt, tmp_dt %>% select(chromosome, cumsum) %>%
                          dplyr::rename(`#CHROM` = chromosome))
  t_snv_dt <- t_snv_dt %>% mutate(abs_pos = cumsum + POS)
  t_snv_dt <- t_snv_dt %>% mutate(early_group = case_when(pEarly > 0.5 ~ "early",
                                                          TRUE ~ "late"))
  
  #chr4
  fft_seg_dt <- ft_seg_dt %>% filter(chromosome == "4")
  ft_snv_dt <- t_snv_dt %>% filter(`#CHROM` == '4')
  g1 <- ggplot(fft_seg_dt)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = majCN, yend=majCN), color="#e26733",
                 size=4)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = minCN, yend=minCN),
                 color="#1f449f",
                 size=4)+
    geom_point(data=ft_snv_dt, aes(x=POS, y=mutCN, color=pEarly),
               size=3)+
    scale_color_gradientn(colors = c("#b6b7b9", "#ed2c2c"),
                          values = c(0, 1),
                          limits = c(0, 1), 
                          na.value = "#b6b7b9")+
    scale_y_continuous(limits=c(0, 2.5), expand = c(0,0))+
    theme_cowplot()+theme(legend.position= "none")+
    labs(x="", y="Copy number")+
    ggtitle("chr4")
  
  #chr5
  fft_seg_dt <- ft_seg_dt %>% filter(chromosome == "5")
  ft_snv_dt <- t_snv_dt %>% filter(`#CHROM` == '5')
  g2 <- ggplot(fft_seg_dt)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = majCN, yend=majCN), color="#e26733",
                 size=4)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = minCN, yend=minCN),
                 color="#1f449f",
                 size=4)+
    geom_point(data=ft_snv_dt, aes(x=POS, y=mutCN, color=pEarly),
               size=3)+
    scale_color_gradientn(colors = c("#b6b7b9", "#ed2c2c"),
                          values = c(0, 1),
                          limits = c(0, 1), 
                          na.value = "#b6b7b9")+
    scale_y_continuous(limits=c(0, 2.5), expand = c(0,0))+
    theme_cowplot()+theme(legend.position= "none")+
    labs(x="", y="Copy number")+
    ggtitle("chr5")
  
  #chr6
  fft_seg_dt <- ft_seg_dt %>% filter(chromosome == "6")
  ft_snv_dt <- t_snv_dt %>% filter(`#CHROM` == '6')
  g3 <- ggplot(fft_seg_dt)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = majCN, yend=majCN), color="#e26733",
                 size=4)+
    geom_segment(aes(x= start.pos, xend= end.pos,
                     y = minCN, yend=minCN),
                 color="#1f449f",
                 size=4)+
    geom_point(data=ft_snv_dt, aes(x=POS, y=mutCN, color=pEarly),
               size=3)+
    scale_color_gradientn(colors = c("#b6b7b9", "#ed2c2c"),
                          values = c(0, 1),
                          limits = c(0, 1), 
                          na.value = "#b6b7b9")+
    scale_y_continuous(limits=c(0, 2.5), expand = c(0,0))+
    theme_cowplot()+theme(legend.position= "none")+
    labs(x="", y="Copy number")+
    ggtitle("chr6")
  
  gg <- plot_grid(g1, g2, g3, nrow=1)
  
  pdf(paste0(out_dir,'/SNU_18_C_chr4to6.pdf'),
      width=9, height=4)
  print(gg)
  dev.off()
  
}
