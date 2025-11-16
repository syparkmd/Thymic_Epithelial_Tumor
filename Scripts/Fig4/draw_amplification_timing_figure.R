library(tidyverse)

meta_dt <- read_tsv('thymoma_meta_table.250826.tsv') %>% as.data.frame

# meta_dt$mut_rate <- meta_dt$n_pointmt/75  # panel size of SureSelectXT HumanAllExon V5+UTRs, 75MB (75000 kb), alignable genome 2800mbps
# already calculted in meta_dt$TMB
# meta_dt$mut_rate <- meta_dt$n_pointmt/meta_dt$bait_size
meta_dt$mut_rate <- meta_dt$TMB
meta_dt$diploid_proportion <- meta_dt$unchanged_prop
meta_dt$mut_rate
meta_dt$n_pointmt/meta_dt$bait_size

#color setting
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'white')
stage_pal=c("#4352F7", "#DE75E6", "#F52075", "#780000", "#800000", "#4CE8A5")
names(stage_pal) <- c('I','II','III','IVa','IVb', '0')

#age - tmb relationship
tmp_df <- meta_dt %>% filter(GTF2I_status2 %in% c('w','m') & final_cellularity >= 0.4 & diploid_proportion >= 0.80)
nrow(tmp_df) #33
model <- lm(mut_rate ~ age_at_diagnosis, data = tmp_df)
summary(model)

lm(formula = mut_rate ~ age_at_diagnosis, data = tmp_df)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.26373 -0.05606 -0.00992  0.10825  0.33960 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)      0.215193   0.178038   1.209   0.2359  
#age_at_diagnosis 0.006123   0.002762   2.216   0.0341 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.159 on 31 degrees of freedom
#Multiple R-squared:  0.1368,	Adjusted R-squared:  0.1089 
#F-statistic: 4.912 on 1 and 31 DF,  p-value: 0.03413

tmp_df %>%
  with(plot(age_at_diagnosis,mut_rate,ylim=c(0,max(mut_rate)),xlim=c(0,90),
            col=histo_pal[histologic_type]))


# plot --------------------------------------------------------------------
par(mar=c(4.1,4.1,2.1,1.1))
ymax=max(meta_dt$mut_rate[meta_dt$GTF2I_status2 %in% c("w","m")])
plot(-100,-100,ylim=c(0,ymax),xlim=c(0,90),ylab = "Number of point mutation/Gbps",xlab="",yaxt='n', bty="L")
mtext("Age at diagnosis",side=1,line=2)
axis(2,at=0:14/10,labels=0:14*100,las=2)
CI <- predict(model, newdata=tmp_df, interval="confidence",level=0.9) %>% as.data.frame
abline(model$coefficients[1],model$coefficients[2],lwd=1.2,col="grey50",lty=2) # pink
polygon(c(sort(tmp_df$age_at_diagnosis),rev(sort(tmp_df$age_at_diagnosis))),
        c(CI$lwr[order(tmp_df$age_at_diagnosis)],rev(CI$upr[order(tmp_df$age_at_diagnosis)])),
        col="#00000020",lty = 0)
meta_dt[meta_dt$GTF2I_status2 %in% c("w","m")&meta_dt$final_cellularity>=0.5&meta_dt$diploid_proportion>=0.8,] %>%
  with(points(age_at_diagnosis,mut_rate))
ymax=max(meta_dt$mut_rate[meta_dt$GTF2I_status2 %in% c("w","m")])
mtext(side = 1,line = -3.7,at = 0,adj = 0,text = "Purity ≥ 0.4 & Diploid % ≥ 0.80 (n=33)",font=1)
mtext(side = 1,line = -2.3,at = 0,adj = 0,text = bquote(R^2~"="~.(round(summary(model)$r.squared,3))~","~italic(p)~"="~.(round(anova(model)$'Pr(>F)'[1],3))))
mtext(side = 1,line = -1.4,at = 0,adj = 0,text = paste0("slope = ",round(1000*model$coefficients["age_at_diagnosis"],2),"/Gbps/year"))


# timing -----------------------------------------------------------------------------------------
snu19proba <- read_tsv('SNU_19_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba',col_types=cols(`#CHROM`="c"))
snu19seg <-   read_tsv('SNU_19_C_wgs_segments.txt.edit.broadAMP',col_types=cols(`#CHROM`="c"))
snu26proba <- read_tsv('SNU_26_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba',col_types=cols(`#CHROM`="c"))
snu26seg <-   read_tsv('SNU_26_C_wgs.026_31.XY.sequenza3.withBP2_segments.txt.edit.broadAMP',col_types=cols(`#CHROM`="c"))

snu14proba <- read_tsv("SNU_14_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba",col_types=cols(`#CHROM`="c"))
snu14seg <-   read_tsv("SNU_14_C_wgs_seqz_re_Sol1_c0.46_p2.2_gF_segments.txt.edit.broadAMP",col_types=cols(`#CHROM`="c"))
snu17proba <- read_tsv("SNU_17_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba",col_types=cols(`#CHROM`="c"))
snu17seg <-   read_tsv("SNU_17_C_wgs_seqz_re_Sol1_c0.23_p2.1_gF_segments.txt.edit.broadAMP",col_types=cols(`#CHROM`="c"))
snu18proba <- read_tsv("SNU_18_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba",col_types=cols(`#CHROM`="c"))
snu18seg <-   read_tsv("SNU_18_C_wgs_seqz_re_Sol1_c0.83_p1.9_gF_segments.txt.edit.broadAMP",col_types=cols(`#CHROM`="c"))
tm_dt1 <- read_tsv('SNU_26_C.026_31.broad_earlysnv')
tm_dt2 <- read_tsv('SNU_19_C.079_21.broad_earlysnv')
tm_dt3 <- read_tsv('SNU_14_C.046_22.broad_earlysnv')
tm_dt4 <- read_tsv('SNU_17_C.023_21.broad_earlysnv')
tm_dt5 <- read_tsv('SNU_18_C.083_19.broad_earlysnv')

tm_dt1$early_SNV/tm_dt1$dist*1E6
tm_dt1$case = "SNU_26_C"
tm_dt2$case = "SNU_19_C"
tm_dt3$case = "SNU_14_C"
tm_dt4$case = "SNU_17_C"
tm_dt5$case = "SNU_18_C"

tm_dt <- bind_rows(tm_dt1,tm_dt2,tm_dt3,tm_dt4,tm_dt5)
tm_dt$id = tm_dt$case


snu19.early.mut.info <- read_rds("~/00_Project/01_thymoma/kjyi/data/snu19.early.mut.info2.Rds")
snu26.early.mut.info <- read_rds("~/00_Project/01_thymoma/kjyi/data/snu26.early.mut.info2.Rds")
snu14.early.mut.info <- read_rds("~/00_Project/01_thymoma/kjyi/data/snu14.early.mut.info2.Rds")
snu17.early.mut.info <- read_rds("~/00_Project/01_thymoma/kjyi/data/snu17.early.mut.info2.Rds")
snu18.early.mut.info <- read_rds("~/00_Project/01_thymoma/kjyi/data/snu18.early.mut.info2.Rds")

snu19.early.mut.info$case = "SNU_19_C"
snu26.early.mut.info$case = "SNU_26_C"
snu14.early.mut.info$case = "SNU_14_C"
snu17.early.mut.info$case = "SNU_17_C"
snu18.early.mut.info$case = "SNU_18_C"

tm_dt_new <- rbind(snu19.early.mut.info,
                   snu26.early.mut.info,
                   snu14.early.mut.info,
                   snu17.early.mut.info,
                   snu18.early.mut.info)
tm_dt_new

for(col in colnames(tm_dt_new)){tm_dt_new[[col]] <- unlist(tm_dt_new[[col]])}
tm_dt_new <- tm_dt_new %>% 
  mutate(mr = est.early.mut.count/(end-start)*1E6,
         mr.upper = est.early.mut.count.CIU/(end-start)*1E6,
         mr.lower = est.early.mut.count.CIL/(end-start)*1E6)
tm_dt_new


# tm_dt$`early_SNVrate/Mb`
# tm_dt$early_SNVrate_range_low
# tm_dt$early_SNVrate_range_high

age_at_diagnosis2 = c("SNU_19_C"=45,
                      "SNU_26_C"=19,
                      "SNU_14_C"=48,
                      "SNU_17_C"=62,
                      "SNU_18_C"=40)
genome_size_in_mb=2.7*10^9/10^6
# total_mut_rate = c("SNU_19_C_wgs"=1138,"SNU_26_C_wgs"=364) / genome_size_in_mb

total_mut_count = c("SNU_19_C"=nrow(snu19proba),
                    "SNU_26_C"=nrow(snu26proba),
                    "SNU_14_C"=nrow(snu14proba),
                    "SNU_17_C"=nrow(snu17proba),
                    "SNU_18_C"=nrow(snu18proba)) 
# nrow(snu19proba)==843
# nrow(snu26proba)==225

tm_dt%>% nrow() # 37 rows
tm_dt_new %>% nrow() # 40 rows


# mutation_rate = 0.006122603  # ----updated, 6.122603  / Gbps/year
mutation_rate = model$coefficients[2]


#A function to add arrows on the chart
par(mar=c(10,11,4,4))
tm_dt <- arrange(tm_dt,desc(case),desc(`early_SNVrate/Mb`))
tm_dt_new <- arrange(tm_dt_new,desc(case),desc(mr))
tm_dt_new = tm_dt_new[tm_dt_new$num.mut > 0,]
tm_dt$total_clonal_SNV_rate
tm_dt$`early_SNVrate/Mb`
total_mut_count/genome_size_in_mb
tm_barplot <- barplot(tm_dt_new$mr,
                      names.arg = tm_dt_new$chr_arm,
                      horiz=T,
                      las=2,
                      # xlim=c(-0.005,mutation_rate*(100)),
                      xlim=c(-0.005,0.653),
                      xaxt='n')



axis(side = 3,at=0:5/10,labels = 0:5/10*1000)
axis(side = 1,at=0:5/10,labels = 0:5/10*1000)
# axis(side=1)
axis(side = 1,line=3,at = mutation_rate*(0:8*10),labels = (0:8*10))

mtext(text = paste0("Estimated age of amplification\nbased on mutation rate ",round(mutation_rate*1000,2),"/Gbps/year"),line = 6,side = 1)
mtext("Number of amplified (early) base substitutions / Gbps",line = 2.5)
axis(side = 2,at = tm_barplot,labels = rep("",length(tm_barplot)))
arrows(x0 = pmax(0,tm_dt_new$mr.lower),
       x1 = pmax(0,tm_dt_new$mr.upper),
       y0 = tm_barplot,
       angle=90, code=3, length=0.025)
# tm_dt$id = c(rep("SNU_26_C",18),rep("SNU_19_C",3))
for(i in unique(tm_dt_new$case)){
  axis(side = 2,at = tm_barplot[c(range(which(tm_dt_new$case==i)))],labels = c("",""),line = 5.1,tck=0.02)
  mtext(text = str_replace(i,"_wgs",""),side = 2,outer = F,line = 5.6,
        at = median(tm_barplot[c(range(which(tm_dt_new$case==i)))]),las=2,
        col="black")
  segments(x0 = age_at_diagnosis2[i]*mutation_rate,
           y0 = tm_barplot[min(which(tm_dt_new$case==i))]-diff(tm_barplot)[1]/2,
           y1 = tm_barplot[max(which(tm_dt_new$case==i))]+diff(tm_barplot)[1]/2,
           lty=2)
  segments(x0 = total_mut_count[i]/genome_size_in_mb,
           y0 = tm_barplot[min(which(tm_dt_new$case==i))]-diff(tm_barplot)[1]/2,
           y1 = tm_barplot[max(which(tm_dt_new$case==i))]+diff(tm_barplot)[1]/2,
           lty=1,col="red")
}
# legend("topright",legend = "age at diagnosis",lty=2,bty="n")
legend("bottomright",legend = c("Age at diagnosis","Total subsitutions / Gbps"),lty=c(2,1),col=c("black","red"),bty="n")


