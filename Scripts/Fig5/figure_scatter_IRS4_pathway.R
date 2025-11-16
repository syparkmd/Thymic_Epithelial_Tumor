#copied from kjyi script

library(tidyverse)
library(scales)
library(circlize)
library(fgsea)
require(GSVA)
# pals
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

# load exp data
meta_dt <- read_tsv('thymoma_meta_table.250826.tsv')
bm_dt <- read_tsv('biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename.txt') %>%
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()
exp_dt <- read_tsv('thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')
exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
rm(bm_dt,exp_dt_pcg,exp_dt)
# load gene set
h_list <- GSEABase::getGmt('h.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c5_list <- GSEABase::getGmt('c5.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c2_list <- GSEABase::getGmt('c2.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
all_list <- all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG","REACTOME")]

all_list["REACTOME_ENERGY_DEPENDENT_REGULATION_OF_MTOR_BY_LKB1_AMPK"]

# run gsva with log10(x+0.01)
l10_exp_dt_mat <- l10_exp_dt %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix
l10_exp_dt_mat <- l10_exp_dt_mat[,meta_dt$id]
all_gs_res <- l10_exp_dt_mat %>% gsva(all_list)
# all_gs_res <- all_gs_res[,meta_dt$id]

l10_exp_dt_mat[1:3,1:3]
all_gs_res[1:3,1:3]
colnames(l10_exp_dt_mat) == colnames(all_gs_res)
lm_p1 = foreach(i = 1:nrow(all_gs_res),.combine=rbind) %do% {
  res = cor.test(y = l10_exp_dt_mat["IRS4",meta_dt$GTF2I_status2=="w"],
                 x = all_gs_res[i,meta_dt$GTF2I_status2=="w"], method = "pearson")
  res2 = cor.test(y = l10_exp_dt_mat["IRS4",meta_dt$GTF2I_status2=="w"],
                  x = all_gs_res[i,meta_dt$GTF2I_status2=="w"], method = "spearman")
  c(res$estimate,res$p.value,res2$estimate,res2$p.value)
}
rownames(lm_p1) = rownames(all_gs_res)
colnames(lm_p1) = c("Pearson coeficient", "Pearson p-value", "Spearman coeficient","Spearman p-value")
View(lm_p1)

lm_p2 = foreach(i = 1:nrow(all_gs_res),.combine=rbind) %do% {
  summary(lm(all_gs_res[i,meta_dt$GTF2I_status2=="w"] ~ 
               meta_dt$final_cellularity[meta_dt$GTF2I_status2=="w"] + 
               l10_exp_dt_mat["IRS4",meta_dt$GTF2I_status2=="w"]))$coefficients -> x
  c(x[2:3,4],x[2:3,1])
}
rownames(lm_p2) = rownames(all_gs_res)
colnames(lm_p2) = c("lm purity p-value", "lm IRS4 p-value", "lm purity coeficient", "lm IRS4 coeficient")
lm_bind = cbind(lm_p1,lm_p2)
if(F) View(lm_bind)
library(rJava)
library(xlsx)
if(F) xlsx::write.xlsx(lm_bind,file = "tables/S.Table3.xlsx")



# plot -------------------------------------------------------------------------
WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w"]
MT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="m"]
# all_gs_res.bak <- all_gs_res
# all_gs_res <- all_gs_res.bak
all_gs_res_WT <- all_gs_res[,WT_samplenames]
all_gs_res_MT <- all_gs_res[,MT_samplenames]
geneexpresionmatrix_WT <- l10_exp_dt_mat[,WT_samplenames]
geneexpresionmatrix_MT <- l10_exp_dt_mat[,MT_samplenames]


tmp_df_WT <- data.frame(IRS4=pmax(0,geneexpresionmatrix_WT["IRS4",]),
                        IRS4_orig=geneexpresionmatrix_WT["IRS4",],
                        STAR=pmax(0,geneexpresionmatrix_WT["STAR",]),
                        STAR_orig=geneexpresionmatrix_WT["STAR",],
                        IGF1R=pmax(0,geneexpresionmatrix_WT["IGF1R",]),
                        IGF1R_orig=geneexpresionmatrix_WT["IGF1R",],
                        IRS2=pmax(0,geneexpresionmatrix_WT["IRS2",]),
                        IRS2_orig=geneexpresionmatrix_WT["IRS2",],
                        PIK3R1=pmax(0,geneexpresionmatrix_WT["PIK3R1",]),
                        PIK3R1_orig=geneexpresionmatrix_WT["PIK3R1",],
                        PIK3CA=pmax(0,geneexpresionmatrix_WT["PIK3CA",]),
                        IGFBPL1=pmax(0,geneexpresionmatrix_WT["IGFBPL1",]),
                        IRS1=pmax(0,geneexpresionmatrix_WT["IRS1",]),
                        IRS1_orig=geneexpresionmatrix_WT["IRS1",],
                        IGFsig=all_gs_res_WT["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                        # IGFsig=all_gs_res["GO_INSULIN_RECEPTOR_SUBSTRATE_BINDING",],
                        IRSbinding=all_gs_res_WT["GO_INSULIN_RECEPTOR_SUBSTRATE_BINDING",],
                        lipid=all_gs_res_WT["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                        PI3Kcascade=all_gs_res_WT["REACTOME_PI3K_CASCADE",],
                        mito=all_gs_res_WT["GO_MITOCHONDRION_LOCALIZATION",],
                        ventri=all_gs_res_WT["GO_VENTRICULAR_SYSTEM_DEVELOPMENT",],
                        pib=all_gs_res_WT["GO_PHOSPHATIDYLINOSITOL_BINDING",],
                        pipb=all_gs_res_WT["GO_PHOSPHATIDYLINOSITOL_PHOSPHATE_BINDING",],
                        # histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],
                        purity=meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)])

tmp_df_MT <- data.frame(IRS4=pmax(0,geneexpresionmatrix_MT["IRS4",]),
                        IRS4_orig=geneexpresionmatrix_MT["IRS4",],
                        STAR=pmax(0,geneexpresionmatrix_MT["STAR",]),
                        STAR_orig=geneexpresionmatrix_MT["STAR",],
                        IGF1R=pmax(0,geneexpresionmatrix_MT["IGF1R",]),
                        IGF1R_orig=geneexpresionmatrix_MT["IGF1R",],
                        IRS2=pmax(0,geneexpresionmatrix_MT["IRS2",]),
                        IRS2_orig=geneexpresionmatrix_MT["IRS2",],
                        PIK3R1=pmax(0,geneexpresionmatrix_MT["PIK3R1",]),
                        PIK3R1_orig=geneexpresionmatrix_MT["PIK3R1",],
                        PIK3CA=pmax(0,geneexpresionmatrix_MT["PIK3CA",]),
                        IGFBPL1=pmax(0,geneexpresionmatrix_MT["IGFBPL1",]),
                        IRS1=pmax(0,geneexpresionmatrix_MT["IRS1",]),
                        IRS1_orig=geneexpresionmatrix_MT["IRS1",],
                        IGFsig=all_gs_res_MT["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                        IRSbinding=all_gs_res_MT["GO_INSULIN_RECEPTOR_SUBSTRATE_BINDING",],
                        lipid=all_gs_res_MT["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                        PI3Kcascade=all_gs_res_MT["REACTOME_PI3K_CASCADE",],
                        mito=all_gs_res_MT["GO_MITOCHONDRION_LOCALIZATION",],
                        ventri=all_gs_res_MT["GO_VENTRICULAR_SYSTEM_DEVELOPMENT",],
                        pib=all_gs_res_MT["GO_PHOSPHATIDYLINOSITOL_BINDING",],
                        pipb=all_gs_res_MT["GO_PHOSPHATIDYLINOSITOL_PHOSPHATE_BINDING",],
                        # histol=histo_pal[meta_dt$histologic_type[match(MT_samplenames,meta_dt$id)]],
                        purity=meta_dt$final_cellularity[match(MT_samplenames,meta_dt$id)])

# summary(lm(IGFsig ~ IRS4_orig+purity, data=tmp_df))
# summary(lm(IGFsig ~ purity, data=tmp_df))
# 
# summary(lm(PIK3R1_orig ~ IRS4_orig+purity, data=tmp_df))
# summary(lm(IRS4_orig ~ PIK3R1_orig+purity, data=tmp_df))
# summary(lm(IRS4_orig ~ PIK3R1_orig+purity, data=tmp_df))



# model1 <- lm(IGFsig ~ IRS4_orig, data=tmp_df)
# model2 <- lm(lipid ~ IRS4_orig, data=tmp_df)
# model3 <- lm(PI3Kcascade ~ IRS4_orig, data=tmp_df)
# model4 <- lm(mito ~ IRS4_orig, data=tmp_df)

model1 <- lm(IGFsig ~ IRS4_orig, data=tmp_df_WT)
model2 <- lm(PI3Kcascade ~ IRS4_orig, data=tmp_df_WT)
model3 <- lm(IGF1R_orig ~ IRS4_orig, data=tmp_df_WT)
model4 <- lm(PIK3R1_orig ~ IRS4_orig, data=tmp_df_WT)


# cor(tmp_df$IRS4_orig,tmp_df$IGFsig)
# cor(tmp_df$IRS4_orig,tmp_df$lipid)
# cor(tmp_df$IRS4_orig,tmp_df$PI3Kcascade)
# cor(tmp_df$IRS4_orig,tmp_df$mito)
CI1 <- predict(model1, newdata=tmp_df_WT, interval="prediction",level=0.9) %>% as.data.frame
CI2 <- predict(model2, newdata=tmp_df_WT, interval="prediction",level=0.9) %>% as.data.frame
CI3 <- predict(model3, newdata=tmp_df_WT, interval="prediction",level=0.9) %>% as.data.frame
CI4 <- predict(model4, newdata=tmp_df_WT, interval="prediction",level=0.9) %>% as.data.frame


# svglite::htmlSVG(width=1150/254,height=1150/254,pointsize = 12*0.7, {
cairo_pdf("figures/IRS4_pathway_cor.2.pdf",width=1150/254,height=1050/254,pointsize = 12*0.7)

par(mar=c(4.5,4.5,0,0),oma=c(0,0,1,1),mfrow=c(2,2))



# 1
plot(1000,xlim=c(0,max(tmp_df_WT$IRS4)),ylim=c(min(tmp_df_WT$IGFsig),max(tmp_df_WT$IGFsig)),
     ylab="",las=2,xaxt="n",xlab="")
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("IGF receptor signaling score",2,2.45)
lines(x=tmp_df_WT$IRS4_orig, y=CI1$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df_WT$IRS4_orig),rev(sort(tmp_df_WT$IRS4_orig))),c(CI1$lwr[order(tmp_df_WT$IRS4_orig)],rev(CI1$upr[order(tmp_df_WT$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df_WT$IRS4,tmp_df_WT$IGFsig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df_WT$IGFsig),
     bquote(r~"="~.(round(summary(model1)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df_WT$IGFsig),
     bquote(
       italic(p)~"="~.(format(anova(model1)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext(expression(italic("IRS4")*" expression (TPM)"),1,2.0)
# legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)


# 2
plot(1000,xlim=c(0,max(tmp_df_WT$IRS4)),ylim=c(min(tmp_df_WT$PI3Kcascade),max(tmp_df_WT$PI3Kcascade)),
     xlab="",ylab="",xaxt="n",las=2)
mtext("PI3K cascade score",2,2.45)
lines(x=tmp_df_WT$IRS4_orig, y=CI2$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df_WT$IRS4_orig),rev(sort(tmp_df_WT$IRS4_orig))),c(CI2$lwr[order(tmp_df_WT$IRS4_orig)],rev(CI2$upr[order(tmp_df_WT$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df_WT$IRS4,tmp_df_WT$PI3Kcascade,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
text(0,max(tmp_df_WT$PI3Kcascade),
     bquote(r~"="~.(round(summary(model2)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df_WT$PI3Kcascade),
     bquote(
       italic(p)~"="~.(format(anova(model2)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext(expression(italic("IRS4")*" expression (TPM)"),1,2.0)
# mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df_WT)),")")))

# 3
plot(1000,xlim=c(0,max(tmp_df_WT$IRS4)),ylim=c(min(tmp_df_WT$IGF1R_orig),max(tmp_df_WT$IGF1R_orig)),
     xlab="",ylab="",xaxt="n",las=2,yaxt='n')
axis(2,at = c(0,log10(5),log10(10),log10(25)),labels = c(0,5,10,25))
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext(expression(italic("IGF1R")*" expression (TPM)"),2,2.45)
lines(x=tmp_df_WT$IRS4_orig, y=CI3$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df_WT$IRS4_orig),rev(sort(tmp_df_WT$IRS4_orig))),c(CI3$lwr[order(tmp_df_WT$IRS4_orig)],rev(CI3$upr[order(tmp_df_WT$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df_WT$IRS4,tmp_df_WT$IGF1R_orig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df_WT$IGF1R_orig),
     bquote(r~"="~.(round(summary(model3)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df_WT$IGF1R_orig),
     bquote(
       italic(p)~"="~.(format(anova(model3)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
# mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df_WT)),")")),outer=T)
# axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext(expression(italic("IRS4")*" expression (TPM)"),1,2.0)


#4

plot(1000,xlim=c(0,max(tmp_df_WT$IRS4)),ylim=c(min(tmp_df_WT$PIK3R1_orig),max(tmp_df_WT$PIK3R1_orig)),
     ylab="",las=2,xaxt="n",yaxt="n",xlab="")
axis(2,at = c(0,1,2,log10(300)),labels = c(0,10,100,300))
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext(expression(italic("PIK3R1")*" expression (TPM)"),2,2.45)
lines(x=tmp_df_WT$IRS4_orig, y=CI4$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df_WT$IRS4_orig),rev(sort(tmp_df_WT$IRS4_orig))),c(CI4$lwr[order(tmp_df_WT$IRS4_orig)],rev(CI4$upr[order(tmp_df_WT$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df_WT$IRS4,tmp_df_WT$PIK3R1_orig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df_WT$PIK3R1_orig),
     bquote(r~"="~.(round(summary(model4)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df_WT$PIK3R1_orig),
     bquote(
       italic(p)~"="~.(format(anova(model4)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext(expression(italic("IRS4")*" expression (TPM)"),1,2.0)
legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)
# })
dev.off()



# linear model comparison for wt, mt, irs1 2 4

fitplot <- function(fit, main="Estimated coefficients for 1-puirity",
                    labels=c("(Intercept)","NK cell score","Neutrophil score",
                             "Dendritic cell score","Cytotoxic score","Exhausted score",
                             "Macrophage score","B cell score","Thymopoiesis score")){

  y = coef(fit)[-1]
  x = rev(seq_along(y))
  ci = confint(fit)[-1,]
  xlim = range(x) + c(-0.5,0.2)
  ylim = range(ci)
  ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
  # ylab = bquote(hat(beta))
  ylab=""
  xlab = "coefficient"
  ss <-summary(fit)
  sss<- ifelse(ss$coefficients[,4] < 0.001, "***", ifelse(ss$coefficients[,4] < 0.05,"*",""))
  sss<-sss[-1]
  par(mar=c(2,8,3,4))
  # plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="",         bty="n", main = main)
  plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n", main = main,cex.main=1.05)
  axis(2, at=x, labels=labels, tick=FALSE,las=1)
  axis(4, at=x, labels=sss, tick=FALSE,las=2)
  abline(v=0, lty=3)
  arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
  summary(fit)
}

par(mfrow=c(4,4))
fitplot(fit = lm(IGFsig ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_WT))),
        main=expression(bold("Estimated coefficients for IGF signaling score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(IGFsig ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_MT))),
        main=expression(bold("Estimated coefficients for IGF signaling score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))

fitplot(fit = lm(IRSbinding ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_WT))),
        main=expression(bold("Estimated coefficients for IRS binding score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(IRSbinding ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_MT))),
        main=expression(bold("Estimated coefficients for IRS binding score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))

fitplot(fit = lm(PI3Kcascade ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig  + purity, data=as.data.frame(scale(tmp_df_WT))),
        main=expression(bold("Estimated coefficients for PI3K cascade score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(PI3Kcascade ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig  + purity, data=as.data.frame(scale(tmp_df_MT))),
        main=expression(bold("Estimated coefficients for PI3K cascade score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))



fitplot(fit = lm(lipid ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for intracellular lipid transport score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(lipid ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for intracellular lipid transport score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(mito ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for mitochondira distribution score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(mito ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for mitochondira distribution score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))

fitplot(fit = lm(ventri ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for vsd score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(ventri ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for vsd score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))


tmp_df_WT$IGFsig

all_list[""]

# candidates
par(mfrow=c(4,2),pty='s')
dim(tmp_df_WT)
fitplot(fit = lm(IRSbinding ~ IRS4_orig + IRS1_orig + IRS2_orig + STAR_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for IRS binding score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","STAR","Tumor cell fraction"))
fitplot(fit = lm(IRSbinding ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for IRS binding score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(PI3Kcascade ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for PI3K cascade score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(PI3Kcascade ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for PI3K cascade score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(lipid ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for intracellular lipid transport score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(lipid ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for intracellular lipid transport score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(mito ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for mitochondira distribution score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(mito ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for mitochondira distribution score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(pib ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_WT[,-18]))),
        main=expression(bold("Estimated coefficients for pib score in "*bolditalic("GTF2I")*" wildtypes")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))
fitplot(fit = lm(pib ~ IRS4_orig + IRS1_orig + IRS2_orig + purity, data=as.data.frame(scale(tmp_df_MT[,-18]))),
        main=expression(bold("Estimated coefficients for pib score in "*bolditalic("GTF2I")*" mutants")),
        labels=c("IRS4","IRS1","IRS2","Tumor cell fraction"))


# supplementary figure 3b
cairo_pdf("figures/IRS124_lm_IRS_PIP.pdf",width=190/25.4,height=80/25.4,pointsize = 12*0.7)
1200/72*25.4
666
par(mfrow=c(2,2),pty="m")
"GO_INSULIN_RECEPTOR_SUBSTRATE_BINDING"
"GO_PHOSPHATIDYLINOSITOL_BINDING"
fitplot(fit = lm(IRSbinding ~ IRS1_orig + IRS2_orig + IRS4_orig + purity, data=as.data.frame(scale(tmp_df_WT))),
        main=expression(bold("Estimated coefficients for IRS binding score in "*bolditalic("GTF2I")*"-wildtypes")),
        labels=c(expression(italic("IRS1")),
                 expression(italic("IRS2")),
                 expression(italic("IRS4")),"Tumor cell fraction"))
fitplot(fit = lm(IRSbinding ~ IRS1_orig + IRS2_orig + IRS4_orig + purity, data=as.data.frame(scale(tmp_df_MT))),
        main=expression(bold("Estimated coefficients for IRS binding score in "*bolditalic("GTF2I")*"-mutants")),
        labels=c(expression(italic("IRS1")),
                 expression(italic("IRS2")),
                 expression(italic("IRS4")),"Tumor cell fraction"))
fitplot(fit = lm(pipb ~ IRS1_orig + IRS2_orig + IRS4_orig + purity, data=as.data.frame(scale(tmp_df_WT))),
        main=expression(bold("Estimated coefficients for PIP binding score in "*bolditalic("GTF2I")*"-wildtypes")),
        labels=c(expression(italic("IRS1")),
                 expression(italic("IRS2")),
                 expression(italic("IRS4")),"Tumor cell fraction"))
fitplot(fit = lm(pipb ~ IRS1_orig + IRS2_orig + IRS4_orig + purity, data=as.data.frame(scale(tmp_df_MT))),
        main=expression(bold("Estimated coefficients for PIP binding score in "*bolditalic("GTF2I")*"-mutants")),
        labels=c(expression(italic("IRS1")),
                 expression(italic("IRS2")),
                 expression(italic("IRS4")),"Tumor cell fraction"))
dev.off()
