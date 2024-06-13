## Heatmaps
library(readr);library(gplots);library(readxl);library(TwoSampleMR);library(data.table)

nmr_metabolites=read_excel("nmr_metabolites_20210610.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]

# I. UKBB
mvmr_harm=read_csv("mvmr_harm.csv")[,-1]
setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.matrix(test)
test=as.data.frame(test)[,grep("beta.exposure", colnames(test))]
ao=available_outcomes()
ao_1=ao[which(ao$year==2020&ao$author=="Borges CM"),]
ao_1=ao_1[which(ao_1$trait%in%nmr_metabolites$`Biomarker name`),]

colnames(test)=sub("beta.exposure.","", colnames(test))
test=apply(test, 2, function(x) as.numeric(as.character(x)))
sapply(test, class)
for (j in 1:nrow(ao_1))
{
  colnames(test)=gsub(ao_1$id[j], ao_1$trait[j], colnames(test))
  colnames(test)=gsub(ao$id[j], ao_1$trait[j], colnames(test))
}
cormatrix=cor(test)
pdf('heatmap_ukbb_snp_redo.pdf', 16, 16)
heatmap.2(cormatrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram='none',
          Rowv=TRUE, Colv=TRUE, keysize=1, margins=c(18.5,20))
graphics.off()

# II. Kett
mvmr_harm=read_csv("mvmr_harm.csv")[,-1]
setDT(mvmr_harm)
test=reshape(mvmr_harm, timevar="id.exposure", idvar=c("SNP"), direction="wide")
test=as.data.frame(test)[,grep("beta.exposure", colnames(test))]
ao=available_outcomes()
ao_1=ao[which(ao$year==2020&ao$author=="Borges CM"),]
colnames(test)=sub("beta.exposure.","", colnames(test))

test=apply(test, 2, function(x) as.numeric(as.character(x)))
sapply(test, class)
for (j in 1:nrow(ao_1))
{
  colnames(test)=gsub(ao_1$id[j], ao_1$trait[j], colnames(test))
  colnames(test)=gsub(ao$id[j], ao_1$trait[j], colnames(test))
}
cormatrix=cor(test)
pdf('heatmap_kett_snp.pdf', 12, 12)
heatmap.2(cormatrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram='none',
          Rowv=TRUE, Colv=TRUE, keysize=1, margins=c(18.5,20))
graphics.off()

# III. Pheno-cov

pheno_cov=read.delim("pcov_nmr_dat_all.txt")
nmr_metabolites_UKBB=read_xlsx("nmr_metabolites_20210610.xlsx")[,-1]
nmr_metabolites_UKBB=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$Include=="yes"),]
colnames(ao)[2]="Biomarker name"
ao1=ao[which(ao$author=="Borges CM"&ao$`Biomarker name`%in%nmr_metabolites_UKBB$`Biomarker name`),]
nmr_metabolites_UKBB=merge(nmr_metabolites_UKBB, ao1, by="Biomarker name")
nmr_metabolites_UKBB$id=sub("met-d-", "", nmr_metabolites_UKBB$id)
pheno_cov=pheno_cov[which(rownames(pheno_cov)%in%nmr_metabolites_UKBB$id), which(colnames(pheno_cov)%in%nmr_metabolites_UKBB$id)]
colnames(pheno_cov)==rownames(pheno_cov)
colnames(pheno_cov)=rownames(pheno_cov)=paste0("met-d-", colnames(pheno_cov))

colnames(pheno_cov)=ao1$`Biomarker name`[match(colnames(pheno_cov), ao1$id)]
rownames(pheno_cov)=ao1$`Biomarker name`[match(rownames(pheno_cov), ao1$id)]
colnames(pheno_cov)=gsub("Total lipids in chylomicrons and extremely large VLDL","Total lipids in QM & extremely large VLDL", colnames(pheno_cov))
rownames(pheno_cov)=gsub("Total lipids in chylomicrons and extremely large VLDL","Total lipids in QM & extremely large VLDL", rownames(pheno_cov))

cormatrix=cor(pheno_cov)

pdf('heatmap_phenocov.pdf', 12, 12)
heatmap.2(cormatrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none", dendrogram='none',
          Rowv=TRUE, Colv=TRUE, keysize=1, margins=c(14,16))
graphics.off()
