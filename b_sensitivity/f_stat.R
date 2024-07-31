## F statistics
# Univariate analysis -- mean F-statistic
# Multivariate analysis -- conditional F-statistic

library(readr);library(TwoSampleMR);library(data.table);library(readxl);library(remotes);library(stringr)
ao=available_outcomes()
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)
pheno_cov=read.delim("pcov_nmr_dat_all.txt")

## UKBB
nmr_metabolites=read_xlsx("nmr_metabolites_20210610.xlsx")[,-1]
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"&nmr_metabolites$Group%in%c("Glycolysis related metabolites", "Amino acids", "Ketone bodies")),]
colnames(ao)[2]="Biomarker name"
ao1=ao[which(ao$author=="Borges CM"&ao$`Biomarker name`%in%nmr_metabolites$`Biomarker name`),]
nmr_metabolites=merge(nmr_metabolites, ao1, by="Biomarker name")
pheno_cov=pheno_cov[which(rownames(pheno_cov)%in%nmr_metabolites$id), which(colnames(pheno_cov)%in%nmr_metabolites$id)]
exposure_data=extract_instruments(nmr_metabolites$id, p1=5e-8, clump=T, r2=0.01)
outcome_dat=read_outcome_data(snps=exposure_data$SNP, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
uvmr_harm=harmonise_data(exposure_data, outcome_dat, action=2)
uvmr_harm$pre_fstat=(uvmr_harm$beta.exposure^2)/(uvmr_harm$se.exposure^2)

models=as.data.frame(read_csv("ukbb_best_model_out_default_3subclasses.csv"))
models$gwas_id=NA
for (i in models$`rf combination`)
{
  models$gwas_id[which(models$`rf combination`==i)]=ao$id[which(ao$`Biomarker name`==i&ao$author=="Borges CM")]
}

mean_fstat_ukbb=data.frame(NA,NA, NA)
colnames(mean_fstat_ukbb)=c("model", "f-stat","f-stat, all SNPs")
for (i in models$gwas_id)
{
  j=which(models$gwas_id==i)
  mean_fstat_ukbb[j,]=c(paste(i), mean(uvmr_harm$pre_fstat[which(uvmr_harm$id.exposure==i)]),
                        round(sum(uvmr_harm$pre_fstat[which(uvmr_harm$id.exposure==i)])/sum(uvmr_harm$pre_fstat),2))
  
}


write.csv(mean_fstat_ukbb,"mean_fstat_ukbb_3subclasses.csv" )

## Kett

nmr_metabolites=read_excel("nmr_metabolites_20210610.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"&nmr_metabolites$Group%in%c("Glycolysis related metabolites", "Amino acids", "Ketone bodies")),]
colnames(nmr_metabolites)[1]="GWAS_id"
nmr_metabolites$GWAS_id=paste0("met-d-", nmr_metabolites$GWAS_id)
name_mismatch="3-hydroxybutyrate"
ao=available_outcomes()
ao=ao[which(ao$year==2016&ao$author=="Kettunen"&ao$category=="Metabolites"&ao$population=="European"),]
ao_1=ao[which(ao$trait%in%name_mismatch),]
ao_1=rbind(ao_1,ao[which(ao$trait%in%nmr_metabolites$`Biomarker name`),])
list_of_datasets=as.data.frame(ao_1)
ket_ids=as.character(list_of_datasets$id)
uvmr_instruments=extract_instruments(ket_ids,p1=5e-8,clump=T,r2=0.01)
outcome_dat=read_outcome_data(snps=uvmr_instruments$SNP, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
uvmr_harm=harmonise_data(uvmr_instruments, outcome_dat, action=2)
uvmr_harm$pre_fstat=(uvmr_harm$beta.exposure^2)/(uvmr_harm$se.exposure^2)
mean_fstat_kett=data.frame(NA,NA)
colnames(mean_fstat_kett)=c("model", "f-stat")

ao_1$trait[which(ao_1$trait==name_mismatch)]="3-Hydroxybutyrate"
models2=as.data.frame(read_csv("kett_best_model_out_default_3subclasses.csv")[,-1])

models2$gwas_id=NA
for (i in models2$`rf combination`)
{
  if (i%in%ao_1$trait){
    models2$gwas_id[which(models2$`rf combination`==i)]=ao_1$id[which(ao_1$trait==i)]
  }
  else
    models2$gwas_id[which(models2$`rf combination`==i)]=NA
}

for (i in na.omit(models2)$gwas_id)
{
  j=which(na.omit(models2)$gwas_id==i)
  mean_fstat_kett[j,]=c(paste(na.omit(models2)$`rf combination`)[j], mean(uvmr_harm$pre_fstat[which(uvmr_harm$id.exposure==i)]))
}

write.csv(mean_fstat_kett,"mean_fstat_kett_3subclasses.csv" )

