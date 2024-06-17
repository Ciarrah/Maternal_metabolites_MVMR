library(dplyr);library(readr);library(TwoSampleMR);library(stringr);library(readxl)
## UKBB
setwd("ukbb_snps")
nmr_metabolites=read_excel("nmr_metabolites_20210610.xlsx")
ao=available_outcomes()
trait_linkage=ao[grep("met-d", ao$id),]
nmr_metabolites$Name=paste0("met-d-", nmr_metabolites$Name)
names(nmr_metabolites)[which(names(nmr_metabolites)=="Name")]="GWAS_id"
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$`Biomarker name`%in%c("Alanine","3-Hydroxybutyrate","Glutamine","Glucose","Isoleucine","Pyruvate")),]
exposure_data=extract_instruments(nmr_metabolites$GWAS_id, p1=5e-8, clump=T, r2=0.01)
outcome_dat=read_outcome_data(snps=exposure_data$SNP, filename="/Volumes/MRC-IEU-research/projects/ieu2/p6/106/working/data/Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
uvmr_harm=harmonise_data(exposure_data, outcome_dat, action=2)
ao=ao[grep("met-d", ao$id),]
ao=ao[which(ao$id%in%uvmr_harm$id.exposure),]
estimates_all=mr(uvmr_harm)
for (j in 1:nrow(ao))
{
  estimates_all$id.exposure=gsub(ao$id[j], ao$trait[j], estimates_all$id.exposure)
}

mr_egger=estimates_all[which(estimates_all$method=="MR Egger"),]
weighted_median=estimates_all[which(estimates_all$method=="Weighted median"),]
weighted_mode=estimates_all[which(estimates_all$method=="Weighted mode"),]

mr_egger[,c("l_ci","u_ci")]=c(mr_egger$b-(1.96*mr_egger$se),mr_egger$b+(1.96*mr_egger$se))
weighted_median[,c("l_ci","u_ci")]=c(weighted_median$b-(1.96*weighted_median$se),weighted_median$b+(1.96*weighted_median$se))
weighted_mode[,c("l_ci","u_ci")]=c(weighted_mode$b-(1.96*weighted_mode$se),weighted_mode$b+(1.96*weighted_mode$se))

mr_egger[,c("b","se","u_ci","l_ci","pval")]=round(mr_egger[,c("b","se","u_ci","l_ci","pval")],3)
weighted_median[,c("b","se","u_ci","l_ci","pval")]=round(weighted_median[,c("b","se","u_ci","l_ci","pval")],3)
weighted_mode[,c("b","se","u_ci","l_ci","pval")]=round(weighted_mode[,c("b","se","u_ci","l_ci","pval")],3)

mr_egger=mr_egger[order(mr_egger$id.exposure),]
weighted_median=weighted_median[order(weighted_median$id.exposure),]
weighted_mode=weighted_mode[order(weighted_mode$id.exposure),]

pdf("ukbb_snps/mr_scatter_plots.pdf")
mr_scatter_plot(mr(uvmr_harm),uvmr_harm)
dev.off()
res_loo=mr_leaveoneout(uvmr_harm)
write.csv(res_loo,"ukbb_snps/res_loo.csv")
pdf("ukbb_snps/mr_res_loo.pdf")
mr_leaveoneout_plot(res_loo)
dev.off()

View(mr_egger)
View(weighted_median)
View(weighted_mode)

write.csv(mr_egger,"ukbb_snps/estimates_egger.csv")
write.csv(weighted_median,"ukbb_snps/estimates_median.csv")
write.csv(weighted_mode,"ukbb_snps/estimates_mode.csv")

#####
#### final plots

for (j in 1:nrow(ao))
{
  uvmr_harm$id.exposure=gsub(ao$id[j], ao$trait[j], uvmr_harm$id.exposure)
}
uvmr_harm$exposure=gsub("\\|| ","",uvmr_harm$exposure)
uvmr_harm$exposure=gsub("id:","",uvmr_harm$exposure)
uvmr_harm$outcome="offspring birthweight"
for (j in 1:nrow(ao))
{
  uvmr_harm$id.exposure=gsub(ao$id[j], ao$trait[j], uvmr_harm$id.exposure)
  uvmr_harm$exposure=gsub(ao$id[j], ao$trait[j], uvmr_harm$exposure)
}
pdf("ukbb_snps/mr_scatter_plots.pdf")
mr_scatter_plot(mr(uvmr_harm, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median",
                                                                "mr_weighted_mode")),uvmr_harm)
dev.off()

pdf("ukbb_snps/mr_scatter_plots.pdf")
mr_scatter_plot(mr(uvmr_harm),uvmr_harm)
dev.off()
res_loo=mr_leaveoneout(uvmr_harm)
res_loo$exposure=gsub("\\|| ","",res_loo$exposure)
res_loo$exposure=gsub("id:","",res_loo$exposure)
res_loo$outcome="offspring birthweight"
for (j in 1:nrow(ao))
{
  res_loo$id.exposure=gsub(ao$id[j], ao$trait[j], res_loo$id.exposure)
  res_loo$exposure=gsub(ao$id[j], ao$trait[j], res_loo$exposure)
}
write.csv(res_loo,"ukbb_snps/res_loo.csv")
pdf("ukbb_snps/mr_res_loo.pdf")
mr_leaveoneout_plot(res_loo)
dev.off()

write.csv(estimates_all,"ukbb_snps/estimates_all_methods.csv")

## Kett
setwd("kett_snps_only")
nmr_metabolites=read_excel("nmr_metabolites_20210610.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
colnames(nmr_metabolites)[1]="GWAS_id"
nmr_metabolites$GWAS_id=paste0("met-d-", nmr_metabolites$GWAS_id)
ao=available_outcomes()
ao=ao[which(ao$year==2016&ao$author=="Kettunen"&ao$category=="Metabolites"&ao$population=="European"),]
ao_1=ao[which(ao$trait%in%nmr_metabolites$`Biomarker name`),]
name_mismatch=c("18:2, linoleic acid (LA)", "Free cholesterol", "3-hydroxybutyrate", "Apolipoprotein A-I", "22:6, docosahexaenoic acid", "Mono-unsaturated fatty acids",
                "Omega-7, omega-9 and saturated fatty acids", "Phosphatidylcholine and other cholines", "Total lipids in chylomicrons and largest VLDL particles", "Serum total triglycerides")
ao_1=rbind(ao_1, ao[which(ao$trait%in%name_mismatch),])
old=c("18:2, linoleic acid", "22:6, docosahexaenoic acid" ,"Omega-7, omega-9 and saturated fatty acids" ,"3-hydroxybutyrate" ,"Apolipoprotein A-I" ,"Concentration of chylomicrons and largest VLDL particles", 
      "Mono-unsaturated fatty acids" ,"Phosphatidylcholine and other cholines", "Serum total triglycerides", "Total phosphoglycerides", "Free cholesterol", "Serum total cholesterol",
      "Total lipids in chylomicrons and largest VLDL particles")
new=c("Linoleic acid", "Docosahexaenoic acid", "Saturated fatty acids", "3-Hydroxybutyrate", "Apolipoprotein A1", "Concentration of chylomicrons and extremely large VLDL particles",
      "Monounsaturated fatty acids", "Phosphatidylcholines", "Total triglycerides", "Phosphoglycerides", "Total free cholesterol", "Total cholesterol", "Total lipids in chylomicrons and extremely large VLDL")
uvmr_harm=read_csv("uvmr_harm.csv")
estimates_all=mr(uvmr_harm)
for (i in 1:length(old))
{
  ao_1$trait[grep(paste(old[i]),ao_1$trait)]=new[i]
}
for (j in 1:nrow(ao_1))
{
  estimates_all$id.exposure=gsub(ao_1$id[j], ao_1$trait[j], estimates_all$id.exposure)
}
write.csv(estimates_all,"/kett_snps_only/estimates_all_methods.csv")
