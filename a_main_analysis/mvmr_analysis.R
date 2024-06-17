## MVMR estimates

# UKBB
library(dplyr);library(readr);library(TwoSampleMR);library(data.table);library(readxl);library(remotes);library(stringr)
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)

ao=available_outcomes()
pheno_cov=read.delim("pcov_nmr_dat_all.txt")
nmr_metabolites_UKBB=read_xlsx("nmr_metabolites_20210610.xlsx")[,-1]
nmr_metabolites_UKBB=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$Include=="yes"),]
colnames(ao)[2]="Biomarker name"
ao1=ao[which(ao$author=="Borges CM"&ao$`Biomarker name`%in%nmr_metabolites_UKBB$`Biomarker name`),]
nmr_metabolites_UKBB=merge(nmr_metabolites_UKBB, ao1, by="Biomarker name")
nmr_metabolites_UKBB$id=sub("met-d-", "", nmr_metabolites_UKBB$id)
pheno_cov=pheno_cov[which(rownames(pheno_cov)%in%nmr_metabolites_UKBB$id), which(colnames(pheno_cov)%in%nmr_metabolites_UKBB$id)]

nmr_metabolites_UKBB$id=paste0("met-d-", nmr_metabolites_UKBB$id)

a=c("Alanine","3-Hydroxybutyrate","Glutamine", "Glucose",  "Isoleucine", "Pyruvate")
b=sort(ao1$id[which(ao1$`Biomarker name`%in%a&ao1$author=="Borges CM")])
pheno_cov_tmp=pheno_cov[sort(rownames(pheno_cov))[which(paste0("met-d-",sort(rownames(pheno_cov)))%in%b)],
                        sort(colnames(pheno_cov))[which(paste0("met-d-",sort(colnames(pheno_cov)))%in%b)]]
pheno_cov_tmp=as.matrix(pheno_cov_tmp)

exposure_data=mv_extract_exposures(b, clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token = ieugwasr::check_access_token(),
                                   find_proxies = TRUE, force_server = FALSE, pval_threshold = 5e-08, pop = "EUR")
snps=pull(exposure_data, SNP)
outcome_dat=read_outcome_data(snps=snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")

mvmr_harm=mv_harmonise_data(exposure_data, outcome_dat, harmonise_strictness=2)

BGX=mvmr_harm$exposure_beta
seBGX=mvmr_harm$exposure_se
colnames(BGX)=gsub("met-d-","",colnames(BGX))
colnames(seBGX)=gsub("met-d-","",colnames(seBGX))

BGY=mvmr_harm$outcome_beta
seBGY=mvmr_harm$outcome_se

correlation=MVMR::phenocov_mvmr(pheno_cov_tmp, seBGX)
formatted=MVMR::format_mvmr(BXGs=BGX, seBXGs=seBGX, seBYG=seBGY, BYG=BGY, RSID=rownames(mvmr_harm$exposure_beta))
MVMR::ivw_mvmr(formatted, correlation)
MVMR::strength_mvmr(formatted, correlation)
MVMR::qhet_mvmr(formatted, pheno_cov_tmp)
MVMR::pleiotropy_mvmr(formatted, correlation)

write_delim(as.data.frame(MVMR::ivw_mvmr(formatted, correlation)), "ukbb_snps/ivw_mvmr_ukbb_est.txt")
write_delim(as.data.frame(MVMR::strength_mvmr(formatted, correlation)), "ukbb_snps/strength_mvmr_ukbb_est.txt")
write_delim(as.data.frame(MVMR::qhet_mvmr(formatted, pheno_cov_tmp)), "ukbb_snps/qhet_mvmr_ukbb_est.txt")
write_delim(as.data.frame(MVMR::pleiotropy_mvmr(formatted, correlation)), "ukbb_snps/pleiotropy_mvmr_ukbb_est.txt")


## model excluding glucose and glutamine

a=c("Alanine","3-Hydroxybutyrate", "Isoleucine", "Pyruvate")
b=sort(ao1$id[which(ao1$`Biomarker name`%in%a&ao1$author=="Borges CM")])
pheno_cov_tmp=pheno_cov[sort(rownames(pheno_cov))[which(paste0("met-d-",sort(rownames(pheno_cov)))%in%b)],
                        sort(colnames(pheno_cov))[which(paste0("met-d-",sort(colnames(pheno_cov)))%in%b)]]
pheno_cov_tmp=as.matrix(pheno_cov_tmp)

exposure_data=mv_extract_exposures(b, clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2, access_token = ieugwasr::check_access_token(),
                                   find_proxies = TRUE, force_server = FALSE, pval_threshold = 5e-08, pop = "EUR")
snps=pull(exposure_data, SNP)
outcome_dat=read_outcome_data(snps=snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")

mvmr_harm=mv_harmonise_data(exposure_data, outcome_dat, harmonise_strictness=2)

BGX=mvmr_harm$exposure_beta
seBGX=mvmr_harm$exposure_se
colnames(BGX)=gsub("met-d-","",colnames(BGX))
colnames(seBGX)=gsub("met-d-","",colnames(seBGX))

BGY=mvmr_harm$outcome_beta
seBGY=mvmr_harm$outcome_se

correlation=MVMR::phenocov_mvmr(pheno_cov_tmp, seBGX)
formatted=MVMR::format_mvmr(BXGs=BGX, seBXGs=seBGX, seBYG=seBGY, BYG=BGY, RSID=rownames(mvmr_harm$exposure_beta))

write_delim(as.data.frame(MVMR::ivw_mvmr(formatted, correlation)), "ukbb_snps/ivw_mvmr_ukbb_est_SA.txt")
write_delim(as.data.frame(MVMR::strength_mvmr(formatted, correlation)), "ukbb_snps/strength_mvmr_ukbb_est_SA.txt")
write_delim(as.data.frame(MVMR::qhet_mvmr(formatted, pheno_cov_tmp)), "ukbb_snps/qhet_mvmr_ukbb_est_SA.txt")
write_delim(as.data.frame(MVMR::pleiotropy_mvmr(formatted, correlation)), "ukbb_snps/pleiotropy_mvmr_ukbb_est_SA.txt")

## Kettunen (15 SNPs)

nmr_metabolites=read_excel("nmr_metabolites_20210610.xlsx")
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
colnames(nmr_metabolites)[1]="GWAS_id"
nmr_metabolites$GWAS_id=paste0("met-d-", nmr_metabolites$GWAS_id)
nmr_metabolites=nmr_metabolites[order(nmr_metabolites$GWAS_id),]
ao=available_outcomes()
ao=ao[which(ao$year==2016&ao$author=="Kettunen"&ao$category=="Metabolites"&ao$population=="European"),]
ao1=ao[which(ao$trait%in%nmr_metabolites$`Biomarker name`),]
name_mismatch=c("18:2, linoleic acid (LA)", "Free cholesterol", "3-hydroxybutyrate", "Apolipoprotein A-I", "22:6, docosahexaenoic acid", "Mono-unsaturated fatty acids",
                "Omega-7, omega-9 and saturated fatty acids", "Phosphatidylcholine and other cholines", "Total lipids in chylomicrons and largest VLDL particles", "Serum total triglycerides")
ao1=rbind(ao1, ao[which(ao$trait%in%name_mismatch),])
list_of_datasets=as.data.frame(ao1)

a=c("Alanine","3-hydroxybutyrate","Glutamine", "Glucose",  "Isoleucine", "Pyruvate")
b=sort(ao1$id[which(ao1$trait%in%a&ao1$author=="Kettunen")])

mvmr_instruments=mv_extract_exposures(b,pval_threshold=5e-8,clump_r2=0.01)
ket_snps=pull(mvmr_instruments, SNP)
expdat=extract_outcome_data(snps = ket_snps, outcome = nmr_metabolites$GWAS_id[which(nmr_metabolites$`Biomarker name`%in%c("Alanine","3-Hydroxybutyrate","Glutamine", "Glucose",  "Isoleucine", "Pyruvate"))])
names(expdat)=gsub("outcome", "exposure", names(expdat))

outcome_dat=read_outcome_data(snps=ket_snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")

mvmr_harm=mv_harmonise_data(expdat, outcome_dat, harmonise_strictness=2)

pheno_cov=read.delim("pcov_nmr_dat_all.txt")
pheno_cov=pheno_cov[which(paste0("met-d-",colnames(pheno_cov))%in%mvmr_harm$expname$id.exposure), which(paste0("met-d-",colnames(pheno_cov))%in%mvmr_harm$expname$id.exposure)]
pheno_cov_tmp=pheno_cov[sort(rownames(pheno_cov))[which(paste0("met-d-",sort(rownames(pheno_cov)))%in%nmr_metabolites$GWAS_id)],
                        sort(colnames(pheno_cov))[which(paste0("met-d-",sort(colnames(pheno_cov)))%in%nmr_metabolites$GWAS_id)]]
pheno_cov_tmp=as.matrix(pheno_cov_tmp)


BGX=mvmr_harm$exposure_beta
seBGX=mvmr_harm$exposure_se
colnames(BGX)=gsub("met-d-","",colnames(BGX))
colnames(seBGX)=gsub("met-d-","",colnames(seBGX))

BGY=mvmr_harm$outcome_beta
seBGY=mvmr_harm$outcome_se

correlation=MVMR::phenocov_mvmr(pheno_cov_tmp, seBGX)
formatted=MVMR::format_mvmr(BXGs=BGX, seBXGs=seBGX, seBYG=seBGY, BYG=BGY, RSID=rownames(mvmr_harm$exposure_beta))
MVMR::ivw_mvmr(formatted, correlation)
MVMR::strength_mvmr(formatted, correlation)
MVMR::qhet_mvmr(formatted, pheno_cov_tmp)
MVMR::pleiotropy_mvmr(formatted, correlation)

write_delim(as.data.frame(MVMR::ivw_mvmr(formatted, correlation)), "kett_snps_only/ivw_mvmr_ket_est.txt")
write_delim(as.data.frame(MVMR::strength_mvmr(formatted, correlation)), "kett_snps_only/strength_mvmr_ket_est.txt")
write_delim(as.data.frame(MVMR::qhet_mvmr(formatted, pheno_cov_tmp)), "kett_snps_only/qhet_mvmr_ket_est.txt")
write_delim(as.data.frame(MVMR::pleiotropy_mvmr(formatted, correlation)), "kett_snps_only/pleiotropy_mvmr_ket_est.txt")

## SA (8 SNPs)

a=c("Alanine","3-hydroxybutyrate","Isoleucine", "Pyruvate")
b=sort(ao$id[which(ao$trait%in%a&ao$author=="Kettunen")])

mvmr_instruments=mv_extract_exposures(b,pval_threshold=5e-8,clump_r2=0.01)
ket_snps=pull(mvmr_instruments, SNP)
expdat=extract_outcome_data(snps = ket_snps, outcome = nmr_metabolites$GWAS_id[which(nmr_metabolites$`Biomarker name`%in%c("Alanine","3-Hydroxybutyrate", "Isoleucine", "Pyruvate"))])
names(expdat)=gsub("outcome", "exposure", names(expdat))

outcome_dat=read_outcome_data(snps=ket_snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")

mvmr_harm=mv_harmonise_data(expdat, outcome_dat, harmonise_strictness=2)

pheno_cov=pheno_cov[which(paste0("met-d-",colnames(pheno_cov))%in%mvmr_harm$expname$id.exposure), which(paste0("met-d-",colnames(pheno_cov))%in%mvmr_harm$expname$id.exposure)]
pheno_cov_tmp=pheno_cov[sort(rownames(pheno_cov))[which(paste0("met-d-",sort(rownames(pheno_cov)))%in%nmr_metabolites$GWAS_id)],
                        sort(colnames(pheno_cov))[which(paste0("met-d-",sort(colnames(pheno_cov)))%in%nmr_metabolites$GWAS_id)]]
pheno_cov_tmp=as.matrix(pheno_cov_tmp)

BGX=mvmr_harm$exposure_beta
seBGX=mvmr_harm$exposure_se
colnames(BGX)=gsub("met-d-","",colnames(BGX))
colnames(seBGX)=gsub("met-d-","",colnames(seBGX))

BGY=mvmr_harm$outcome_beta
seBGY=mvmr_harm$outcome_se

correlation=MVMR::phenocov_mvmr(pheno_cov_tmp, seBGX)
formatted=MVMR::format_mvmr(BXGs=BGX, seBXGs=seBGX, seBYG=seBGY, BYG=BGY, RSID=rownames(mvmr_harm$exposure_beta))

write_delim(as.data.frame(MVMR::ivw_mvmr(formatted, correlation)), "kett_snps_only/ivw_mvmr_ket_est_SA.txt")
write_delim(as.data.frame(MVMR::strength_mvmr(formatted, correlation)), "kett_snps_only/strength_mvmr_ket_est_SA.txt")
write_delim(as.data.frame(MVMR::qhet_mvmr(formatted, pheno_cov_tmp)), "kett_snps_only/qhet_mvmr_ket_est_SA.txt")
write_delim(as.data.frame(MVMR::pleiotropy_mvmr(formatted, correlation)), "kett_snps_only/pleiotropy_mvmr_ket_est_SA.txt")

