## MVMR estimates

# UKBB
library(dplyr);library(readr);library(TwoSampleMR);library(data.table);library(readxl);library(remotes);library(stringr)
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force=T)

ao=available_outcomes()
pheno_cov=read.delim("pcov_nmr_dat_all.txt")
nmr_metabolites_UKBB=read_xlsx("/Volumes/MRC-IEU-research/projects/ieu2/p6/106/working/data/data_v2/nmr_metabolites_20210610.xlsx")[,-1]
nmr_metabolites_UKBB=nmr_metabolites_UKBB[which(nmr_metabolites_UKBB$Include=="yes"),]
colnames(ao)[2]="Biomarker name"
ao1=ao[which(ao$author=="Borges CM"&ao$`Biomarker name`%in%nmr_metabolites_UKBB$`Biomarker name`),]
nmr_metabolites_UKBB=merge(nmr_metabolites_UKBB, ao1, by="Biomarker name")
nmr_metabolites_UKBB$id=sub("met-d-", "", nmr_metabolites_UKBB$id)
pheno_cov=pheno_cov[which(rownames(pheno_cov)%in%nmr_metabolites_UKBB$id), which(colnames(pheno_cov)%in%nmr_metabolites_UKBB$id)]

nmr_metabolites_UKBB$id=paste0("met-d-", nmr_metabolites_UKBB$id)

exposure_data=read_csv("mvmr_instruments.csv")[,-1]
snps=pull(exposure_data, SNP)
outcome_dat=read_outcome_data(snps=snps, filename="Maternal_Effect_European_meta_NG2019.txt",
                              snp_col="RSID", beta_col="beta", effect_allele_col="ea", other_allele_col="nea", eaf_col="eaf", pval_col="p")
a=c("Alanine","3-Hydroxybutyrate","Glutamine", "Glucose",  "Isoleucine", "Pyruvate")
b=sort(ao1$id[which(ao1$`Biomarker name`%in%a&ao1$author=="Borges CM")])
c=exposure_data[which(exposure_data$id.exposure%in%b),]

pheno_cov_tmp=pheno_cov[sort(rownames(pheno_cov))[which(paste0("met-d-",sort(rownames(pheno_cov)))%in%b)],
                        sort(colnames(pheno_cov))[which(paste0("met-d-",sort(colnames(pheno_cov)))%in%b)]]

mvmr_harm=mv_harmonise_data(c, outcome_dat, harmonise_strictness=2)
#mv_multiple(mvmr_harm, pval_threshold = 5e-08)

BGX=mvmr_harm$exposure_beta[,b]
seBGX=mvmr_harm$exposure_se[,b]
colnames(BGX)=gsub("met-d-","",colnames(BGX))
colnames(seBGX)=gsub("met-d-","",colnames(seBGX))

BGY=mvmr_harm$outcome_beta
seBGY=mvmr_harm$outcome_se

correlation=MVMR::phenocov_mvmr(pheno_cov_tmp, seBGX)
formatted=MVMR::format_mvmr(BXGs=BGX, seBXGs=seBGX, seBYG=seBGY, BYG=BGY, RSID=rownames(mvmr_harm$exposure_beta))

write_delim(formatted, "mvmr_formatted_data_CB.txt")

MVMR::ivw_mvmr(formatted, correlation)
