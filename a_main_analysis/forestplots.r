setwd("ukbb_snps")
devtools::install_github("NightingaleHealth/ggforestplot")
library(readr); library(ggplot2); library(data.table);library(TwoSampleMR);library(ggforestplot);library(readxl)

# UKBB plots

nmr_metabolites=read_xlsx("nmr_metabolites_20210610.xlsx")[,-1]
nmr_metabolites=nmr_metabolites[which(nmr_metabolites$Include=="yes"),]
ao=available_outcomes()
colnames(ao)[2]="Biomarker name"
ao1=ao[which(ao$author=="Borges CM"&ao$`Biomarker name`%in%nmr_metabolites$`Biomarker name`),]
nmr_metabolites=merge(nmr_metabolites, ao1, by="Biomarker name")
nmr_metabolites$id=sub("met-d-", "", nmr_metabolites$id)
nmr_metabolites=nmr_metabolites[order(nmr_metabolites$id),]
ukbb_uvmr=read_csv("ukbb_snps/naive_uvmr_results.csv")[,-1]
ukbb_mvmr=read.csv("ukbb_snps/ivw_mvmr_ukbb_est.txt", sep="")

colnames(ukbb_mvmr)=c("b", "se", "tval", "pval")
ukbb_mvmr$exposure=c("Alanine","3-Hydroxybutyrate","Glutamine", "Glucose",  "Isoleucine", "Pyruvate")
ukbb_mvmr$id=paste0("met-d-",nmr_metabolites$id[which(nmr_metabolites$`Biomarker name`%in%ukbb_mvmr$exposure)])

ukbb_uvmr$group=ukbb_mvmr$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  if (i%in%ukbb_uvmr$trait)
  {
    ukbb_uvmr$group[which(ukbb_uvmr$trait==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
  if (i%in%ukbb_mvmr$exposure)
  {
    ukbb_mvmr$group[which(ukbb_mvmr$exposure==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}
ukbb_uvmr=ukbb_uvmr[which(ukbb_uvmr$trait%in%ukbb_mvmr$exposure),]
ukbb_uvmr=ukbb_uvmr[order(ukbb_uvmr$trait),]

ukbb_mvmr$Dataset="MVMR"
ukbb_uvmr$Dataset="UVMR"

dataframe=rbind(ukbb_uvmr[,intersect(colnames(ukbb_uvmr), colnames(ukbb_mvmr))],
                ukbb_mvmr[,intersect(colnames(ukbb_uvmr), colnames(ukbb_mvmr))])
dataframe=dataframe[order(dataframe$id),]
dataframe$trait=NA

for (i in 1:nrow(dataframe))
{
  dataframe$trait[i]=nmr_metabolites$`Biomarker name`[which(nmr_metabolites$id==gsub("met-d-","",dataframe$id[i]))]
}

colnames(dataframe)[grep("b",colnames(dataframe))]="Effect estimate (95% CI)"

png("forest_plot_uvmr_mvmr_ukbb_ungrouped.png", width = 1000, height = 1400)
ggforestplot::forestplot(df=dataframe, estimate=`Effect estimate (95% CI)`, pvalue = pval, name=trait, 
                         colour=Dataset, shape=Dataset)+ #ggforce::facet_col(facets=~group, scales="free_y", space="free")+
  theme(legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        text=element_text(size=20),) 
dev.off()

## SA UKBB

ukbb_mvmr=read.csv("ukbb_snps/ivw_mvmr_ukbb_est_SA.txt", sep="")

colnames(ukbb_mvmr)=c("b", "se", "tval", "pval")
ukbb_mvmr$exposure=c("Alanine","3-Hydroxybutyrate", "Isoleucine", "Pyruvate")
ukbb_mvmr$id=paste0("met-d-",nmr_metabolites$id[which(nmr_metabolites$`Biomarker name`%in%ukbb_mvmr$exposure)])
ukbb_uvmr=ukbb_uvmr[which(ukbb_uvmr$trait%in%ukbb_mvmr$exposure),]
ukbb_mvmr$Dataset="MVMR"

ukbb_mvmr$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{

  if (i%in%ukbb_mvmr$exposure)
  {
    ukbb_mvmr$group[which(ukbb_mvmr$exposure==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}

dataframe=rbind(ukbb_uvmr[,intersect(colnames(ukbb_uvmr), colnames(ukbb_mvmr))],
                ukbb_mvmr[,intersect(colnames(ukbb_uvmr), colnames(ukbb_mvmr))])
dataframe=dataframe[order(dataframe$id),]
dataframe$trait=NA

for (i in 1:nrow(dataframe))
{
  dataframe$trait[i]=nmr_metabolites$`Biomarker name`[which(nmr_metabolites$id==gsub("met-d-","",dataframe$id[i]))]
}
colnames(dataframe)[grep("b",colnames(dataframe))]="Effect estimate (95% CI)"

png("forest_plot_uvmr_mvmr_ukbb_SA_ungrouped.png", width = 1000, height = 1400)
ggforestplot::forestplot(df=dataframe, estimate=`Effect estimate (95% CI)`, pvalue = pval, name=trait, 
                         colour=Dataset, shape=Dataset)+#ggforce::facet_col(facets=~group, scales="free_y", space="free")+
  theme(legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        text=element_text(size=20),) 
dev.off()

## Kettunen plots
setwd("kett_snps_only")

kett_uvmr=read_csv("kett_snps_only/naive_uvmr_results.csv")[,-1]
kett_mvmr=read.csv("kett_snps_only/ivw_mvmr_ket_est.txt", sep="")

old=c("18:2, linoleic acid", "22:6, docosahexaenoic acid" ,"Omega-7, omega-9 and saturated fatty acids" ,"3-hydroxybutyrate" ,"Apolipoprotein A-I" ,"Concentration of chylomicrons and largest VLDL particles", 
      "Mono-unsaturated fatty acids" ,"Phosphatidylcholine and other cholines", "Serum total triglycerides", "Total phosphoglycerides", "Free cholesterol", "Serum total cholesterol",
      "Total lipids in chylomicrons and largest VLDL particles")
new=c("Linoleic acid", "Docosahexaenoic acid", "Saturated fatty acids", "3-Hydroxybutyrate", "Apolipoprotein A1", "Concentration of chylomicrons and extremely large VLDL particles",
      "Monounsaturated fatty acids", "Phosphatidylcholines", "Total triglycerides", "Phosphoglycerides", "Total free cholesterol", "Total cholesterol", "Total lipids in chylomicrons and extremely large VLDL")

for (i in 1:length(old))
{
  kett_uvmr$trait[grep(paste(old[i]),kett_uvmr$trait)]=new[i]
}

colnames(kett_mvmr)=c("b", "se", "tval", "pval")
kett_mvmr$exposure=c("Alanine","3-Hydroxybutyrate","Glutamine", "Glucose",  "Isoleucine", "Pyruvate")
kett_mvmr$id=paste0("met-d-",nmr_metabolites$id[which(nmr_metabolites$`Biomarker name`%in%kett_mvmr$exposure)])

kett_uvmr$group=kett_mvmr$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  if (i%in%kett_uvmr$trait)
  {
    kett_uvmr$group[which(kett_uvmr$trait==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
  if (i%in%kett_mvmr$exposure)
  {
    kett_mvmr$group[which(kett_mvmr$exposure==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}
kett_uvmr=kett_uvmr[which(kett_uvmr$trait%in%kett_mvmr$exposure),]
kett_uvmr=kett_uvmr[order(kett_uvmr$trait),]

kett_mvmr$Dataset="MVMR"
kett_uvmr$Dataset="UVMR"

dataframe=rbind(kett_uvmr[,intersect(colnames(kett_uvmr), colnames(kett_mvmr))],
                kett_mvmr[,intersect(colnames(kett_uvmr), colnames(kett_mvmr))])
dataframe=dataframe[order(dataframe$id),]
dataframe$trait=NA

for (i in 1:nrow(dataframe))
{
  dataframe$trait[i]=ao$`Biomarker name`[which(ao$id==dataframe$id[i])]
}

colnames(dataframe)[grep("b",colnames(dataframe))]="Effect estimate (95% CI)"
dataframe$trait=gsub("3-hy","3-Hy", dataframe$trait)

png("forest_plot_uvmr_mvmr_kett_ungrouped.png", width = 1000, height = 1400)
ggforestplot::forestplot(df=dataframe, estimate=`Effect estimate (95% CI)`, pvalue = pval, name=trait, 
                         colour=Dataset, shape=Dataset)+#ggforce::facet_col(facets=~group, scales="free_y", space="free")+
  theme(legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        text=element_text(size=20),) 
dev.off()

## SA kett

kett_mvmr=read.csv("kett_snps_only/ivw_mvmr_ket_est_SA.txt", sep="")

colnames(kett_mvmr)=c("b", "se", "tval", "pval")
kett_mvmr$exposure=c("Alanine","3-Hydroxybutyrate",  "Isoleucine", "Pyruvate")
kett_mvmr$id=paste0("met-d-",nmr_metabolites$id[which(nmr_metabolites$`Biomarker name`%in%kett_mvmr$exposure)])

kett_mvmr$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  if (i%in%kett_mvmr$exposure)
  {
    kett_mvmr$group[which(kett_mvmr$exposure==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}
kett_mvmr$group=NA
for (i in nmr_metabolites$`Biomarker name`)
{
  
  if (i%in%kett_mvmr$exposure)
  {
    kett_mvmr$group[which(kett_mvmr$exposure==i)]=nmr_metabolites$Group[which(nmr_metabolites$`Biomarker name`==i)]
  }
}
kett_mvmr$Dataset="MVMR"
dataframe=rbind(kett_uvmr[,intersect(colnames(kett_uvmr), colnames(kett_mvmr))],
                kett_mvmr[,intersect(colnames(kett_uvmr), colnames(kett_mvmr))])
dataframe=dataframe[order(dataframe$id),]
dataframe$trait=NA
for (i in 1:nrow(dataframe))
{
  dataframe$trait[i]=ao$`Biomarker name`[which(ao$id==dataframe$id[i])]
}
colnames(dataframe)[grep("b",colnames(dataframe))]="Effect estimate (95% CI)"
dataframe$trait=gsub("3-hy","3-Hy", dataframe$trait)

png("forest_plot_uvmr_mvmr_kett_SA.png", width = 1000, height = 1400)
ggforestplot::forestplot(df=dataframe, estimate=`Effect estimate (95% CI)`, pvalue = pval, name=trait, 
                         colour=Dataset, shape=Dataset)+ggforce::facet_col(facets=~group, scales="free_y", space="free")+
  theme(legend.title = element_text(size = 24),
        legend.text = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        text=element_text(size=20),) 
dev.off()
