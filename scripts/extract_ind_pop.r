library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)

# read in PanUKBB file
pan <- fread("all_pops_non_eur_pruned_within_pop_pc_covs.tsv", header=T, sep="\t")
pan <- pan[!duplicated(pan$s),]

# read in PanUKBB ID file
pan_id <- read.table("ukb48818bridge31063.txt", header=F, sep="\t")
colnames(pan_id)=c("FID","s") # or c("f.eid", "s")

# read in genotype .fam file
fam <- fread("ukb_allCHR.fam", header=F, sep="\t")
colnames(fam) <- c("FID","IID","IID_fa","IID_ma","Sex","PhenoValue")

# combine pan with pan_id
pan <- pan %>% inner_join(pan_id, by= "s")
table(pan$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6751    996   9064   2782 426881   1622

# subset csa and afr individuals
csa <- pan[pan$pop=="CSA",]
afr <- pan[pan$pop=="AFR",]


# combine csa with genotype .fam file
keep_csa <- csa %>% inner_join(fam, by= "FID")
# combine afr with genotype .fam file
keep_afr <- afr %>% inner_join(fam, by= "FID")


csa_fam <- subset(keep_csa,select=c(FID,IID,IID_fa,IID_ma,Sex,PhenoValue))
afr_fam <- subset(keep_afr,select=c(FID,IID,IID_fa,IID_ma,Sex,PhenoValue))

write.table(csa_fam,"keep_sample.CSA.fam",col.names=F,row.names=F,quote=F,sep="\t")
write.table(afr_fam,"keep_sample.AFR.fam",col.names=F,row.names=F,quote=F,sep="\t")
