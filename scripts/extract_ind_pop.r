library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)

# read in PanUKBB file
pan <- fread("all_pops_non_eur_pruned_within_pop_pc_covs.tsv", header=T, sep="\t")
dim(pan)
#[1] 448216     28
pan <- pan[!duplicated(pan$s),]
dim(pan)
#[1] 448119     28

# read in PanUKBB ID file
pan_id <- read.table("ukb48818bridge31063.txt", header=F, sep="\t")
colnames(pan_id)=c("FID","s") # or c("f.eid", "s")
dim(pan_id)
#[1] 502460      2

# read in genotype .fam file
fam <- fread("ukb_allCHR.fam", header=F, sep="\t")
colnames(fam) <- c("FID","IID","IID_fa","IID_ma","Sex","PhenoValue")
dim(fam)
#[1] 487409      6

# combine pan with pan_id
pan <- pan %>% inner_join(pan_id, by= "s")
dim(pan)
#[1] 448096     29
table(pan$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6751    996   9064   2782 426881   1622

# subset csa and afr individuals
csa <- pan[pan$pop=="CSA",]
dim(csa)
#[1] 9064   29
afr <- pan[pan$pop=="AFR",]
dim(afr)
#[1] 6751   29

# combine csa with genotype .fam file
keep_csa <- csa %>% inner_join(fam, by= "FID")
dim(keep_csa)
#[1] 9013   34
keep_afr <- afr %>% inner_join(fam, by= "FID")
dim(keep_afr)
#[1] 6750   34

csa_fam <- subset(keep_csa,select=c(FID,IID,IID_fa,IID_ma,Sex,PhenoValue))
afr_fam <- subset(keep_afr,select=c(FID,IID,IID_fa,IID_ma,Sex,PhenoValue))

write.table(csa_fam,"keep_sample.CSA.fam",col.names=F,row.names=F,quote=F,sep="\t")
write.table(afr_fam,"keep_sample.AFR.fam",col.names=F,row.names=F,quote=F,sep="\t")
