'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library(data.table)
library(rsq)

#################################################################################
###################### Read data: full population ###############################
#################################################################################
# read in relatedness file
rel <- fread("ukb_rel_a48818_s488117.dat",header=T,sep=" ")
rel <- as_tibble(rel)
rel$ID1 <- as.character(rel$ID1)
rel$ID2 <- as.character(rel$ID2)

# read in PCA covariates file
pc <- fread("all_pops_non_eur_pruned_within_pop_pc_covs.tsv",header=T,sep="\t")
pc <- as_tibble(pc)
pc$s <- as.character(pc$s)
table(pc$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6805    996   9109   2783 426901   1622

# read in sample ID file
sampleID <- fread("ukb48818bridge31063.txt",header=F,sep="\t")
sampleID <- as_tibble(sampleID)
colnames(sampleID) <- c("f.eid","s")
sampleID$f.eid <- as.character(sampleID$f.eid)
sampleID$s <- as.character(sampleID$s)

# read in withdrwan samples
withdrawn <- fread("withdrawn_w48818_2023-04-25.csv",header=F)
withdrawn <- as_tibble(withdrawn)
colnames(withdrawn) <- "f.eid"
withdrawn$f.eid <- as.character(withdrawn$f.eid)

# read in phenotype file, including 249 metabolites and covariates
pheno <- fread("pheno_249metabolites.tsv",header=T, sep="\t")
pheno <- as_tibble(pheno)
pheno$f.eid <- as.character(pheno$f.eid)

# read in phenotype file ukb34137.tab: contiansTownsend index f.189.0.0
ukb34137 <- fread("ukb34137.tab", header=TRUE, sep="\t")
ukb34137 <- as_tibble(ukb34137)

#################################################################################
####################### Read ancestry% data: AFR population #####################
#################################################################################
# read in admixture results (K=2)
q_metrics <- fread("/scratch/hx37930/project/admixture/AFR/result/supervised/k2/run1/studyAFR.2.Q",header=F,sep=" ")
q_metrics <- as_tibble(q_metrics)

fam <- fread("/scratch/hx37930/project/admixture/AFR/data/genotype/AFR_bed/ukb_allCHR.AFR.select.fam",header=F,sep=" ")
fam <- as_tibble(fam)

ancestry_prop <- cbind(fam[,1],q_metrics)
colnames(ancestry_prop) <- c("f.eid","pop1","pop2")
ancestry_prop <- as_tibble(ancestry_prop)
ancestry_prop$f.eid <- as.character(ancestry_prop$f.eid)

#################################################################################
#################### phenotype pre-processing: full population ##################
#################################################################################
# 1. genotyping array
pheno$f.22000.0.0=ifelse(pheno$f.22000.0.0>0,1,0)

# 2. Townsend index (TSI)
townsend <- subset(ukb34137,select=c(f.eid,f.189.0.0))
townsend$f.eid <- as.character(townsend$f.eid)
pheno <- pheno %>% left_join(townsend,by="f.eid")

# 3. Assessment_centres: code 10003 is missed in AFR participants,code 11023 is missed in CSA participants
pheno$f.54.0.0 <- mapvalues(as.character(pheno$f.54.0.0), c(11012, 11021,11011,11008,11003,11020,11005,11004,11018,11010,11016,11001,11017,11009,11013,11002,11007,11014,11006,11022,11023,10003), c("a11012","a11021","a11011","a11008","a11003","a11020",'a11005',"a11004","a11018","a11010","a11016","a11001","a11017","a11009","a11013","a11002","a11007","a11014","a11006","a11022","a11023","a10003"))
pheno$Assessment_centres_10003=ifelse(pheno$f.54.0.0=='a10003',1,0)
pheno$Assessment_centres_11001=ifelse(pheno$f.54.0.0=='a11001',1,0)
pheno$Assessment_centres_11002=ifelse(pheno$f.54.0.0=='a11002',1,0)
#non-England
pheno$Assessment_centres_11003=ifelse(pheno$f.54.0.0=='a11003',1,0)
#non-England
pheno$Assessment_centres_11004=ifelse(pheno$f.54.0.0=='a11004',1,0)
#non-England
pheno$Assessment_centres_11005=ifelse(pheno$f.54.0.0=='a11005',1,0)
pheno$Assessment_centres_11006=ifelse(pheno$f.54.0.0=='a11006',1,0)
pheno$Assessment_centres_11007=ifelse(pheno$f.54.0.0=='a11007',1,0)
pheno$Assessment_centres_11008=ifelse(pheno$f.54.0.0=='a11008',1,0)
pheno$Assessment_centres_11009=ifelse(pheno$f.54.0.0=='a11009',1,0)
pheno$Assessment_centres_11010=ifelse(pheno$f.54.0.0=='a11010',1,0)
pheno$Assessment_centres_11011=ifelse(pheno$f.54.0.0=='a11011',1,0)
pheno$Assessment_centres_11012=ifelse(pheno$f.54.0.0=='a11012',1,0)
pheno$Assessment_centres_11013=ifelse(pheno$f.54.0.0=='a11013',1,0)
pheno$Assessment_centres_11014=ifelse(pheno$f.54.0.0=='a11014',1,0)
pheno$Assessment_centres_11016=ifelse(pheno$f.54.0.0=='a11016',1,0)
pheno$Assessment_centres_11017=ifelse(pheno$f.54.0.0=='a11017',1,0)
pheno$Assessment_centres_11018=ifelse(pheno$f.54.0.0=='a11018',1,0)
pheno$Assessment_centres_11020=ifelse(pheno$f.54.0.0=='a11020',1,0)
pheno$Assessment_centres_11021=ifelse(pheno$f.54.0.0=='a11021',1,0)
#non-England
pheno$Assessment_centres_11022=ifelse(pheno$f.54.0.0=='a11022',1,0)
#non-England
pheno$Assessment_centres_11023=ifelse(pheno$f.54.0.0=='a11023',1,0)

# 4. merge sampleID; match s (PCA) with eid (phenotype)
sample_pca <- unique(subset(pc,select=c(s, pop)))
table(sample_pca$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6752    996   9065   2783 426901   1622

sampleID_merge <- sample_pca %>% inner_join(sampleID,by="s")
table(sampleID_merge$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6751    996   9064   2782 426881   1622

pheno_mergeID <- sampleID_merge %>% inner_join(pheno,by="f.eid")
table(pheno_mergeID$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6749    994   9063   2782 426810   1622

write.table(pheno_mergeID, "pheno_249metabolites.clean.fullPop.tsv",sep="\t",col.names=T,row.names=F,quote=F)

##################################################################################
# phenotype QC (full population): 1) remove individuals with withdrawn info; 2) remove individuals with 3rd-degree relatives or closer; 3) mismatched information between phenotypic and genetix sex; 4) Sex chromosomes aneuploidy; 5) outlier for het and missing genotype rate
##################################################################################
# 1. remove samples with withdrawn info (full population)
pheno_rmWithdrwan <- subset(pheno_mergeID, ! f.eid %in% withdrawn$f.eid)
table(pheno_rmWithdrwan$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6749    994   9063   2782 426810   1622

# 2. remove individuals with 3rd-degree relatives or closer
# create a list of samples need to be removed due to 3rd-degree relatives or closer (full population)
rmID_Kinship <- unique(rel[rel$Kinship >= 0.0442,2])
colnames(rmID_Kinship) <- "f.eid"

# remove individuals with 3rd-degree relatives or closer
pheno_rmKinship <- subset(pheno_rmWithdrwan, ! f.eid %in% rmID_Kinship$f.eid)
table(pheno_rmKinship$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6174    949   8356   2668 351650   1536

# 3. mismatched information between phenotypic (f.31.0.0) and genetic sex (f.22001.0.0)
pheno_rmMismatchedSex <-  pheno_rmKinship[pheno_rmKinship$f.31.0.0 == pheno_rmKinship$f.22001.0.0,]
table(pheno_rmMismatchedSex$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6172    949   8354   2666 351390   1533

# 4. Sex chromosomes aneuploidy
pheno_rmSexAneuploidy <- pheno_rmMismatchedSex %>%filter(is.na(f.22019.0.0)) 
table(pheno_rmSexAneuploidy$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6168    947   8347   2666 351047   1533

# 5. outlier for het and missing genotype rate
pheno_rmHetMiss <- pheno_rmSexAneuploidy %>%filter(is.na(f.22027.0.0))
table(pheno_rmHetMiss$pop)
#AFR    AMR    CSA    EAS    EUR    MID
#6167    944   8296   2648 350602   1530

write.table(pheno_rmHetMiss, "pheno_249metabolites.clean.fullPop.QCed.tsv",sep="\t",col.names=T,row.names=F,quote=F)

#################################################################################
########################### Subset AFR population ###############################
#################################################################################
# Combine ancestry% with QCed phenotype
df_combo <- ancestry_prop %>% left_join(pheno_rmHetMiss,by="f.eid")
df_combo <- as_tibble(df_combo)

write.table(df_combo, "pheno_249metabolites.clean.QCed.AFR.tsv",sep="\t",col.names=T,row.names=F,quote=F)

##################################################################################
####################  rank-based inverse transformation ##########################
##################################################################################
inversenormal <- function(x) { 
    # inverse normal if you have missing data 
    return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))) 
} 

df_combo_INT <- df_combo
for (i in 17:265){
    df_combo_INT[,i] <- inversenormal(df_combo[,i])
}

write.table(df_combo, "pheno_249metabolites.clean.QCed.AFR.INT.tsv",sep="\t",col.names=T,row.names=F,quote=F)
