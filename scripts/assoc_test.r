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
df_combo <- fread("pheno_249metabolites.clean.QCed.AFR.INT.tsv")

##################################################################
# metabolite ~ ancestry% association
##################################################################
pops <- c("pop1","pop2")

# model1: metabolite ~ ancestry% + sex (f.31.0.0) + age (f.21003.0.0) + Genotyping_array (f.22000.0.0) + Assessment_centres
header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","pop_rsq","sex_rsq","age_rsq")
write.table(t(as.data.frame(header)),file="result_Model1Ancestry.k2.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
for (i in 17:265){
	for (j in 1:length(pops)){
		sampleSize <- sum(!is.na(df_combo[,i]))
		fitModel <- glm(unlist(df_combo[,i])~unlist(df_combo[,pops[j]])+f.31.0.0+f.21003.0.0+f.22000.0.0+Assessment_centres_11023+Assessment_centres_11001+Assessment_centres_11002+Assessment_centres_11003+Assessment_centres_11004+Assessment_centres_11005+Assessment_centres_11006+Assessment_centres_11007+Assessment_centres_11008+Assessment_centres_11009+Assessment_centres_11010+Assessment_centres_11011+Assessment_centres_11012+Assessment_centres_11013+Assessment_centres_11014+Assessment_centres_11016+Assessment_centres_11017+Assessment_centres_11018+Assessment_centres_11020+Assessment_centres_11021+Assessment_centres_11022,data=df_combo)
		ml_summary <- summary(fitModel)
		pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
		sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
		age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
		rsqValue <- rsq.partial(fitModel, adj = TRUE)
		pop_rsq <- rsqValue$partial.rsq[1]
		sex_rsq <- rsqValue$partial.rsq[2]
		age_rsq <- rsqValue$partial.rsq[3]
		result <- t(as.data.frame(c(sampleSize,colnames(df_combo[i]),pops[j],pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,pop_rsq,sex_rsq,age_rsq)))
		write.table(result,file="result_Model1Ancestry.k2.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
		pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"
	}
}

# model 2: add townsend index (sensitivity analysis)
# metabolite ~ ancestry% + sex (f.31.0.0) + age (f.21003.0.0) + Genotyping_array (f.22000.0.0) + Townsend index (f.189.0.0) + Assessment_centres
header <- c("sampleSize","phenotype","population","pop_coef","pop_se","pop_pvalue","sex_coef","sex_se","sex_pvalue","age_coef","age_se","age_pvalue","TSI_coef","TSI_se","TSI_pvalue","pop_rsq","sex_rsq","age_rsq","TSI_rsq")
write.table(t(as.data.frame(header)),file="result_Model2AncestryAddTSI.k2.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
for (i in 17:265){
        for (j in 1:length(pops)){
                sampleSize <- sum(!is.na(df_combo[,i]))
		fitModel <- glm(unlist(df_combo[,i])~unlist(df_combo[,pops[j]])+f.31.0.0+f.21003.0.0+f.189.0.0+f.22000.0.0+Assessment_centres_11023+Assessment_centres_11001+Assessment_centres_11002+Assessment_centres_11003+Assessment_centres_11004+Assessment_centres_11005+Assessment_centres_11006+Assessment_centres_11007+Assessment_centres_11008+Assessment_centres_11009+Assessment_centres_11010+Assessment_centres_11011+Assessment_centres_11012+Assessment_centres_11013+Assessment_centres_11014+Assessment_centres_11016+Assessment_centres_11017+Assessment_centres_11018+Assessment_centres_11020+Assessment_centres_11021+Assessment_centres_11022,data=df_combo)
                ml_summary <- summary(fitModel)
                pop_coef <- ml_summary$coefficients[2,1]; pop_se <- ml_summary$coefficients[2,2]; pop_pvalue <- ml_summary$coefficients[2,4]
                sex_coef <- ml_summary$coefficients[3,1]; sex_se <- ml_summary$coefficients[3,2]; sex_pvalue <- ml_summary$coefficients[3,4]
                age_coef <- ml_summary$coefficients[4,1]; age_se <- ml_summary$coefficients[4,2]; age_pvalue <- ml_summary$coefficients[4,4]
		TSI_coef <- ml_summary$coefficients[5,1]; TSI_se <- ml_summary$coefficients[5,2]; TSI_pvalue <- ml_summary$coefficients[5,4]
                rsqValue <- rsq.partial(fitModel, adj = TRUE)
                pop_rsq <- rsqValue$partial.rsq[1]
                sex_rsq <- rsqValue$partial.rsq[2]
                age_rsq <- rsqValue$partial.rsq[3]
		TSI_rsq <- rsqValue$partial.rsq[4]
                result <- t(as.data.frame(c(sampleSize,colnames(df_combo[i]),pops[j],pop_coef,pop_se,pop_pvalue,sex_coef,sex_se,sex_pvalue,age_coef,age_se,age_pvalue,TSI_coef,TSI_se,TSI_pvalue,pop_rsq,sex_rsq,age_rsq,TSI_rsq)))
		write.table(result,file="result_Model2AncestryAddTSI.k2.txt",col.names = FALSE, append = TRUE,row.names = F, quote = FALSE, na = "-",sep='\t')
                pop_coef="NA"; pop_se="NA";pop_pvalue="NA";sex_coef="NA";sex_se="NA";sex_pvalue="NA";age_coef="NA";age_se="NA";age_pvalue="NA";TSI_coef <- "NA";TSI_se <- "NA";TSI_pvalue <- "NA";pop_rsq <- "NA"; sex_rsq <- "NA"; age_rsq <- "NA"; TSI_rsq<- "NA"
        }
}

