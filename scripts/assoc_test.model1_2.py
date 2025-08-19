#!bin/python3

# load modules
import pandas as pd
import numpy as np
from scipy.stats import rankdata, norm, shapiro
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.formula.api import ols

# pheno (QCed)
pheno_file = "/scratch/hx37930/project/admixture/AFR/data/pheno_249metabolites.clean.QCed.AFR.tsv"
pheno_qced = pd.read_csv(pheno_file, sep="\t")

# output directory
outDir = "/scratch/hx37930/project/admixture/AFR/data"

##############################################
############### covariate QC #################
##############################################
# function of ranked-based inverse normalization
def rank_based_inverse_norm(x):
    # Mask the NaN values (they won't be included in the ranking)
    mask = x.isna()
    # Rank the data while keeping NaNs in place
    rank = rankdata(x[~mask], method='average')
    # Create a copy of the original data to preserve NaNs
    inverse_normalized = np.full_like(x, np.nan, dtype=float)
    # Compute the inverse normal transformation
    # (rank - 0.5) / count of non-NaN elements
    inverse_normalized[~mask] = norm.ppf((rank - 0.5) / np.sum(~np.isnan(x)))
    # Return the normalized values
    return inverse_normalized

pheno_int = pheno_qced.copy()
for i in range(16,265):
    pheno_int.iloc[:,i] = rank_based_inverse_norm(pheno_int.iloc[:,i])

#pheno_int.to_csv(f'{outDir}/pheno.QCed.int.AFR.04252025.txt',sep='\t',index=False,header=True,na_rep="na")

# get required covariates: dietary info AHEI_SCORE_m was not included due to high missing rate (only ~ 600 participants with dietary scores)
col_covar_ID_m1 = ["pop2","f.31.0.0","f.21003.0.0","f.22000.0.0","Assessment_centres_11023","Assessment_centres_11001","Assessment_centres_11002","Assessment_centres_11003","Assessment_centres_11004","Assessment_centres_11005","Assessment_centres_11006","Assessment_centres_11007","Assessment_centres_11008","Assessment_centres_11009","Assessment_centres_11010","Assessment_centres_11011","Assessment_centres_11012","Assessment_centres_11013","Assessment_centres_11014","Assessment_centres_11016","Assessment_centres_11017","Assessment_centres_11018","Assessment_centres_11020","Assessment_centres_11021","Assessment_centres_11022"]
col_covar_index_m1 = []
for i in col_covar_ID_m1:
    col_index = pheno_int.columns.get_loc(i)
    col_covar_index_m1.append(col_index)

col_covar_ID_m2 = ["pop2","f.31.0.0","f.21003.0.0","f.22000.0.0", "f.189.0.0","Assessment_centres_11023","Assessment_centres_11001","Assessment_centres_11002","Assessment_centres_11003","Assessment_centres_11004","Assessment_centres_11005","Assessment_centres_11006","Assessment_centres_11007","Assessment_centres_11008","Assessment_centres_11009","Assessment_centres_11010","Assessment_centres_11011","Assessment_centres_11012","Assessment_centres_11013","Assessment_centres_11014","Assessment_centres_11016","Assessment_centres_11017","Assessment_centres_11018","Assessment_centres_11020","Assessment_centres_11021","Assessment_centres_11022"]
col_covar_index_m2 = []
for i in col_covar_ID_m2:
    col_index = pheno_int.columns.get_loc(i)
    col_covar_index_m2.append(col_index)

#################################################################
############### perform association test ########################
#################################################################
def lm_assoc_test(y_index,x_index,toFactor=False):  
    single_meta  = pheno_int.iloc[:,np.append(x_index,y_index)]
    # change data type of specific column to factor/category
    if toFactor == True:
        single_meta['f.20117.0.0'] = single_meta['f.20117.0.0'].replace(-3, np.nan)
        single_meta['f.20117.0.0'] = single_meta['f.20117.0.0'].astype('category')
        single_meta['f.22032.0.0'] = single_meta['f.22032.0.0'].astype('category')
    # rename column name by replacing '.' with '_'
    for c in single_meta.columns:
        single_meta.rename(columns={c: c.replace('.','_')},inplace=True)
    # metabolite
    y_name = single_meta.columns[-1]
    # covariates
    x_names = single_meta.columns[0:-1]     
    # create result table
    result = {'trait':y_name}
    # remove missing value 
    single_meta.dropna(inplace=True)
    # full model
    formula_full_model = y_name + ' ~ ' + ' + '.join(x_names)
    full_model = ols(formula = formula_full_model, data = single_meta).fit()
    sample_size = int(full_model.nobs)
    result['sampleSize'] = sample_size
    # output all covariates except 21 Assessment_centres
    for i in range(len(full_model.params)-21):
        result[f'coef_{full_model.params.index[i]}'] = full_model.params.iloc[i]
        result[f'se_{full_model.params.index[i]}'] = full_model.bse.iloc[i]
        result[f'p_{full_model.pvalues.index[i]}'] = full_model.pvalues.iloc[i]
    # get rsquared of full model
    result['r2_full'] = full_model.rsquared
    # reduce model (leave-one-out covariate)
    for i in x_names[:-21]:
        x_names_reduced = x_names.tolist().copy()
        x_names_reduced.remove(i)
        formula_reduced_model = y_name + ' ~  ' + ' + '.join(x_names_reduced)
        reduced_model = ols(formula = formula_reduced_model, data = single_meta).fit()
        partial_r2 = (full_model.rsquared - reduced_model.rsquared) / (1 - reduced_model.rsquared)
        result[f'r2_{i}'] = partial_r2 
    return pd.DataFrame([result])

# run association test for each metabolite
result_m1 = result_m2 = pd.DataFrame()
for i in range(16,265):
    result_m1 = pd.concat([result_m1, lm_assoc_test(i,col_covar_index_m1,toFactor=False)])
    result_m1.to_csv(f'{outDir}/result_Model1Ancestry.k2.py.txt',sep='\t',index=False,header=True,na_rep="na")
    result_m2 = pd.concat([result_m2, lm_assoc_test(i,col_covar_index_m2,toFactor=False)])
    result_m2.to_csv(f'{outDir}/result_Model2AncestryAddTSI.k2.py.txt',sep='\t',index=False,header=True,na_rep="na")
