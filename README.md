# Benchmarking Mendelian Randomization methods for causal inference using genome‐wide association study summary statistics

## Install required packages
#install.packages("devtools")
#install.packages("remotes")

devtools::install_github("gqi/MRMix")

devtools::install_github("xue-hr/MRcML")

devtools::install_github("jean997/cause@v1.2.0")

devtools::install_github("rondolab/MR-PRESSO")

install.packages("MendelianRandomization")

devtools::install_github("YangLabHKUST/MR-APSS")

devtools::install_github("QingCheng0218/MR.CUE@main")

remotes::install_github("MRCIEU/TwoSampleMR")

devtools::install_github("qingyuanzhao/mr.raps")

install.packages(“robustbase”)

## R code
### Step 1：data-preprocessing
The R codes for data preprocessing for each dataset are available [].

### Step 2: run MR Methods
The R codes for running the 15 MR methods for each dataset are available in [main_run_MR_methods.R]().
To run the codes of *main_run_MR_methods.R*, you must load the required packages and the R functions in the folder [Rfuncs]. 

## Datasets
The four datasets used in the MR benchmarking study can be downloaded here:  
[Dataset1]();  
[Dataset2]();  
[Dataset3]();  
[Dataset4](). 

The GWAS summary data files (after data preprocessing) and the summary statistics for selected IVs are available for downloading.

## Results

## Reference
Benchmarking Mendelian Randomization methods for causal inference using genome‐wide association study summary statistics

## Contact information
Please feel free to contact Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
