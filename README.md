# Benchmarking Mendelian Randomization methods for causal inference using genome‐wide association study summary statistics
## Introduction

## Datasets
The four datasets used in the MR benchmarking study can be downloaded here:  
[Dataset1](https://gohkust-my.sharepoint.com/:u:/r/personal/maxhu_ust_hk/Documents/MR-benchmarking-data/Dataset1.zip?csf=1&web=1&e=x1eJbM);  
[Dataset2](https://gohkust-my.sharepoint.com/:u:/r/personal/maxhu_ust_hk/Documents/MR-benchmarking-data/Dataset2.zip?csf=1&web=1&e=S3lGxS);  
[Dataset3](https://gohkust-my.sharepoint.com/:u:/r/personal/maxhu_ust_hk/Documents/MR-benchmarking-data/Dataset3.zip?csf=1&web=1&e=7rue4N);  
[Dataset4](https://gohkust-my.sharepoint.com/:u:/r/personal/maxhu_ust_hk/Documents/MR-benchmarking-data/Dataset4.zip?csf=1&web=1&e=MOQUXH). 

In each dataset, the GWAS summary data files for each exposure or outcome trait (after data preprocessing) and the summary statistics after IV selection (step 1) are available.

## R code
### Step 0: install required packages
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

### Step 1: IV selection 
In this step, we perform IV selection for each trait pair in each dataset.
The R code for IV selection is available in [IV_selection.R](https://github.com/YangLabHKUST/MRbenchmarking/blob/main/IV_selection.R).


### Step 2: run MR Methods
In this step, we applied each compared method to each dataset.
The R codes for running the 15 MR methods for each dataset are available in [main_run_MR_methods.R](https://github.com/YangLabHKUST/MRbenchmarking/blob/main/main_run_MR_methods.R).
To run the codes of *main_run_MR_methods.R*, you must load the required packages and the R functions in the folder [Rfuncs](https://github.com/YangLabHKUST/MRbenchmarking/tree/main/Rfuncs). 

## Results

## Reference
Benchmarking Mendelian Randomization methods for causal inference using genome‐wide association study summary statistics

## Contact information
Please feel free to contact Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
