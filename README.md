# Benchmarking Mendelian Randomization methods for causal inference using genome‐wide association study summary statistics
## The experimental design for benchmarking MR methods
We present a benchmarking analysis of MR methods for causal inference with real-world genetic datasets. Our focus is on MR methods that utilize GWAS summary statistics as input, as they do not require access to individual-level GWAS data and are widely applicable. Specifically, we consider 16 MR methods, including the standard IVW (fixed), IVW (random) and 14 other advanced MR methods: dIVW, Egger, RAPS,  Weighted-median, Weighted-mode, MR-PRESSO, MRMix, cML-MA, MR-Robust, MR-Lasso, MR-CUE, CAUSE, MRAPSS and MR-ConMix (Figure A). The procedure for running the MR methods is outlined in Figure B.  To assess the performance of these MR methods, we utilized real-world datasets and focused on three key aspects: type I error control, the accuracy of causal effect estimates, replicability, and power (Figure C).
![My Image](design.png)

## Datasets
The datasets used in the MR benchmarking study can be downloaded here.

### Dataset 1: GWASATLAS Dataset for evaluation of type I error control in confounding scenario (a): Population stratification
[GWASs for exposures](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EblwNZLAorRAvCLbvYugudEBaGtWc72q2HrhAvGkCHTmaA?e=DvW21W);
[GWASs for outcomes](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EQdc-MQQeLZKn7G2Oj9N30sBuuZCZ4c3LCjgi92NokmAOw?e=ei6ubC);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset-GWASATLAS-negativecontrol.zip?download=1);

### Dataset 2: the Neal Lab Dataset for evaluation of type I error control in confounding scenario (a): Population stratification  
[GWASs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EU1e5jC9jS9DptDaVKUJwlsB6BgCAQ8OWntCWGzP7zWzdA?e=r6B2D2);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset-NealeLab-negativecontrol.zip?download=1).

### Dataset 3: the Pan UKBB Dataset for evaluation of type I error control in confounding scenario (a): Population stratification  
[GWASs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/ERu0_x-u0FhDsjbFHhmB1BsBPzTnc6VvPBclwnU2QeEs5g?e=NyffU2);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset-PanUKBB-negativecontrol.zip?download=1).

### Dataset 4: the dataset for evaluation of type I error control in confounding scenario (b): Pleiotropy  
[GWASs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EVGtR-AH6WBCvmleRgAmZJIBYDK8tty61YxbeFobnMRCRg?e=6nL2d0);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset-Pleiotropy-negativecontrol.zip?download=1)

### Dataset 5: the dataset for evaluation of type I error control in confounding scenario (c): Family-level confounders  
[GWASs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/Efflau1WW7VAhgyBEaZsw2IBh59CUv7HLdbCE-cAPJuesw?e=N9uflB);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset-familylevelconf-negativecontrol.zip?download=1)

### Dataset 6: the dataset for evaluation of the accuracy of causal effect estimates  
[GWASs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EVGtR-AH6WBCvmleRgAmZJIBYDK8tty61YxbeFobnMRCRg?e=6nL2d0);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset_ukb-ukb.zip?download=1);

### Dataset 7: the dataset for evaluation of replicability  
[GWASs](https://gohkust-my.sharepoint.com/:u:/g/personal/maxhu_ust_hk/EcfHZhJfqrxLiBiIV8W5BWgBJgIBklOJcc0ebggGqCD4wg?e=iLuN8l);
[Formatted IV data for MR analysis](https://zenodo.org/records/10929572/files/dataset-LDL-CAD.zip?download=1);

Note: 
(1) "GWASs" refers to the formatted GWAS summary-level data files after quality control;   
(2) "Formatted IV data for MR analysis" contains the following three types of files:    
     "Tested Trait pairs": the exposure-outcome trait pairs to be analyzed;    
     "MRdat": refers to the summary statistics of LD clumped IV sets for each trait pair tested which can be directed used for MR analysis;   
     "bg_paras": refers to the estimated background parameters "Omega" and "C" which will be used for MR estimation in MR-APSS.    
(3) The details on data preprocessing including quality control of GWAS summary statistics, formatting GWASs, and  LD clumping for IV selection can be found in the supplementary note of our paper[1].  
    Implementation details on data preprocessing can be found in the [MR-APSS software tutorial]((https://github.com/YangLabHKUST/MR-APSS/blob/master/MRAPSS_Rpackage_Tutorial.pdf)) on MR-APSS [GitHub](https://github.com/YangLabHKUST/MR-APSS) website.


## R code
### Install required packages
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

### Run MR Methods
We perform IV selection for each trait pair in each dataset.
The R code for IV selection is available in [IV_selection.R](https://github.com/YangLabHKUST/MRbenchmarking/blob/main/IV_selection.R).

We then applied each compared method using the dataset after IV selection.
The R codes for running the 15 MR methods for each dataset are available in [main_run_MR_methods.R](https://github.com/YangLabHKUST/MRbenchmarking/blob/main/main_run_MR_methods.R).
To run the codes of *main_run_MR_methods.R*, you must load the required packages and the R functions in the folder [Rfuncs](https://github.com/YangLabHKUST/MRbenchmarking/tree/main/Rfuncs). 

## Results of MR methods
[Results for dataset 1](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EU2TpRPP8wdHk8Fd0UhxmaUBOHfbxJvFJ9yoEcTdPIJmKw?e=834Hog);  
[Results for dataset 2](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EcZooOkcIABOj72yHxgns2kB9AOxH63CHkxzCSMSFJU1Fg?e=fz1vnq);  
[Results for dataset 3](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EcF4WbLwCtxMv8tV1VQtaq8BkTCh0lJND7BmPGjgKWfozA?e=G91CeG);  
[Results for dataset 4](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EYpIpzxmn8tMixfwdNFpDQYB51IRf3UzA9b_y3v_2ZO6MQ?e=J6ZpM9);  
[Results for dataset 5](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EacRr7UW-l1KlSAsRVXPenkBHbRiZM2bSFnfukEEbBMggg?e=0xRcTf);  
[Results for dataset 6](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EZK18RG4tlxPo1L-eEaYyMwBhHfEEw69PDskAQHH9lcwqw?e=BLF3j2);  
[Results for dataset 7](https://gohkust-my.sharepoint.com/:x:/g/personal/maxhu_ust_hk/EeTc1DeCfhpKiT72cH9FXf4Buxqfo5291qkjQDw4sr8hqw?e=ZQKGdm). 

## Reference
Xianghong Hu, Mingxuan Cai, Jiashun Xiao, Xiaomeng Wan, Zhiwei Wang, Hongyu Zhao, and Can Yang, Benchmarking Mendelian Randomization methods for causal inference using genome‐wide association study summary statistics, https://medrxiv.org/cgi/content/short/2024.01.03.24300765v1.

## Contact information
Please feel free to contact Xianghong Hu (maxhu@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any questions.
