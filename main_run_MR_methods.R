source("./Rfuncs/Run_Methods_func.R", echo=TRUE)
source("./Rfuncs/ReadSummaryStat_thisstudy.R", echo=TRUE)
source("./Rfuncs/EstRho_thisstudy.R", echo=TRUE)
source("./Rfuncs/MRConMix.R", echo=TRUE)
source("./Rfuncs/mr_lasso.R", echo=TRUE)

setwd("MR-benchmarking-data/Dataset1")    # Change the directory for different datasets



##### Run MR-APSS ##### 
library(MRAPSS)

IV.Threshold = c(5e-05,5e-06,5e-07,5e-08)

# The exposure and outcome pairs to be analyzed
pairs=read.table("TestedPairs", header=T)

for( i in 1:nrow(pairs)){
  
    exposure = as.character(pairs[i, "exposure"])
    outcome = as.character(pairs[i, "outcome"])
    
    cat("Pair: ", exposure,"~", outcome,"\n")
    
    # read in GWAS summary data for IVs
    clumped = try(read.table(paste0("./MRdata/",exposure,"~",outcome), header = T))
    
    # read in background parameters Omega and C
    C = try(as.matrix(read.table(paste0("./bg_paras/", exposure, "~", outcome, "_C"),
                                 header = F)))
    
    
    Omega = try(as.matrix(read.table(paste0("./bg_paras/", exposure, "~", outcome, "_Omega"),
                                     header = F)))
    
    if(inherits(clumped , 'try-error')) next
    if(inherits(C, 'try-error')) next
    if(inherits(Omega, 'try-error')) next
    if(nrow(clumped) < 3 ) next
    Omega = matrix(Omega[nrow(Omega), ], 2, 2)
    C = matrix(C[nrow(C), ], 2,2)
    
    # The p-value threshold for selection bias correction
    Threshold =  ifelse(IV.Threshold==5e-05, unique(clumped$Threshold), IV.Threshold)
    
    # MR-APSS
    cat("Run MR-APSS ... \n")
    res = run_APSS_func(clumped = clumped,
                        exposure = exposure,
                        outcome = outcome,
                        C = C,
                        Omega=Omega,
                        IV.Threshold = IV.Threshold,
                        Threshold = Threshold,
                        Cor.SelectionBias = T)
    
    # saving resuts
    write.table(res, "MRAPSS.MRres", quote=F, col.names = F, append = T,row.names = F)

}

##### run 11 methods ##### 
set.seed(1234)
methods = c("IVW_fe", "IVW", "RAPS","Egger", "MRMix","Weighted-median", "Weighted-mode",  "MR-Robust","MR-Lasso", "cML-MA", "MR-PRESSO")

# The exposure and outcome pairs to be analyzed
pairs=read.table("TestedPairs", header=T)

for(Thresh in c(5e-08,5e-07,5e-06,5e-05)){

for( i in 1:nrow(pairs)){
  
    exposure = as.character(pairs[i, "exposure"])
    outcome = as.character(pairs[i, "outcome"])
  
      # read in GWAS summary data for IVs
      clumped = try(read.table(paste0("./MRdat/1kgRef_", exposure,"~",outcome), header = T))
      ldmat = NULL
      
      if(inherits(clumped,"try-error")) next
      
      for(i in 1:length(methods)){
        
        res = runMRmethods(dat=clumped, ldmat=ldmat, Threshold=Thresh, 
                           exposure=exposure, outcome=outcome, methods.list = methods[i])
        
        
        write.table(res, "MRmethods.MRres", quote=F, col.names = F, append = T,row.names = F)
 
    }
  }
}



##### Run MR-ConMix  ##### 
set.seed(1234)

# The exposure and outcome pairs to be analyzed
pairs=read.table("TestedPairs", header=T)

for(Thresh in c(5e-08,5e-07,5e-06,5e-05)){

for( i in 1:nrow(pairs)){
  
    exposure = as.character(pairs[i, "exposure"])
    outcome = as.character(pairs[i, "outcome"])
      
      # read in GWAS summary data for IVs
      clumped = try(read.table(paste0("./MRdat/1kgRef_", exposure,"~",outcome), header = T))
      ldmat = NULL
      
      if(inherits(clumped,"try-error")) next
      
      res = runMRmethods(dat=clumped, ldmat=ldmat, Threshold=Thresh, 
                         exposure=exposure, outcome=outcome, ConMixonly = T)
      
      write.table(res, "ConMix.MRres", quote=F, col.names = F, append = T,row.names = F)
      
  }
  
}


##### Run CAUSE  ##### 
library(readr)
library(cause)

start = proc.time()

# The exposure and outcome pairs to be analyzed
pairs=read.table("TestedPairs", header=T)

for(Thresh in c(5e-08,5e-07,5e-06,5e-05)){

for( i in 1:nrow(pairs)){
  
    exposure = as.character(pairs[i, "exposure"])
    outcome = as.character(pairs[i, "outcome"])
    if(exposure==outcome) next
    
    # read GWAS summary statistics
    cat(exposure,"~",outcome,"\n")
    
    fileexp = paste0("./Formatted/exposure.set", as.character(exposure))
    fileout = paste0("./Formatted/outcome.set", as.character(outcome))
    
    X1 = suppressMessages(readr::read_delim(fileexp, delim="\t",
                                            escape_double = FALSE, trim_ws = TRUE, progress = F))
    
    X2 = suppressMessages(readr::read_delim(fileout, delim="\t",
                                            escape_double = FALSE, trim_ws = TRUE, progress = F))
    
    X1$b = X1$Z/sqrt(X1$N)
    X2$b = X2$Z/sqrt(X2$N)
    X1$se = 1/sqrt(X1$N)
    X2$se = 1/sqrt(X2$N)
    
    X <- try(gwas_merge(X1, X2, 
                        snp_name_cols = c("SNP", "SNP"),
                        beta_hat_cols = c("b", "b"),
                        se_cols = c("se", "se"),
                        A1_cols = c("A1", "A1"),
                        A2_cols = c("A2", "A2")))
    
    if(inherits(X , 'try-error')) next
    
    d0 = X1[, c("SNP", "P")]
    colnames(d0) = c("snp", "pval.exp")
    X0 = merge(X, d0, by="snp")
    
    SNPs = try(read.table(paste0("./MRdat/1kgRef_", exposure,"~",outcome), header = T))$SNP
    clumped0 = subset(X0, snp %in% SNPs)
    
    varlist <- with(X, sample(snp, size=min(nrow(X), 1000000), replace=FALSE))
    params <- try(est_cause_params(X, varlist))
    if(inherits(params , 'try-error')) next
    if(inherits(clumped0 , 'try-error')) next
    
    for(Thresh in c(5e-05, 5e-06,5e-07,5e-08)){
      clumped = subset(clumped0, pval.exp <= Thresh)
      # cause
      if(!is.null(clumped)){
        top_ldl_pruned_vars =intersect(as.character(X$snp), as.character(clumped$snp))
        
        cause_res <- try(cause(X=X, variants = top_ldl_pruned_vars , param_ests = params, force=TRUE))
        
        if(inherits( cause_res , 'try-error')) next
        
        res_elpd <- data.frame(exposure,
                               outcome,
                               Thresh,
                               length(top_ldl_pruned_vars),
                               cause_res$elpd)
        
        res.cause.est = summary(cause_res, ci_size=0.95)
        
        res_est = data.frame(exposure, outcome,
                             Thresh,length(top_ldl_pruned_vars),
                             matrix(c(res.cause.est$quants[[2]][,1],
                                      res.cause.est$quants[[2]][,2],
                                      res.cause.est$quants[[2]][,3]), nrow=1))
        
        colnames(res_elpd) = c("exposure","outcome","Threshold","nsnp","model1","model2","delta_elpd", "se_delta_elpd", "Z")
        colnames(res_est) = c("exposure","outcome","Threshold","nsnp", "beta.hat","b_l","b_u","eta","eta_l","eta_u","q","q_l","q_u")
        res_elpd = unique(subset(res_elpd, model1=="sharing"&model2=="causal"))
        res_elpd$pval = pnorm(res_elpd$Z)
        res_est = unique(res_est[, c("exposure","outcome","Threshold","nsnp", "beta.hat","b_l","b_u")])
        res_est$se = (res_est$b_u - res_est$b_l)/2/1.96
        cause_res = unique(merge(unique(res_elpd[, c("exposure","outcome","pval")]),
                                 res_est[, c("exposure","outcome","Threshold","nsnp", "beta.hat","se")],
                                 by=c("exposure","outcome")))
        cause_res$Method = "CAUSE"
        write.table(cause_res, file="CAUSE.MRres", append=T, col.names = F, row.names = F, quote = F)
        
        rm(top_ldl_pruned_vars)
        rm(res_est)
        rm(res_elpd)
        rm(res.cause.est)
        rm(cause_res)
      
    }
    
  }
  
}
print(proc.time()-start)


##### run MR-CUE ##### 
library(MRAPSS)
load("snps.hm3.RData")
set.seed(1234)
library(MR.CUE)
library(ggplot2)
filepan <- vector("list", 22)
NumChr = 22
for(i in 1:NumChr){
  filepan[[i]] <- paste0("MR-benchmarking-data/dataset_mrcue/UK10KCHR", i, "LDhm3.RDS")
}
snpinfo = "MR-benchmarking-data/dataset_mrcue/UK10Ksnpinforhm3.RDS"

ld_r2_thresh = 0.001
lambda = 0.85
pth = 1.96
  
# The exposure and outcome pairs to be analyzed
pairs=read.table("TestedPairs", header=T)

for( i in 1:nrow(pairs)){
  
    exposure = as.character(pairs[i, "exposure"])
    outcome = as.character(pairs[i, "outcome"])
    
    cat(exposure, "~", outcome, " ")
    fileexp = paste0("./Formatted/", as.character(exposure))
    fileout = paste0("./Formatted/", as.character(outcome))
    
    datexp = read.table(fileexp, header=T)
    datout = read.table(fileout, header = T)
    
    paras = est_paras(dat1 = datexp[, !names(datexp) %in% c("chr", "BP")],
                      dat2 = datout[, !names(datout) %in% c("chr", "BP")],
                      trait1.name = exposure,
                      trait2.name = outcome,
                      h2.fix.intercept = F,
                      LDSC = F,
                      ldscore.dir = "dat/eur_w_ld_chr")
    
    dat = merge(paras$dat, snps.hm3[, c("SNP", "chr", "BP")], by="SNP")
    datexp = dat[, c("SNP", "chr", "BP", "A1", "A2", "b.exp", "se.exp", "pval.exp")]
    datout = dat[, c("SNP", "chr", "BP", "A1", "A2", "b.out", "se.out", "pval.out")]
    colnames(datexp) = colnames(datout) = c("SNP", "chr", "BP", "A1", "A2", "beta", "se", "pvalue")
    
    write.table(datexp[, c("SNP", "chr", "BP", "A1", "A2", "beta", "se", "pvalue")], file = paste0(exposure, ".cuedat"), quote = F, col.names = T, row.names = F, sep = "\t")
    write.table(datout[, c("SNP", "chr", "BP", "A1", "A2", "beta", "se", "pvalue")], file = paste0(outcome, ".cuedat"), quote = F, col.names = T, row.names = F, sep = "\t")
    
    
    fileexp = paste0(exposure, ".cuedat")
    fileout = paste0(outcome, ".cuedat")
    
    rho=0
    RhoEst = try(EstRho_thisstudy(fileexp, fileout, filepan, snpinfo, ld_r2_thresh, lambda, pth))
    if(!inherits(RhoEst,"try-error")) {
      rho = mean(RhoEst$Rhores)
    }
    
    for(pva_cutoff  in c(5e-05, 5e-06,5e-07, 5e-08)){
      
      cat(pva_cutoff, "\n")
      data <- try(ReadSummaryStat_thisstudy(fileexp, fileout, filepan, snpinfo, pva_cutoff, lambda))
      if(inherits(data,"try-error")) next
      
      F4gammah <- data$ResF4gammah;
      F4Gammah <- data$ResF4Gammah;
      F4se1 <- data$ResF4se1;
      F4se2 <- data$ResF4se2;
      F4Rblock <- data$ResF4Rblock;
      F4SNPs <- data$ResF4SNPchr
      
      opt = list(agm = 0, bgm = 0, atau1 = 0, btau1 = 0,
                 atau2 = 0, btau2 = 0,
                 a = 2, b = length(F4gammah), maxIter = 5000, thin = 10, burnin = 5000);
      
      RealRes = try(MRCUE(F4gammah, F4Gammah, F4se1, F4se2, F4Rblock, rho=rho, coreNum = 30, opt=opt))
      
      if(!inherits(RealRes,"try-error")) {
        cat("MR-CUE: beta.hat=", RealRes$beta.hat, " se=", RealRes$beta.se, "pval=", RealRes$beta.p.value, "\n")
        
        
        write.table(data.frame(exposure = exposure, outcome = outcome, method = "MR-CUE", Threshold = pva_cutoff,
                               nsnp = length(unlist(F4gammah)), beta = RealRes$beta.hat, se = RealRes$beta.se, pval = RealRes$beta.p.value),
                    file = "MRCUE.MRres",
                    quote=F, col.names = F, append = T,row.names = F)
        
        
      }
    }
    
    unlink(paste0(exposure, "~", outcome, "_exp", ".cuedat"))
    unlink(paste0(exposure, "~", outcome, "_out", ".cuedat"))
}


  ### Run dIVW(lambda=0)
  ### run dIVW
library(mr.divw)
set.seed(1234)
pairs=read.table("TestedPairs", header=T)

for(Thresh in c(5e-08,5e-07,5e-06,5e-05)){

for( i in 1:nrow(pairs)){
  
    exposure = as.character(pairs[i, "exposure"])
    outcome = as.character(pairs[i, "outcome"])
    
    # read in GWAS summary data for IVs
    clumped = try(read.table(paste0("./clumped/UKB600/MRdat/UKB226/1kgRef_", exposure,"~",outcome), header = T))
    if(inherits(clumped,"try-error")) next
    
    for(Thresh in c(5e-08,5e-07,5e-06,5e-05)){
      
      data = subset(clumped, pval.exp<=Thresh)
      res = mr.divw(beta.exposure=data$b.exp, beta.outcome=data$b.out, 
                    se.exposure=data$se.exp, se.outcome=data$se.out, 
                    diagnostics=F, over.dispersion=T)
      
      divw.res = data.frame(exposure = exposure, outcome=outcome,
                            method = "dIVW", Threshold = Thresh,
                            nsnp = res$n.IV , beta = res$beta.hat, se = res$beta.se, 
                            pval = pchisq(res$beta.hat^2/res$beta.se^2, 1, lower.tail = F))
      
      write.table(divw.res, "/import/home/share/xhu/MR-Benchmarking-AJHG-revison-r1/add-dIVW/AtlasUKB_NCO_1kgRef_dIVW.MRres",
                  quote=F, col.names = F, append = T,row.names = F)
      
    }
  }
  
}

