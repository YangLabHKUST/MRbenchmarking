
runMRmethods <- function(dat, ldmat=NULL, Threshold=5e-08, exposure="exposure", outcome="outcome",
                         methods.list = c("IVW","dIVW", "RAPS","Egger", "MRMix",
                                          "Weighted-median", "Weighted-mode", "BWMR","MR-PRESSO","MR-Lasso"),
                         ConMixonly=FALSE){
  set.seed(1234)
  cat("Pair: ", exposure,"~", outcome,"\n")
  dat = subset(dat, b.exp!=0)
  
  indx = which(dat$pval.exp <= Threshold)
  
  nsnp = length(indx)
  
  res= NULL
  
  
  if(nsnp>3){
    
    # IVW-fixed
    
    if("IVW_fe" %in% methods.list & ConMixonly==FALSE){
      cat("IVW \n")
      res.IVW <- try(TwoSampleMR::mr_ivw_fe(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      
      if(!inherits(res.IVW,"try-error")){
        cat("IVW: beta.hat=", res.IVW$b, " se=", res.IVW$se, "pval=", res.IVW$pval, "\n")
        
        IVW_res <- data.frame(exposure = exposure, outcome = outcome, method =  "IVW_fe", Threshold = Threshold,
                              nsnp = nsnp, beta =  res.IVW$b, se = res.IVW$se, pval = res.IVW$pval)
        res = rbind(res, IVW_res)
      }
    }
    
    
    # IVW-random
    if("IVW" %in% methods.list & ConMixonly==FALSE){
      cat("IVW \n")
      res.IVW <- try(TwoSampleMR::mr_ivw(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      
      if(!inherits(res.IVW,"try-error")){
        cat("IVW: beta.hat=", res.IVW$b, " se=", res.IVW$se, "pval=", res.IVW$pval, "\n")
        
        IVW_res <- data.frame(exposure = exposure, outcome = outcome, method =  "IVW", Threshold = Threshold,
                              nsnp = nsnp, beta =  res.IVW$b, se = res.IVW$se, pval = res.IVW$pval)
        res = rbind(res, IVW_res)
      }
    }
    
 
    # RAPS
    if("RAPS" %in% methods.list & ConMixonly==FALSE){
      raps_res = NULL
      raps.df = data.frame(beta.exposure = dat$b.exp[indx], beta.outcome = dat$b.out[indx],
                           se.exposure = dat$se.exp[indx], se.outcome = dat$se.out[indx])
      
      res.raps <- try(suppressWarnings(mr.raps::mr.raps(raps.df, diagnostics=F)))
      
      if(!inherits(res.raps,"try-error")){
        
        cat("RAPS: beta.hat=", res.raps$beta.hat, " se=", res.raps$beta.se, "\n")
        
        raps_res <- data.frame(exposure = exposure, outcome = outcome, method =  "RAPS", Threshold = Threshold,
                               nsnp = nsnp, beta =  res.raps$beta.hat, se = res.raps$beta.se,
                               pval = pchisq(res.raps$beta.hat^2/res.raps$beta.se^2, 1, lower.tail = F))
      }
      res = rbind(res, raps_res)
      
    }
    
    
    # Egger
    if("Egger" %in% methods.list & ConMixonly==FALSE){
      
      res.egger <- try(TwoSampleMR::mr_egger_regression(dat$b.exp[indx], dat$b.out[indx],
                                                        dat$se.exp[indx], dat$se.out[indx]))
      if(!inherits(res.egger,"error")){
        
        cat("Egger: beta.hat=", res.egger$b, " se=", res.egger$se, "pval=", res.egger$pval, "\n")
        
        egger_res <- data.frame(exposure = exposure, outcome = outcome, method =  "Egger", Threshold = Threshold,
                                nsnp = nsnp, beta =  res.egger$b, se = res.egger$se, pval = res.egger$pval)
        res = rbind(res, egger_res)
      }
      
    }
    
    # MRmix 
    if("MRMix" %in% methods.list & ConMixonly==FALSE){
      
      res.MRMix <- try(MRMix::MRMix(dat[indx,]$b.exp, dat[indx,]$b.out, dat[indx,]$se.exp, dat[indx,]$se.out))
      
      if(!inherits(res.MRMix,"try-error")){
        
        
        cat("MRMix: beta.hat=", res.MRMix$theta, " se=", res.MRMix$SE_theta, "pval=", res.MRMix$pvalue_theta, "\n")
        
        
        MRMix_res <- data.frame(exposure = exposure, outcome = outcome, method = "MRMix", Threshold = Threshold,
                                nsnp = nsnp, beta = res.MRMix$theta, se = res.MRMix$SE_theta, pval = res.MRMix$pvalue_theta)
        
        res = rbind(res, MRMix_res)
        
      }
    }
    
    # Weighted-median
    if("Weighted-median" %in% methods.list & ConMixonly==FALSE){
      
      res.median <- try(TwoSampleMR::mr_weighted_median(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      cat("Weighted-median: beta.hat=", res.median$b, " se=", res.median$se, "pval=", res.median$pval, "\n")
      
      if(!inherits(res.median,"try-error")) {
        median_res <- data.frame(exposure = exposure, outcome = outcome, method =  "Weighted-median", Threshold = Threshold,
                                 nsnp = nsnp, beta = res.median$b, se = res.median$se, pval = res.median$pval)
        res = rbind(res, median_res)
      }
      
    }
    
    # Weighted-mode
    if("Weighted-mode" %in% methods.list & ConMixonly==FALSE){
      
      res.mode <- try(TwoSampleMR::mr_weighted_mode(dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx]))
      
      cat("Weighted-mode: beta.hat=", res.mode$b, " se=", res.mode$se, "pval=", res.mode$pval, "\n")
      
      if(!inherits(res.mode,"try-error")){
        mode_res <- data.frame(exposure = exposure, outcome = outcome, method =  "Weighted-mode", Threshold = Threshold,
                               nsnp = nsnp, beta = res.mode$b, se = res.mode$se, pval = res.mode$pval)
        res = rbind(res, mode_res)
      }
      
      
    }
    

    # MR-PRESSO
    if("MR-PRESSO" %in% methods.list & ConMixonly==FALSE){
      
      NbDistribution = 1000
      if(nrow(dat[indx,]) >= NbDistribution ) NbDistribution = 1.1 * nrow(dat[indx,])
      res.presso <- try(MRPRESSO::mr_presso(BetaOutcome = "b.out", BetaExposure = "b.exp",
                                            SdOutcome = "se.out", SdExposure = "se.exp",
                                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                                            data = dat[indx,], NbDistribution = NbDistribution,  
                                            SignifThreshold = 0.05, seed=1234))
      save(res.presso, file="CTSD.PRESSO")
      if(!inherits(res.presso, "try-error")){
        
        if(is.na(res.presso[[1]][2,3])){
          presso_res <- data.frame(exposure = exposure, outcome = outcome, method =  "MR-PRESSO", Threshold = Threshold,
                                   nsnp = nsnp, beta = res.presso[[1]][1,3], se = res.presso[[1]][1,4], pval = res.presso[[1]][1,6]) 
          cat("MR-PRESSO: beta.hat=", res.presso[[1]][1,3], " se=", res.presso[[1]][1,4], "pval=", res.presso[[1]][1,6], "\n")
          
        }else{
          presso_res <- data.frame(exposure = exposure, outcome = outcome, method =  "MR-PRESSO", Threshold = Threshold,
                                   nsnp = nsnp, beta = res.presso[[1]][2,3], se = res.presso[[1]][2,4], pval = res.presso[[1]][2,6]) 
          cat("MR-PRESSO: beta.hat=", res.presso[[1]][2,3], " se=", res.presso[[1]][2,4], "pval=", res.presso[[1]][2,6], "\n")
          
        }
        
        res = rbind(res, presso_res)
      }
    }
    
    # MR-Robust
    if("MR-Robust"  %in% methods.list & ConMixonly==FALSE){
      
      rob =  try(robustbase::lmrob(dat$b.out[indx]~dat$b.exp[indx]-1, 
                                   weights=dat$se.out[indx]^-2, k.max=500))
      
      if(!inherits(rob,"try-error")) {
        fitrob = summary(rob)
        betaIVW.robust = fitrob$coef[1]
        # sebetaIVW.robust.fixed = summary(lmrob(betaYG~betaXG-1, weights=sebetaYG^-2, k.max=500))$coef[1,2]/
        #   summary(lmrob(betaYG~betaXG-1, weights=sebetaYG^-2, k.max=500))$sigma
        sebetaIVW.robust.random = fitrob$coef[1,2]/min(fitrob$sigma,1)
        IVW.robust.random_res <- data.frame(exposure = exposure, outcome = outcome, 
                                            method =  "MR-Robust", Threshold = Threshold,
                                            nsnp = nsnp, beta =  betaIVW.robust, se = sebetaIVW.robust.random,
                                            pval = pchisq(betaIVW.robust^2/sebetaIVW.robust.random^2, 1, lower.tail = F))
        res = rbind(res, IVW.robust.random_res)
        
      }
    }
    
    # MR-Lasso
    if("MR-Lasso"  %in% methods.list & ConMixonly==FALSE){
      
      Lasso_fit = try(mr_lasso(snps=dat$SNP[indx],dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx], lambda=NULL))
      
      if(!inherits(Lasso_fit,"try-error")) {
        cat("MR-Lasso: beta.hat=", Lasso_fit$Estimate, " pval=", Lasso_fit$Pvalue,"\n")
        
        Lasso_res = data.frame(exposure = exposure, outcome = outcome, method = "MR-Lasso",Threshold = Threshold,
                               nsnp = nsnp, beta =  Lasso_fit$Estimate, se = Lasso_fit$StdError, pval = Lasso_fit$Pvalue)
        
        res = rbind(res,  Lasso_res)
      }
    }
    
    
    # cML-MA-DP
    if("cML-MA-DP" %in% methods.list& ConMixonly==FALSE){
      final_dat = dat[indx,]
      n= min(median(1/final_dat$se.exp^2), median(1/final_dat$se.out^2))
      res_cML = try(suppressWarnings(MRcML::mr_cML_DP(final_dat$b.exp,
                                                      final_dat$b.out,
                                                      final_dat$se.exp,
                                                      final_dat$se.out,
                                                      n = n,
                                                      random_start = 10,
                                                      random_start_pert = 10,
                                                      random_seed = 1,
                                                      num_pert = 200)))
      cat("cML-MA: beta.hat=", res_cML$MA_BIC_DP_theta, " se=", res_cML$MA_BIC_DP_se, "pval=", res_cML$MA_BIC_DP_p, "\n")
      
      
      cML_res = data.frame(exposure = exposure, outcome = outcome, method =  "cML-MA", Threshold = Threshold,
                           nsnp = nsnp, beta =  res_cML$MA_BIC_DP_theta, se = res_cML$MA_BIC_DP_se, pval = res_cML$MA_BIC_DP_p)
      
      res = rbind(res, cML_res)
    }   
    
    # cML-MA
    if("cML-MA" %in% methods.list & ConMixonly==FALSE){
      final_dat = dat[indx,]
      n= min(median(1/final_dat$se.exp^2), median(1/final_dat$se.out^2))
      res_cML = try(suppressWarnings(MRcML::mr_cML(final_dat$b.exp,
                                                   final_dat$b.out,
                                                   final_dat$se.exp,
                                                   final_dat$se.out,
                                                   n = n,
                                                   random_start = 10,
                                                   random_seed = 1)))
      cat("cML-MA: beta.hat=", res_cML$MA_BIC_theta, " se=", res_cML$MA_BIC_se, "pval=", res_cML$MA_BIC_p, "\n")
      
      
      cML_res = data.frame(exposure = exposure, outcome = outcome, method =  "cML-MA", Threshold = Threshold,
                           nsnp = nsnp, beta =  res_cML$MA_BIC_theta, se = res_cML$MA_BIC_se, pval = res_cML$MA_BIC_p)
      
      res = rbind(res, cML_res)
    }   
    
    
    # MR-ConMix
    if(ConMixonly){
      
      ConMix_fit = try(mr_conmix(snps=dat$SNP[indx],dat$b.exp[indx], dat$b.out[indx], dat$se.exp[indx], dat$se.out[indx], CIMin  = -1.5, CIMax  = 1.5))
      
      if(!inherits(ConMix_fit,"try-error")) {
        cat("MR-ConMix: beta.hat=", ConMix_fit$Estimate, " pval=", ConMix_fit$Pvalue," CL=", ConMix_fit$CL, " CU=", ConMix_fit$CU,"\n")
        
        res = data.frame(exposure = exposure, outcome = outcome, method =  "MR-ConMix", Threshold = Threshold,
                         nsnp = nsnp, beta =  ConMix_fit$Estimate, pval = ConMix_fit$Pvalue, CL=ConMix_fit$CL, CU=ConMix_fit$CU)
      }
    }
  }
  
  return(res)
}


### Func for running MR-APSS 
run_APSS_func <- function(clumped=NULL,
                          exposure=NULL,
                          outcome=NULL,
                          C = diag(2),
                          Omega = matrix(0, 2, 2),
                          IV.Threshold = 5e-08, # IV selection threshold
                          Threshold = 5e-08,    # threshold for correcting for selection bias
                          Cor.SelectionBias=T){
  
  
  res=NULL
  
  for(i in 1:length(Threshold)){
    
    start = proc.time()
    if(!is.null(clumped)){
      
      test = subset(clumped, pval.exp <= IV.Threshold[i])
      if(nrow(test) < 4 ) next
      test$Threshold = Threshold[i]
      cat("IV selection threshold:", IV.Threshold[i] ,"\n")
      
      
      MRres = try(MRAPSS::MRAPSS(test,
                                 exposure=exposure,
                                 outcome=outcome,
                                 C= C,
                                 Omega= Omega,
                                 Cor.SelectionBias = Cor.SelectionBias))
      
      print(proc.time()-start)
      
      if(inherits(MRres, 'try-error')) {
        MRres=NULL
      }
    }
    
    
    res0 = data.frame(exposure = MRres$exposure,
                      outcome = MRres$outcome,
                      nSNP = nrow(MRres$MRdat),
                      pi0 = MRres$pi0,
                      nvalid = nrow(MRres$MRdat)*MRres$pi0,
                      sigma.sq= MRres$sigma.sq,
                      tau.sq= MRres$tau.sq,
                      beta = MRres$beta,
                      se = MRres$beta.se,
                      pval= MRres$pvalue,
                      method = MRres$method,
                      Threshold = Threshold[i],
                      IVstrength = MRres$IVsignal.sum
    )
    
    if(nrow(res0)!=0){
      res0$IV.Threshold = IV.Threshold[i]
    }
    
    res = rbind(res, res0)
    
  }
  
  return(res)
  
}




