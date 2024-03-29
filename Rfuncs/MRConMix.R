mr_conmix <- function(snps, Bx ,By, Bxse, Byse, psi    = 0, CIMin  = NA, CIMax  = NA, CIStep = 0.01,
                      alpha = 0.05)
{
  
  nsnps = length(Bx)
  
  ratio = By/Bx; ratio.se = abs(Byse/Bx);
  if (is.na(CIMin)) { CIMin = min((By-2*Byse)/Bx) }
  if (is.na(CIMax)) { CIMax = max((By+2*Byse)/Bx) }
  if (psi < 0 | psi == 0) {   psi = 1.5*sd(ratio)  }
  theta = seq(from = CIMin, to = CIMax, by = CIStep)
  iters = length(theta)
  lik=NULL
  for (j1 in 1:iters) {
    lik.inc = -(theta[j1]-ratio)^2/2/ratio.se^2 - log(sqrt(2*pi*ratio.se^2))
    lik.exc = -ratio^2/2/(psi^2+ratio.se^2) - log(sqrt(2*pi*(psi^2+ratio.se^2)))
    valid = (lik.inc>lik.exc)*1
    lik[j1] = sum(c(lik.inc[valid==1], lik.exc[valid==0]))
    if (which.max(lik)==length(lik)) { valid.best = valid }
  }
  phi = ifelse(sum(valid.best)<1.5, 1,
               max(sqrt(sum(((ratio[valid.best==1]-weighted.mean(ratio[valid.best==1],
                                                                 ratio.se[valid.best==1]^-2))^2*
                               ratio.se[valid.best==1]^-2))/(sum(valid.best)-1)), 1))
  loglik = lik
  
  lik.inc0 = -ratio^2/2/ratio.se^2 - log(sqrt(2*pi*ratio.se^2))
  lik.exc0 = -ratio^2/2/(psi^2+ratio.se^2) - log(sqrt(2*pi*(psi^2+ratio.se^2)))
  valid = (lik.inc0>lik.exc0)*1
  loglik0 = sum(c(lik.inc0[valid==1], lik.exc0[valid==0]))
  
  
  whichin = which(2*loglik>(2*max(loglik)-qchisq(1-alpha, df=1)*phi^2))
  # provides an index of estimate values in the 95% confidence interval
  betaConMix = CIMin+CIStep*(which.max(loglik)-1)
  # modal estimate
  CIRange    = CIMin+CIStep*(whichin-1);
  CILower <- c(min(CIRange), CIRange[which(diff(CIRange)>1.01*CIStep)+1])
  CIUpper <- c(CIRange[which(diff(CIRange)>1.01*CIStep)], max(CIRange))
  
  Pvalue = pchisq(2*(max(loglik)-loglik0)*phi^2, df=1, lower.tail=FALSE)
  
  return(list(Psi = as.numeric(psi),
              Estimate = as.numeric(betaConMix),
              CIRange  = as.numeric(CIRange),
              CILower  = as.numeric(CILower),
              CIUpper  = as.numeric(CIUpper),
              CL = min(CIRange),
              CU = max(CIRange),
              CIMin    = as.numeric(CIMin),
              CIMax    = as.numeric(CIMax),
              CIStep   = as.numeric(CIStep),
              Valid    = as.numeric(which(valid.best==1)),
              ValidSNPs= as.character(snps[which(valid.best==1)]),
              Pvalue   = as.numeric(Pvalue),
              SNPs = nsnps,
              Alpha = alpha))
}

