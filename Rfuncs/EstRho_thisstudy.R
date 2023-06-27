EstRho_thisstudy <- function(fileexp, fileout, filepan, snpinfo, ld_r2_thresh, lambad, pth) {
  
  daexp = read.table(fileexp, header = T)
  daout = read.table(fileout, header = T)
  
  daexp$a1 = 0
  daexp$a2 = 0
  daexp$a1[daexp$A1=="C"] = 1
  daexp$a1[daexp$A1=="G"] = 1
  daexp$a2[daexp$A2=="C"] = 1
  daexp$a2[daexp$A2=="G"] = 1
  daexp = daexp[daexp$a1+daexp$a2==1,]
  daexp[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daexp$chr = as.numeric(daexp$chr)
  daexp = daexp[!is.na(daexp$chr),]
  daexp = unique(daexp)
  rownames(daexp) = daexp$SNP
  
  daout = daout[daout$SNP %in% daexp$SNP, ]
  daout = unique(daout)
  rownames(daout) = daout$SNP
  daout$a1 = 0
  daout$a2 = 0
  daout$a1[daout$A1=="C"] = 1
  daout$a1[daout$A1=="G"] = 1
  daout$a2[daout$A2=="C"] = 1
  daout$a2[daout$A2=="G"] = 1
  daout = daout[daout$a1+daout$a2==1,]
  bp = daout$BP;
  daout[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daout$chr = as.numeric(daout$chr)
  daout = daout[!is.na(daout$chr),]
  
  daexp$SNP = as.character(daexp$SNP)
  daout$SNP = as.character(daout$SNP)
  #daexp = daexp[daout$SNP,]
  #daout = daout[daout$a1==daexp$a1|daout$a1==daexp$a2,]
  #daexp = daexp[daout$SNP,]
  # ===================================================== #
  # Remove MHC region.
  QCres = RemoveMHCCfun(mhcstart, mhcend, daexp$beta, daout$beta, daexp$se, daout$se, bp, daout$chr,
                        daout$SNP,  Inf, Inf)
  rsnameRMHC = QCres$rsnamenew;
  
  if(length(rsnameRMHC)!=dim(daout)[1]){
    daexp = daexp[daexp$SNP%in%rsnameRMHC, ];
    daout = daout[daout$SNP%in%rsnameRMHC, ];
  }
  
  
  # ===================================================== #
  daRefsnpinfo = readRDS(snpinfo);
  daRefsnpinfosame = daRefsnpinfo[daRefsnpinfo$SNP%in%daexp$SNP, ];
  daRefsnpinfosame$a1 = 0
  daRefsnpinfosame$a1[daRefsnpinfosame$A1=="C"] = 1
  daRefsnpinfosame$a1[daRefsnpinfosame$A1=="G"] = 1
  daRefsnpinfosame[, "A1"] = NULL
  
  # ===================================================== #
  # Random select SNPs in case that the number of SNPs are too large.
  if(dim(daRefsnpinfosame)[1] > 1e6){
    idxrand = sort(sample(1:dim(daRefsnpinfosame)[1], 1000000, replace = FALSE));
    daRefsnpinfosame = daRefsnpinfosame[idxrand, ];
    daexp = daexp[daexp$SNP%in%daRefsnpinfosame$SNP, ];
    daout = daexp[daout$SNP%in%daRefsnpinfosame$SNP, ];
  }
  
  # ===================================================== #
  ResF4gamma   = vector("list", 0)
  ResF4gammase = vector("list", 0)
  ResF4Gamma   = vector("list", 0)
  ResF4Gammase = vector("list", 0)
  # ResF4SNPchr  = vector("list", 0)
  # ResF4R       = vector("list", 0)
  # Res_id  = c()
  
  for (chrk in 1:22) {
    if (is.null(filepan[[chrk]])) next
    if (!file.exists(filepan[[chrk]])) next
    dapan = readRDS(filepan[[chrk]])
    dapan[] <- lapply(dapan, function(x) if (is.factor(x)) as.character(x) else {x})
    
    daexp_k = daexp[daexp$chr==chrk,]
    daout_k = daout[daexp_k$SNP,]
    
    # ===================================================== #
    RefSNP_k = daRefsnpinfosame[daRefsnpinfosame$CHR==chrk, ]
    snpint = intersect(union(dapan$SNP1, dapan$SNP2), daout_k$SNP);
    stopifnot(sum(RefSNP_k$SNP%in%snpint)==length(snpint))
    # ===================================================== #
    
    dapan = dapan[dapan$SNP1 %in% daout_k$SNP,]
    dapan = dapan[dapan$SNP2 %in% daout_k$SNP,]
    
    snp = union(dapan$SNP1, dapan$SNP2)
    daexp_k = daexp_k[daexp_k$SNP %in% snp, ]
    daout_k = daout_k[daexp_k$SNP,]
    
    snppan = union(dapan$SNP1, dapan$SNP2)
    if(length(snppan)!=length(snpint)){
      sigleSNP = snpint[!snpint%in%snppan];
      Sigpan = data.frame(rep(chrk, length(sigleSNP)),paste0(chrk*10, max(dapan$BlockID)+seq(1:length(sigleSNP))), sigleSNP, sigleSNP, rep(1, length(sigleSNP)))
      colnames(Sigpan) = colnames(dapan)
      dapan = rbind(dapan, Sigpan)
    }
    
    
    arma_blc_id = unique(dapan$BlockID)
    
    nid = length(arma_blc_id)
    F4gamma   = vector("list", nid)
    F4gammase = vector("list", nid)
    F4Gamma   = vector("list", nid)
    F4Gammase = vector("list", nid)
    # F4SNPchr  = vector("list", nid)
    # F4R       = vector("list", nid)
    snpset = rep(FALSE, nid)
    for (k in 1:nid) {
      uid = arma_blc_id[k]
      da2 = dapan[dapan$BlockID==uid,]
      comsnp = union(da2$SNP1, da2$SNP2)
      
      daRef = RefSNP_k[RefSNP_k$SNP%in%comsnp, ];
      da0 = daexp_k[comsnp, ]; 
      da1 = daout_k[comsnp, ];
      # ===================================================== #
      # match Ref panel data.
      da0$beta[da0$a1!=daRef$a1] = -1*da0$beta[da0$a1!=daRef$a1];
      da1$beta[da0$a1!=daRef$a1] = -1*da1$beta[da0$a1!=daRef$a1];
      # ===================================================== #
      
      M = diag(rep(1,length(comsnp)))
      rownames(M) = comsnp
      colnames(M) = comsnp
      for (j in 1:nrow(da2)) {
        M[da2$SNP1[j], da2$SNP2[j]] = da2$r[j]*lambad
        M[da2$SNP2[j], da2$SNP1[j]] = da2$r[j]*lambad
      }
      rownames(M) = NULL
      colnames(M) = NULL
      
      id = LDclump(M*1.0, ld_r2_thresh) + 1
      
      F4gamma[[k]] = da0$beta[id]
      F4gammase[[k]] = da0$se[id]
      F4Gamma[[k]] = da1$beta[id]
      F4Gammase[[k]] = da1$se[id]
      # F4SNPchr[[k]] = comsnp[id]
      # F4R[[k]] = M[id,id]
    }
    
    ResF4gamma    = append(ResF4gamma, F4gamma)
    ResF4gammase  = append(ResF4gammase, F4gammase)
    ResF4Gamma    = append(ResF4Gamma, F4Gamma)
    ResF4Gammase = append(ResF4Gammase, F4Gammase)
    # ResF4SNPchr   = append(ResF4SNPchr, F4SNPchr)
    # ResF4R        = append(ResF4R, F4R)
    # Res_id        = c(Res_id, arma_blc_id)
  }
  
    try(
    rm(daexp, daout, dapan, daexp_k, daout_k, da0, da1, da2)
  )
  
  bh1_ind = unlist(ResF4gamma)
  bh2_ind = unlist(ResF4Gamma)
  se1_ind = unlist(ResF4gammase)
  se2_ind = unlist(ResF4Gammase)
  
  z1_ind = bh1_ind / se1_ind;
  z2_ind = bh2_ind / se2_ind;
  
  maxIter = 4000;
  thin = 10;
  burnin = 1000;
  
  nsave = maxIter / thin;
  
  if(length(pth)==1){
    Rhores = rep(NA, nsave);
    a = rep(-pth, 2);
    b = rep(pth, 2);
    z1_new = z1_ind [which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    z2_new = z2_ind[which(abs(z1_ind) < pth&abs(z2_ind) < pth)];
    rhores = truncEstfun(a, b, z1_new, z2_new, maxIter, burnin, thin)
    rhohat = mean(rhores);
    p1 = length(z1_new);
    pvalue = testR(rhohat, p1);
    Rhores = rhores;
    pres = p1;
  }else{
    rhohat = pvalue = rep(NA, length(pth));
    Rhores = matrix(NA, nrow = nsave, ncol = length(pth));
    pres = rep(NA, length(pth));
    for(i in 1:length(pth)){
      pth1 = pth[i];
      a = rep(-pth1, 2);
      b = rep(pth1, 2);
      z1_new = z1_ind [which(abs(z1_ind) < pth1&abs(z2_ind) < pth1)];
      z2_new = z2_ind[which(abs(z1_ind) < pth1&abs(z2_ind) < pth1)];
      rhores = truncEstfun(a, b, z1_new, z2_new, 4000, 1000, 10)
      rhohat[i] = mean(rhores);
      p1 = length(z1_new);
      pvalue[i] = testR(rhohat[i], p1);
      Rhores[, i] = rhores;
      pres[i] = p1;
    }
  }
  # ---------------------------------------------------------
  return(list(rhohat = rhohat, pvalue = pvalue, pres = pres, Rhores = Rhores));
}