ReadSummaryStat_thisstudy <- function(fileexp, fileout, filepan, snpinfo, pva_cutoff, lambad) {
  daexp = read.table(fileexp, header = T)
  daout = read.table(fileout, header = T)
  
  daexp = daexp[daexp$pvalue < pva_cutoff, ]
  daexp$a1 = 0
  daexp$a2 = 0
  daexp$a1[daexp$A1=="C"] = 1
  daexp$a1[daexp$A1=="G"] = 1
  daexp$a2[daexp$A2=="C"] = 1
  daexp$a2[daexp$A2=="G"] = 1
  daexp = daexp[daexp$a1+daexp$a2==1,]
  daexp[,c("A1", "A2", "a2", "pvalue", "BP")] = NULL
  daexp$chr = as.numeric(daexp$chr)
  daexp = unique(daexp)
  rownames(daexp) = daexp$SNP
  
  daout = daout[daout$SNP %in% daexp$SNP, ]
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
  daout = unique(daout)
  
  # daexp = daexp[daout$SNP,]
  # daout = daout[daout$a1==daexp$a1,]
  # daexp = daexp[daout$SNP,]
  
  daexp$SNP = as.character(daexp$SNP)
  daout$SNP = as.character(daout$SNP)
  # ===================================================== #
  # Remove MHC region.
  QCres = RemoveMHCCfun(mhcstart, mhcend, daexp$beta, daout$beta, daout$beta, daout$se, bp, daout$chr,
                        daout$SNP,  Inf, Inf)
  rsnameRMHC = QCres$rsnamenew
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
  
  ResF4gamma   = vector("list", 0)
  ResF4gammase = vector("list", 0)
  ResF4Gamma   = vector("list", 0)
  ResF4Gammase = vector("list", 0)
  ResF4SNPchr  = vector("list", 0)
  ResF4R       = vector("list", 0)
  Res_id  = c()
  
  for (chrk in unique(daexp$chr)[order(unique(daexp$chr))]) {
    if (is.null(filepan[[chrk]])) next
    if (!file.exists(filepan[[chrk]])) next
    dapan = readRDS(filepan[[chrk]])
    dapan[] <- lapply(dapan, function(x) if (is.factor(x)) as.character(x) else {x})
    
    daexp_k = daexp[daexp$chr==chrk,]
    if(nrow(daexp_k)==0) next
    daout_k = daout[daexp_k$SNP,]
    # ===================================================== #
    RefSNP_k = daRefsnpinfosame[daRefsnpinfosame$CHR==chrk, ]
    snpint = as.character(intersect(union(dapan$SNP1, dapan$SNP2), daout_k$SNP));
    stopifnot(sum(RefSNP_k$SNP%in%snpint)==length(snpint))
    # ===================================================== #
    
    dapan = dapan[dapan$SNP1 %in% daout_k$SNP,]
    dapan = dapan[dapan$SNP2 %in% daout_k$SNP,]
    daexp_k = daexp_k[daexp_k$SNP %in% snpint, ]
    daout_k = daout_k[daexp_k$SNP,]
    snppan = union(dapan$SNP1, dapan$SNP2)
    
    
    if(length(snppan)!=length(snpint)){
      cat("this is the case", chrk, "\n")
      sigleSNP = snpint[!snpint%in%snppan];
      Sigpan = data.frame(rep(chrk, length(sigleSNP)),paste0(chrk*10, rep(max(as.numeric(dapan$BlockID)),length(sigleSNP))+seq(1:length(sigleSNP))), 
                          sigleSNP, sigleSNP, rep(1, length(sigleSNP)))
      colnames(Sigpan) = colnames(dapan)
      dapan = rbind(dapan, Sigpan)
    }
    
    
    arma_blc_id = unique(dapan$BlockID)

    nid = length(arma_blc_id)
    
    if(nid==0) next
    
    F4gamma   = vector("list", nid)
    F4gammase = vector("list", nid)
    F4Gamma   = vector("list", nid)
    F4Gammase = vector("list", nid)
    F4SNPchr  = vector("list", nid)
    F4R       = vector("list", nid)
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
      
      F4gamma[[k]] = da0$beta
      F4gammase[[k]] = da0$se
      F4Gamma[[k]] = da1$beta
      F4Gammase[[k]] = da1$se
      F4SNPchr[[k]] = comsnp
      
      M = diag(rep(1,length(comsnp)))
      rownames(M) = comsnp
      colnames(M) = comsnp
      if(length(comsnp)!=1){
        for (j in 1:nrow(da2)) {
          M[da2$SNP1[j], da2$SNP2[j]] = da2$r[j]*lambad
          M[da2$SNP2[j], da2$SNP1[j]] = da2$r[j]*lambad
        }
      }
      rownames(M) = NULL
      colnames(M) = NULL
      F4R[[k]] = M
    }
    
    ResF4gamma    = append(ResF4gamma, F4gamma)
    ResF4gammase  = append(ResF4gammase, F4gammase)
    ResF4Gamma    = append(ResF4Gamma, F4Gamma)
    ResF4Gammase = append(ResF4Gammase, F4Gammase)
    ResF4SNPchr   = append(ResF4SNPchr, F4SNPchr)
    ResF4R        = append(ResF4R, F4R)
    Res_id        = c(Res_id, arma_blc_id)
  }
  
  return(
    list(
      ResF4gammah = ResF4gamma,
      ResF4se1 = ResF4gammase,
      ResF4Gammah = ResF4Gamma,
      ResF4se2 = ResF4Gammase,
      ResF4SNPchr = ResF4SNPchr,
      ResF4Rblock = ResF4R,
      arma_blc_id = Res_id
    )
  )
}
