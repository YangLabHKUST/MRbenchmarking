setwd("MR-benchmarking-data/dataset1")

library(readr)
library(MRAPSS)
ts1= load("outcome_set.RData")
ts2 = load("outcome_set.RData")

for( exposure in ts1){
  for( outcome in ts2){
    # Start the clock!
    start = proc.time()
    
    # read in formatted GWAS data
    trait1.dir = paste0("./Formatted/exposure.set/", exposure)
    trait2.dir = paste0("./Formatted/outcome.set/", outcome)
    
    cat(exposure,"~",outcome,"\n")
    
    trait1 = readr::read_delim(trait1.dir, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
    
    trait2 = readr::read_delim(trait2.dir, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = F)
    
    # estimate parameters for the background model
    paras = est_paras(dat1 = trait1,
                      dat2 = trait2,
                      trait1.name = exposure,
                      trait2.name = outcome,
                      h2.fix.intercept = F,
                      LDSC = T,
                      ldscore.dir = "MR-benchmarking-data/eur_w_ld_chr")
    
    if(inherits(paras, 'try-error')) next
    
    write.table(matrix(as.vector(paras$Omega), nrow=1), paste0("./bg_paras/1kgRef_", exposure, "~", outcome, "_Omega"),
                row.names = F, col.names = F, quote = F)
    
    write.table(matrix(as.vector(paras$C), nrow=1), paste0("./bg_paras/1kgRef_", exposure, "~", outcome, "_C"),
                row.names = F, col.names = F, quote = F)
    
    # clumping      
    cat("Begin clumping ...\n ")
    
    
    # For convenience, users can perform LD clumping through API and do not need to provide a reference panel, i.e.,
    clumped =  clump(paras$dat,
                     IV.Threshold = 5e-05,
                     SNP_col = "SNP",
                     pval_col = "pval.exp",
                     clump_kb = 1000)
    
    
    if(inherits(clumped , 'try-error')) next
    
    
    write.table(clumped, file = paste0("./MRdat/1kgRef_", exposure,"~",outcome), col.names=T, row.names=F, quote=F)
    
    # Stop the clock
    proc.time() - start
    
  }
  
}

