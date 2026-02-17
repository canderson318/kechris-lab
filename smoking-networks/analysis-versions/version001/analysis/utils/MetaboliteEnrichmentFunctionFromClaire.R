## Function created by Claire Guo based on code 
## in file "sub_pathway_enrichment.R" from Suneeta Godbole
## metadata is the metabolite metadata, must have columns named "SUB_PATHWAY"
## and "SUPER_PATHWAY"

metab_pathway_fisher <- function(metadata,
                                 signif_metabids){
  
  out <- list()
  pathcols <- c("SUB_PATHWAY","SUPER_PATHWAY")
  for(p in 1:length(pathcols)){
    signif_pathways <- metadata[metadata$metabid %in% signif_metabids , pathcols[p]]
    
    exp <- as.data.frame(table(metadata[,pathcols[p]]))
    obs <- as.data.frame(table(Var1=signif_pathways))
    obs <- merge(exp, obs, by = "Var1", all.x = TRUE)
    names(obs) <- c(pathcols[p],  "exp.freq" , "obs.freq")
    obs[is.na(obs$obs.freq),"obs.freq"] <- 0
    obs$fisher.p <- NA
    
    
    for(sub in obs[,pathcols[p]]){
      cntg_tab <- matrix(ncol = 2, nrow = 2, 
                         c(obs$obs.freq[obs[,pathcols[p]] == sub],
                           (sum(obs$obs.freq) - obs$obs.freq[obs[,pathcols[p]] == sub]),
                           obs$exp.freq[obs[,pathcols[p]] == sub] - obs$obs.freq[obs[,pathcols[p]] == sub],      
                           
                           (sum(obs$exp.freq) - sum(obs$obs.freq) - (obs$exp.freq[obs[,pathcols[p]] == sub] - obs$obs.freq[obs[,pathcols[p]] == sub]))))
      obs[obs[,pathcols[p]] == sub,"fisher.p"] <- fisher.test(cntg_tab, simulate.p.value = FALSE, alternative = "greater")$p.value
      
    }
    
    obs <- obs[order(obs$fisher.p),]
    obs$fdr.p <- p.adjust(obs$fisher.p, method = "fdr")
    obs <- obs[order(obs$fdr.p),]
    enrich <- obs[obs$fdr.p <.05,]
    
    out[[pathcols[p]]] <- obs
  }
  
  out
}