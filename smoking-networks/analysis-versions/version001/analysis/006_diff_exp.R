options(width = 200, max.print = 200)

pacman::p_load(tidyverse, magrittr, zeallot)

rm(list = ls()); gc()

"results/006" %>% {if(!dir.exists(.)) dir.create(.)}

#\\\
#\\\
# Load Data
#\\\
#\\\

# load adjacency matrices
c(adj_form, adj_curr) %<-% lapply(c("form", "curr"), function(nm){
    m = read.csv(sprintf('results/005/%s-adj.csv',nm))
    rnames = m[,1]
    m = m[,-1]
    m = data.matrix(m)
    dimnames(m) = list(rnames, rnames)
    return(m)
})


# metadata
rowData = read.csv('processed-data/001/rowData.csv' )
rowData[,1] = NULL



# \\\
# \\\
# Differential Expression with fisher test
# \\\
# \\\

diff_exp = function(signif_metab_ids,metadata){
    vars = c("sub_pathway", 'super_pathway')
    out = list()
    for(var in vars){
        signif_pathways <- metadata[metadata$metab_id %in% signif_metab_ids , var]
        exp <- as.data.frame(table(metadata[,var]))
        obs <- as.data.frame(table(Var1=signif_pathways))
        obs <- merge(exp, obs, by = "Var1", all.x = TRUE)
        names(obs) <- c(var,  "exp_freq" , "obs_freq")
        obs[is.na(obs$obs_freq),"obs_freq"] <- 0
        obs$fisher.p <- NA
        obs[[var]] = as.character(obs[[var]])
        for(sub in obs[,var]){
            observed = obs$obs_freq[obs[,var] == sub]
            observed_tot = sum(obs$obs_freq) 
            expected = obs$exp_freq[obs[,var] == sub]
            expected_tot = sum(obs$exp_freq)

            # E.G.
            ##             In Pathway    Not in Pathway
            ## Significant        215                95
            ## Non-significant    148               300
            tab <- matrix(ncol = 2, nrow = 2,  c(
                                    observed,                            (observed_tot - observed),
                            expected - observed, (expected_tot - observed_tot) - (expected - observed)))
            obs[obs[,var] == sub,"fisher.p"] <- fisher.test(tab, simulate.p.value = FALSE, alternative = "greater")$p.value
        }
        obs <- obs[order(obs$fisher.p),]
        obs$fdr.p <- p.adjust(obs$fisher.p, method = "fdr")
        obs <- obs[order(obs$fdr.p),]
        enrich <- obs[obs$fdr.p <.05,]
        out[[var]] <- obs
    }
    return(out)
}

signif_metab_ids = unique(rownames(adj_curr), rownames(adj_form))
metadata = rowData
var = "super_pathway"
sub = "Lipid"

res_form = diff_exp(rownames(adj_form), rowData)
res_curr = diff_exp(rownames(adj_curr), rowData)
res_intersect = diff_exp(unique(rownames(adj_curr), rownames(adj_form)), rowData)

"results/006/diff-exp" %>% {if(!dir.exists(.)) dir.create(.)}


write.csv( res_form$super_pathway , "results/006/diff-exp/former-superpathway.csv")
write.csv( res_form$sub_pathway , "results/006/diff-exp/former-subpathway.csv")

write.csv(res_curr$super_pathway, "results/006/diff-exp/current-superpathway.csv")
write.csv(res_curr$sub_pathway, "results/006/diff-exp/current-subpathway.csv")

write.csv(res_intersect$super_pathway, "results/006/diff-exp/intersection-superpathway.csv")
write.csv(res_intersect$sub_pathway, "results/006/diff-exp/intersection-subpathway.csv")

