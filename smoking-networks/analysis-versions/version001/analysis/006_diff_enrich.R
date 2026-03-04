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

# load components
load_metabs  <-  function(list){
    lapply(list, function(S){
        M = unlist(str_split(S, "\\|\\|"))
        m = lapply(M, function(s){
            str_split(s, '•') %>% 
                lapply(., `[`, 3)  %>% # grab third element, chemical name
                lapply(., function(x){
                    str_trim(x)
                })
        }) %>% 
        unlist()
    })
}

get_metab_id <- function(comp_list){
    lapply(comp_list, function(comp){
        sapply(comp, function(nm){
            rowData$metab_id[rowData$chemical_name == nm]
        }, USE.NAMES = FALSE)
    })
}

comp_curr = readLines("results/005/curr-metabs/comp-metabs.txt") %>% 
    load_metabs() %>% 
    get_metab_id()
comp_form = readLines("results/005/form-metabs/comp-metabs.txt") %>% 
    load_metabs() %>% 
    get_metab_id()

node_unique_curr = readLines("results/005/curr-metabs/unique-nodes.txt") %>% 
    load_metabs() %>% 
    get_metab_id()
node_unique_form = readLines("results/005/form-metabs/unique-nodes.txt") %>% 
    load_metabs() %>% 
    get_metab_id()


# \\\
# \\\
# Differential Enrichment with fisher test
# \\\
# \\\

enrichment = function(signif_metab_ids,metadata){
    vars = c("sub_pathway", 'super_pathway')
    out = list()
    for(var in vars){
        signif_pathways <- metadata[metadata$metab_id %in% signif_metab_ids , var]
        exp <- as.data.frame(table(metadata[,var]))
        res <- as.data.frame(table(Var1=signif_pathways))
        res <- merge(exp, res, by = "Var1", all.x = TRUE)
        names(res) <- c(var,  "exp" , "obs")
        res[is.na(res$obs),"obs"] <- 0
        
        res$OR <- NA
        res$fisher.p <- NA
        res$obs_in <- NA
        res$obs_out <- NA
        res$unobs_in <- NA
        res$unobs_out <- NA
        res[[var]] = as.character(res[[var]])

        for(sub in res[,var]){
            observed = res$obs[res[,var] == sub] # count metabs observed in pathway
            observed_tot = sum(res$obs) # Total observed metabolites
            observable = res$exp[res[,var] == sub] # Count metabolites observable in pathay
            observable_tot = sum(res$exp) # total observable metabolites

            # E.G.
            ##             In Pathway    Not in Pathway
            ## observed           215                95
            ## unobserved         148               300
            tab <- matrix(ncol = 2, nrow = 2,  byrow = TRUE,
                c(observed,            (observed_tot - observed),
                  observable - observed, (observable_tot - observed_tot) - (observable - observed)
                 ))

            res[res[,var] == sub, c("obs_in", "obs_out", "unobs_in", "unobs_out") ] = tab[matrix(c(1,1,1,2,2,1,2,2),ncol = 2, byrow = TRUE)]
            
            # fill dataframe with p values
            fish_test = fisher.test(tab, simulate.p.value = FALSE, alternative = "greater")
            res[res[,var] == sub,"fisher.p"] <- signif(fish_test$p.value,3)
            res[res[,var] == sub,"OR"] <- as.numeric(fish_test$estimate) %>% round(3)
        }
        res <- res[order(res$fisher.p),]
        res$fdr.p <- signif(p.adjust(res$fisher.p, method = "fdr"),3)
        res$logOR <- log(res$OR) %>% round(3)
        res <- res[order(res$fdr.p),]
        res <- res[, c(var, "exp", "obs", "obs_in", "obs_out", "unobs_in", "unobs_out", "OR", "logOR", "fisher.p", "fdr.p")]
        out[[var]] <- res
    }
    return(out)
}
# signif_metab_ids =node_unique_form[[1]]
# metadata = rowData
# var = "sub_pathway"
# sub = "Corticosteroids"
# a <- tab[1,1] # "Cases WITH the thing I'm testing for"
# b <- tab[1,2] # "Cases WITHOUT the thing"
# c <- tab[2,1] # "Controls WITH the thing"
# d <- tab[2,2] # "Controls WITHOUT the thing"
# OR <- (a * d) / (b * c)


# overall comparisons between each groups metabs to background
res_node_form = enrichment(rownames(adj_form), rowData)
res_node_curr = enrichment(rownames(adj_curr), rowData)
res_node_intersect = enrichment(unique(rownames(adj_curr), rownames(adj_form)), rowData)

# Metabs unique to each group
res_node_unique_form =enrichment(node_unique_form[[1]], rowData)
res_node_unique_curr =enrichment(node_unique_curr[[1]], rowData)

# Components
res_node_comp_form  <-  list()
for(i in seq(comp_form)){
    component = comp_form[[i]]
    res_node_comp_form[[i]] =enrichment(component, rowData)
}
res_node_comp_curr  <-  list()
for(i in seq(comp_curr)){
    component = comp_curr[[i]]
    res_node_comp_curr[[i]] =enrichment(component, rowData)
}

res_node_comp_curr[[3]]$sub_pathway[1,]
res_node_comp_form[[3]]$sub_pathway[1,]

# \\\\
# \\\\
# Save results
# \\\\
# \\\\

"results/006/enrichment" %>% {if(!dir.exists(.)) dir.create(.)}

save_res <- function(res, dir,nm){
    super = res$super_pathway %>% 
        filter(fisher.p <= .05)  
    write.csv(super , sprintf("%s%s-superpathway.csv", dir, nm), row.names = FALSE)

    sub = res$sub_pathway %>% 
        filter(fisher.p <= .05) 
    write.csv( sub , sprintf("%s%s-subpathway.csv",dir, nm), row.names = FALSE)
}

save_res(res_node_curr, dir="results/006/enrichment/", nm = 'current')
save_res(res_node_form, dir="results/006/enrichment/", nm = 'former')
save_res(res_node_intersect, dir="results/006/enrichment/", nm = 'intersection')
save_res(res_node_unique_curr, dir="results/006/enrichment/", nm = 'curr-unique')
save_res(res_node_unique_form, dir="results/006/enrichment/", nm = 'form-unique')



"results/006/enrichment/curr-components" %>% {if(!dir.exists(.)) dir.create(.)}
for(i in seq(res_node_comp_curr)){
    save_res(res_node_comp_curr[[i]], dir = 'results/006/enrichment/curr-components/', nm = i)
}

"results/006/enrichment/form-components" %>% {if(!dir.exists(.)) dir.create(.)}
for(i in seq(res_node_comp_form)){
    save_res(res_node_comp_form[[i]], dir = 'results/006/enrichment/form-components/', nm = i)
}

