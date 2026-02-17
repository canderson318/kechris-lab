options(width = 200, max.print = 200)

pacman::p_load(tidyverse, magrittr, zeallot)

rm(list = ls()); gc()

"results/007" %>% 
    {if(!dir.exists(.)) dir.create(.)}
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

# load hubs

extract_metabs  <-  function(list){
    lapply(list, function(S){
        M = unlist(str_split(S, "\\|\\|"))
        m = lapply(M, function(s){
            str_split(s, '•') %>% 
                lapply(., `[`, 3)  %>% 
                lapply(., function(x){
                    str_trim(x)
                })
        }) %>% 
        unlist()
    })
}

hub_metabs_curr = extract_metabs(readLines('results/005/curr-hub-metabs/metabs.txt') )
hub_metabs_form = extract_metabs(readLines('results/005/form-hub-metabs/metabs.txt') )

# load counts 
counts_curr = data.matrix(read.csv('processed-data/002/separate/current.csv', header = FALSE))
counts_form = data.matrix(read.csv('processed-data/002/separate/former.csv', header = FALSE))
dim(counts_form)
dim(counts_curr)
colnames(counts_curr) =  rowData$metab_id
colnames(counts_form) =  rowData$metab_id

#\\\
#\\\
# look at correlations
#\\\
#\\\

dir = "results/007/hub-metab-cor-hms/"

dir %>% 
    {if(!dir.exists(.)) dir.create(.)}

# install.packages("raster")
# install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
for(j in 1:length(mats)){
    for(i in 1:length(hub_metabs)){
        mats = list("Current" = counts_curr, "Former" = counts_form)
        metabs = hub_metabs[[i]]
        if(is.na(metabs)[1] | (length(metabs) == 0)){
            next 
        }
        counts = mats[[j]]
        hub_metab_ids = rowData$metab_id[rowData$chemical_name %in% metabs]
        sub = counts[, which(colnames(counts) %in% hub_metab_ids)]
        colnames(sub) <-  metabs
        corr = cor(sub)
        pdf(sprintf("%s%s-%s.pdf", dir, names(mats)[j], i), height = 10, width = 10)
        hm  = ComplexHeatmap::Heatmap(corr, name = "rho", column_title = sprintf("%s smokers, hub %s",names(mats)[j], i) )
        ComplexHeatmap::draw(hm)
        dev.off()
    }
}
