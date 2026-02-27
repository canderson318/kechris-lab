options(width = 200, max.print = 200)

pacman::p_load(tidyverse, magrittr, zeallot, ComplexHeatmap)

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

# load components
extract_metabs  <-  function(list){
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

all_comp_metabs = list("Current" = extract_metabs(readLines('results/005/curr-metabs/comp-metabs.txt') ), 
                    "Former" = extract_metabs(readLines('results/005/form-metabs/comp-metabs.txt') ))
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

dir = "results/007/comp-metab-cor-hms/"

dir %>% 
    {if(!dir.exists(.)) dir.create(.)}

# install.packages("raster")
# install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

mats = list("Current" = counts_curr, "Former" = counts_form)
hms = list()

# make  Heatmaps
for(j in 1:length(mats)){
    group = names(mats)[j]
    # grab groups components
    comp_metabs  = all_comp_metabs[[group]]
    
    # for each compoonent plot correlations
    for(i in 1:length(comp_metabs)){
        metabs = comp_metabs[[i]]
        
        if(is.na(metabs)[1] | (length(metabs) == 0)) next

        counts = mats[[j]]
        comp_metab_ids = rowData$metab_id[rowData$chemical_name %in% metabs]
        sub = counts[, which(colnames(counts) %in% comp_metab_ids)]
        colnames(sub) <-  metabs
        corr = cor(sub)
        hm  = Heatmap(corr, name = "rho", column_title = sprintf("%s smokers, comp %s",group, i), cluster_columns = FALSE,cluster_rows = FALSE, col = viridis::plasma(100))
        hms[[group]][[i]] <- hm
    }
}


# Draw Heatmaps
for(j in 1:length(mats)){
    group = names(mats)[j]
    # for each compoonent plot correlations
    for(i in 1:length(comp_metabs)){
        if(is.na(metabs)[1] | (length(metabs) == 0)) next
        pdf(sprintf("%s%s-%s.pdf", dir, group, i), height = 10, width = 10)
        ComplexHeatmap::draw(hms[[group]][[i]])
        dev.off()
    }
}

# #\\\
# #\\\
# # Plot Adjacency at component nodes
# #\\\
# #\\\

# get_metab_id <- function(comp_list){
#     lapply(comp_list, function(comp){
#         sapply(comp, function(nm){
#             rowData$metab_id[rowData$chemical_name == nm]
#         }, USE.NAMES = TRUE)
#     })
# }


# all_comp_metab_ids = sapply(all_comp_metabs, get_metab_id, simplify = FALSE)
# metabs = all_comp_metab_ids$Current[[1]]

# m = adj_curr[which(rownames(adj_curr) %in%metabs), which(colnames(adj_curr) %in%metabs) ]
# dimnames(m)  <- list(names(metabs)[match(metabs, colnames(m))], names(metabs)[match(metabs, colnames(m))])

# Heatmap(m, cluster_rows = FALSE, col = viridis::plasma(100))

# library(igraph)
# g = graph_from_adjacency_matrix(m, mode = 'undirected') 

# plot(g)
