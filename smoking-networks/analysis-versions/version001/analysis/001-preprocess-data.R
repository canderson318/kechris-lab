
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!'SummarizedExperiment' %in% installed.packages())
  BiocManager::install('SummarizedExperiment')

if(!'qs' %in% installed.packages())
  install.packages("qs")


pacman::p_load(
  dplyr, readr, stringr, tidyr, ggplot2, magrittr, zeallot, tibble, SummarizedExperiment, qs, janitor,
               install = F)

rm(list= ls()); gc()

setwd('/Users/canderson/Documents/school/local-kechris-lab/rotation-project/analysis-versions/version001')

#///
#///
# Load Data
#///
#///

system("cat ../../info/info.txt")

# This file contains the metadata of each metabolite
rowData <- read_csv('../../raw-data/COPDGene_P2_MetaboliteInformation_20211021.csv') %>% 
  clean_names()

# This file contains the relevant clinical data from Phase 2 for these subjects. 
colData <- read_csv('../../raw-data/COPDGene_P1P2P3_SM_NS_Long_Oct23_ShareWithChristian.csv') %>% 
  clean_names()

# This file contains metabolite measurements from 995 metabolites for 1125 subjects from Phase 2 of COPDGene
counts <- read_csv('../../raw-data/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021_ShareWithChristian.csv') %>% 
  clean_names()

# data dictionary
dat_dict <- openxlsx::read.xlsx("../../raw-data/DataDictionary.xlsx")

# make ids lowercase like counts names
rowData$metab_id <- tolower(rowData$metab_id)

# make sex field
colData$sexf <- colData$gender==2


#///
#///
# Process data
#///
#///
#### utils/DataPreProcessCodes.R that contains a function which can allow you to pre-process the data by log-transformation, 
#### adjusting for necessary covariates, and standardization of the features. 
#### For your analysis, you can adjust the metabolite measurements for 
#### age, sex, BMI, clinic center, total leukocyte counts (wbc), percentage of each leukocyte type, and hemoglobin. 
#### I have included these variables in the clinical data, apart from FEV1 and GOLD stage variables.

source('analysis/utils/DataPreProcessCodes.R')

covariates <- c( "ccenter", "age_visit", "sexf", "bmi", "wbc", "neutrophl_pct", "eosinphl_pct", "lymphcyt_pct", "monocyt_pct", "hemoglobin")

c(adjusted_matrix,
  subject_ids,
  covariates_used,
  n_subjects,
  n_features,
  n_removed_missing) %<-% preprocess_and_adjust(count_df = counts, 
                                                clinical_df = colData, 
                                                subject_col = "sid",
                                                covariates = covariates,
                                                log_transform = TRUE,
                                                log_offset = 0,
                                                scale_features = TRUE
                                                )
dim(adjusted_matrix) # subs x metab
length(subject_ids)

# make counts p x n
adjusted_matrix = t(adjusted_matrix)


# plot(density(adjusted_matrix[1,]),col = rgb(0,0,0,.1), ylim = c(0, 1), main = 'Transformed and Adjusted Metabolite Counts')
# 
# apply(adjusted_matrix[-1,], 1, function(x){
#   lines(density(x),col = rgb(0,0,0,.1))
# })



#///
#///
# Format summarized experiment object
# Format summarized experiment object
#///
#///

# filter rowData for metabs in counts
rowData = rowData[rowData$metab_id %in% rownames(adjusted_matrix)  ,  ]

# filter colData for sids in counts
colData = colData[ colData$sid %in% colnames(adjusted_matrix) , ]

if(!(nrow(rowData) == nrow(adjusted_matrix)) & !(nrow(colData) == ncol(adjusted_matrix)) ){
  stop("Dims Don't Match!!!")
}

if(!all(identical(colData$sid, colnames(adjusted_matrix)), identical(rowData$metab_id, rownames(adjusted_matrix))) ){
  stop("Dimnames do not match row/col data!!!")
}

# include raw counts too
counts_filtered <- data.matrix(t(counts[, -1]))
colnames(counts_filtered) <- counts$sid
counts_filtered <- counts_filtered[rownames(counts_filtered)%in%rownames(adjusted_matrix), 
                                   colnames(counts_filtered)%in%colnames(adjusted_matrix)]


if(!all(identical(rownames(counts_filtered), rownames(adjusted_matrix)), 
        identical(colnames(counts_filtered), colnames(adjusted_matrix)) ) ){
  stop("Raw counts and adjusted do not match!!!")
}


se <- SummarizedExperiment(assays = list(counts = counts_filtered, adjusted_logcounts = adjusted_matrix),
                           rowData = rowData, 
                           colData = colData)

# add datadict
metadata(se) <- list(data_dictionary = dat_dict)


# ///
# ///
# Filter out unamed metabolites
# ///
# ///

# unamed have 'X-\d+' pattern
unnamed <- rowData(se)$chemical_name %>% 
  grep("^X-\\d+", ., value = TRUE)

unnamed_metab_ids <- rowData(se)$metab_id[rowData(se)$chemical_name %in% unnamed]

se <- se[!rownames(se) %in% unnamed_metab_ids, ]

nrow(adjusted_matrix)- length(unnamed) # 758
dim(se) # 758 1117


#///
#///
# Save Data
#///
#///

"processed-data/001/" %>% 
  {if(!dir.exists(.)) dir.create(.)}


qsave(x = se, file = "processed-data/001/se001.qs", preset = 'high')

# each table separately 
write.csv(assays(se)$counts, 'processed-data/001/raw_counts.csv')
write.csv(assays(se)$adjusted_logcounts, 'processed-data/001/adjusted_logcounts.csv')
write.csv(rowData(se), 'processed-data/001/rowData.csv')
write.csv(colData(se), 'processed-data/001/colData.csv')

