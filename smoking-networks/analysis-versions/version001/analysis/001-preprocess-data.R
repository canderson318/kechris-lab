
pacman::p_load(
  dplyr, readr, stringr, tidyr, ggplot2, magrittr, zeallot, tibble, janitor, openxlsx,
               install = F)

rm(list= ls()); gc()


# wd_path = '/Users/canderson/Documents/school/local-kechris-lab/rotation-project/analysis-versions/version001'
wd_path = '/projects/canderson2@xsede.org/kechris-lab/smoking-networks/analysis-versions/version001'

tryCatch(
  setwd(wd_path),
  error = function(e) {
    message("Failed to set working directory: ", e$message)
    quit(status = 1)
  }
)


#///
#///
# Load Data
#///
#///

# system("cat ../../info/info.txt")

# This file contains the metadata of each metabolite
suppressMessages(
  rowData <- read_csv('../../raw-data/COPDGene_P2_MetaboliteInformation_20211021.csv') %>% 
    clean_names()
)

# This file contains the relevant clinical data from Phase 2 for these subjects. 
suppressMessages(
  colData <- read_csv('../../raw-data/COPDGene_P1P2P3_SM_NS_Long_Oct23_ShareWithChristian.csv') %>% 
    clean_names()
)

# This file contains metabolite measurements from 995 metabolites for 1125 subjects from Phase 2 of COPDGene
suppressMessages(
  counts <- read_csv('../../raw-data/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021_ShareWithChristian.csv') %>% 
    clean_names()
)

# data dictionary
dat_dict <- read.xlsx("../../raw-data/DataDictionary.xlsx")



#///
#///
# Process data
#///
#///

#\\
# prefilter data
#\\

# make ids lowercase like counts names
rowData$metab_id <- tolower(rowData$metab_id)

# make sex field
colData$sexf <- colData$gender==2

# where smokers
subj_smoked = colData$sid[colData$smoking_status > 0] # 1060 subjects

# where gold 0
subj_gold0  = colData$sid[colData$finalgold_visit == 0] # 467 subjects

# subjects that meet both conditions
subj_keep = intersect(subj_gold0,subj_smoked) # 448 

colData = colData[colData$sid %in% subj_keep , ] 
nrow(colData) # 447

counts = counts[counts$sid %in% subj_keep , ]
nrow(counts) # 447


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


#///
#///
# Check dims match between matrices
#///
#///

# filter rowData for metabs in counts
rowData = rowData[rowData$metab_id %in% rownames(adjusted_matrix)  ,  ]

# filter colData for sids in counts
colData = colData[ colData$sid %in% colnames(adjusted_matrix) , ]

if(!(nrow(rowData) == nrow(adjusted_matrix)) & !(nrow(colData) == ncol(adjusted_matrix)) )   stop("Dims Don't Match!!!")
if(!all(identical(colData$sid, colnames(adjusted_matrix)), identical(rowData$metab_id, rownames(adjusted_matrix))) )  stop("Dimnames do not match row/col data!!!")

# include raw counts too
counts_filtered <- data.matrix(t(counts[, -1]))
colnames(counts_filtered) <- counts$sid
counts_filtered <- counts_filtered[rownames(counts_filtered)%in%rownames(adjusted_matrix), 
                                   colnames(counts_filtered)%in%colnames(adjusted_matrix)]


if(!all(identical(rownames(counts_filtered), rownames(adjusted_matrix)), identical(colnames(counts_filtered), colnames(adjusted_matrix)) ) ) stop("Raw counts and adjusted do not match!!!")

# ///
# ///
# Filter out unamed metabolites
# ///
# ///

# unamed have 'X-\d+' pattern
unnamed <- rowData$chemical_name %>% 
  grep("^X-\\d+", ., value = TRUE)

unnamed_metab_ids <- rowData$metab_id[rowData$chemical_name %in% unnamed]

adjusted_matrix <- adjusted_matrix[!rownames(adjusted_matrix) %in% unnamed_metab_ids, ]
rowData <- rowData[!rowData$metab_id %in% unnamed_metab_ids, ]

nrow(adjusted_matrix)- length(unnamed) # 758
dim(adjusted_matrix) # 758 1117


#///
#///
# Save Data
#///
#///

"processed-data/001/" %>% 
  {if(!dir.exists(.)) dir.create(.)}


# each table separately 
# write.csv(counts_filtered, 'processed-data/001/raw_counts.csv')
write.csv(adjusted_matrix, 'processed-data/001/adjusted_logcounts.csv')
write.csv(rowData,         'processed-data/001/rowData.csv')
write.csv(colData,         'processed-data/001/colData.csv')

