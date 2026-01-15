## a function that accepts the count dataframe (e.g. metabolic measurements) 
## with subject IDs as a column and a clinical dataframe. 
## The function log-transforms the count data, adjusts the count data
## using the names of the covariates passed, and then standardizes the residuals

preprocess_and_adjust <- function(count_df,
                                  clinical_df,
                                  subject_col = "sid",
                                  covariates,
                                  log_transform = TRUE,
                                  log_offset = 0,
                                  scale_features = TRUE) {
  
  ## -----------------------------
  ## 0. Validate inputs
  ## -----------------------------
  if (!subject_col %in% colnames(count_df))
    stop("subject_col '", subject_col, "' not found in count_df")
  if (!subject_col %in% colnames(clinical_df))
    stop("subject_col '", subject_col, "' not found in clinical_df")
  
  missing_covars <- setdiff(covariates, colnames(clinical_df))
  if (length(missing_covars) > 0)
    stop("Covariates not found in clinical_df: ", 
         paste(missing_covars, collapse = ", "))
  
  ## -----------------------------
  ## 1. Match subjects
  ## -----------------------------
  common_ids <- intersect(count_df[[subject_col]],
                          clinical_df[[subject_col]])
  
  if (length(common_ids) == 0)
    stop("No overlapping subject IDs found.")
  
  # Subset to common subjects and match order
  count_df <- count_df[count_df[[subject_col]] %in% common_ids, ]
  clinical_df <- clinical_df[clinical_df[[subject_col]] %in% common_ids, ]
  
  count_df <- count_df[match(common_ids, count_df[[subject_col]]), ]
  clinical_df <- clinical_df[match(common_ids, clinical_df[[subject_col]]), ]
  
  ## -----------------------------
  ## 2. Remove subjects with missing covariates
  ## -----------------------------
  # Now that order matches, we can use complete.cases and drop subjects with missing covariates
  covar_df <- clinical_df[, covariates, drop = FALSE]
  cc <- complete.cases(covar_df)
  
  if (sum(cc) == 0)
    stop("No subjects with complete covariate data.")
  
  if (sum(!cc) > 0) {
    message("Removing ", sum(!cc), " subjects with missing covariate data")
  }
  
  count_df <- count_df[cc, ]
  clinical_df <- clinical_df[cc, ]
  covar_df <- covar_df[cc, , drop = FALSE]
  
  ## -----------------------------
  ## 3. Extract feature matrix
  ## -----------------------------
  feature_mat <- as.matrix(
    count_df[, setdiff(colnames(count_df), subject_col)]
  )
  rownames(feature_mat) <- count_df[[subject_col]]
  
  ## -----------------------------
  ## 4. Log transform
  ## -----------------------------
  if (log_transform) {
    feature_mat <- log(feature_mat + log_offset)
  }
  
  ## -----------------------------
  ## 5. Regress out covariates
  ## -----------------------------
  adjusted_mat <- matrix(NA,
                         nrow = nrow(feature_mat),
                         ncol = ncol(feature_mat))
  colnames(adjusted_mat) <- colnames(feature_mat)
  rownames(adjusted_mat) <- rownames(feature_mat)
  
  # Create design matrix (automatically handles factors)
  design <- model.matrix(~ ., data = covar_df)
  
  for (j in seq_len(ncol(feature_mat))) {
    fit <- lm.fit(design, feature_mat[, j])
    adjusted_mat[, j] <- fit$residuals
  }
  
  ## -----------------------------
  ## 6. Standardize features
  ## -----------------------------
  if (scale_features) {
    adjusted_mat <- scale(adjusted_mat)
  }
  
  ## -----------------------------
  ## 7. Return results
  ## -----------------------------
  return(list(
    adjusted_matrix = adjusted_mat,
    subject_ids = rownames(adjusted_mat),
    covariates_used = covariates,
    n_subjects = nrow(adjusted_mat),
    n_features = ncol(adjusted_mat),
    n_removed_missing = sum(!cc)
  ))
}

# ## using the function for pre-processing the metabolic measurements
# out <- preprocess_and_adjust(
#   count_df    = countsFiltered,
#   clinical_df = clinicalFiltered,
#   subject_col = "sid",
#   covariates  = c("age_visit", "gender", "ccenter", "BMI")
# )
# 
# ## extracting the pre-processed data
# adjusted_features <- out$adjusted_matrix
