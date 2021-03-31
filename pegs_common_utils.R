### Common utilitarian constants and functions for analyzing the PEGS datasets ###
## Creator/Author: Farida Akhtari (farida.akhtari@nih.gov)


### Load required packages ###

# A package to load packages :)
if (!"librarian" %in% rownames(installed.packages())) {
  install.packages("librarian")
}
library(librarian)
librarian::shelf(data.table, naniar, tidyverse, wrapr)


### Common constants ###

# special codes
PEGS_SP_CODES <- qc(.M, .S, -444444, -555555, -666666, -777777, -888888, -999999)


### Common functions ###

# Convert columns in the PEGS dataframe to appropriate variable types as
# specified in the metadata file. Special codes are converted to NA.
# Arguments:-
# pegs_df       : The PEGS dataframe to convert
# pegs_df_meta  : The accompanying PEGS metadata dataframe
# Returns:-
# The converted PEGS dataframe
pegs_convert_type <- function(pegs_df, pegs_df_meta) {
  # All special codes need to be replaced with NA before type conversion
  pegs_df <- pegs_df %>% replace_with_na_all(condition = ~ .x %in% PEGS_SP_CODES)

  for (col in colnames(pegs_df)) {
    print(col)
    true_class <- pegs_df_meta %>%
      filter(long_variable_name == col) %>%
      pull(true_class)

    pegs_df[[col]] <- switch(true_class,
      "numeric" = as.numeric(pegs_df[[col]]),
      "binary" = as.factor(pegs_df[[col]]),
      "character" = as.character(pegs_df[[col]]),
      "factor" = as.factor(pegs_df[[col]]),
      "ordered factor" = as.factor(pegs_df[[col]]),
      pegs_df[[col]]
    )
  }
  return(pegs_df)
}

# For backward compatibility
epr_convert_type <- function(pegs_df, pegs_df_meta) pegs_convert_type(pegs_df, pegs_df_meta)


# Define the specified phenotype in the UK Biobank data. Identify cases
# and controls for the specified disease/phenotype based on survey data and
# disease-specific inclusion/exclusion criteria.
# Arguments:-
# ukb_df      : The UK Biobank data frame
# phenotype   : The phenotype string
# Returns:-
# A dataframe (f.eid = participant ID, Y = phenotype). All participants to be
# excluded have their phenotype marked as NA.
prepare_ukb_phenotype <- function(ukb_data, phenotype) {

  ## Extract disease columns
  # TODO: include cancer columns

  # Non-cancer illness code, self-reported
  ukb_illness <- ukb_data[, c(
    paste0("f.20002.0.", 0:33), paste0("f.20002.1.", 0:33),
    paste0("f.20002.2.", 0:33), paste0("f.20002.3.", 0:33)
  )]

  # ICD9 codes
  ukb_ICD9 <- ukb_data[, c(
    paste0("f.41203.0.", 0:27), paste0("f.41205.0.", 0:29),
    paste0("f.41271.0.", 0:46)
  )]

  # ICD10 codes
  ukb_ICD10 <- ukb_data[, c(
    paste0("f.40001.", 0:1, ".0"), paste0("f.40002.0.", 0:13),
    paste0("f.40002.1.", 0:13), paste0("f.41202.0.", 0:65),
    paste0("f.41204.0.", 0:183), paste0("f.41270.0.", 0:212)
  )]


  ## Identify cases and controls for specified phenotype:

  # Type 2 diabetes (T2D)
  if (phenotype == "T2D") {
    # Identify individuals with different subtypes:
    # (gets row numbers in the UK Biobank data for the matches)

    # Any diabetes
    ind_diabetes <- sort(unique(unlist(c(
      lapply(ukb_illness, function(x) which(x %in% 1220:1223)),
      lapply(ukb_ICD9, function(x) which(substr(x, 1, 3) == 250)),
      lapply(ukb_ICD10, function(x) which(substr(x, 1, 3) %in% paste0("E", 10:14)))
    ))))

    # Type 1 diabetes (T1D)
    ind_T1D <- sort(unique(unlist(c(
      lapply(ukb_illness, function(x) which(x == 1222)),
      lapply(ukb_ICD9, function(x) which(x %in% c(25001, 25011, 25021, 25091))),
      lapply(ukb_ICD10, function(x) which(substr(x, 1, 3) == "E10"))
    ))))

    # Type 2 diabetes
    ind_T2D <- sort(unique(unlist(c(
      lapply(ukb_illness, function(x) which(x == 1223)),
      lapply(ukb_ICD9, function(x) which(x %in% c(25000, 25010, 25020, 25090))),
      lapply(ukb_ICD10, function(x) which(substr(x, 1, 3) == "E11"))
    ))))

    print(paste("#ppts with any type of diabetes =", length(ind_diabetes)))
    print(paste("#ppts with T1D =", length(ind_T1D)))
    print(paste("#ppts with T2D =", length(ind_T2D)))

    ## Mark cases and controls

    # Initialize phenotype vector
    Y <- rep(0, nrow(ukb_data))
    # Exclude individuals with any diabetes subtype (from controls)
    Y[ind_diabetes] <- NA
    # Mark individuals with T2D as cases
    Y[ind_T2D] <- 1
    # Exclude individuals with T1D (from cases and controls)
    Y[ind_T1D] <- NA

    cat("Phenotype = ", phenotype, "\n")
    print(table(Y, exclude = NULL))

    # Create cases & controls dataframe
    ukb_pheno_prepd <- data.frame(f.eid = ukb_data$f.eid, Y = Y)

    return(ukb_pheno_prepd)
  }
  else {
    stop("Unrecognized phenotype string : 'phenotype'")
  }
}


# Define the phenotype in the PEGS/EPR data. Identify cases and controls for
# the specified disease/phenotype based on survey data and disease-specific
# inclusion/exclusion criteria. Missing vs skipped responses for child
# questions are correctly accounted for by examining the response for the
# parent question in the exclusion criteria.
# Arguments:-
# pegs_data   : The converted PEGS data frame (as done by pegs_convert_type())
# phenotype   : The phenotype string
# Returns:-
# A dataframe (epr_number = participant ID, Y = phenotype). All participants to
# be excluded have their phenotype marked as NA.
prepare_pegs_phenotype <- function(pegs_data, phenotype) {
  pheno_data <- NULL

  # Type 2 diabetes (T2D)
  if (phenotype == "T2D") {
    # Exclude participants with likely different disease etiology, i.e. those
    # with age at onset < 20 (these are likely Type 1 diabetes or MODY ppts) or
    # those with gestational diabetes only
    pegs_data <- pegs_data %>%
      mutate(onset_age_under_20 = (if_else(he_c022e_diabetes_age_CHILDQ < 20, 1, 0)))

    diabetes_exclusions <- pegs_data %>%
      filter(he_c022_diabetes_PARQ == 1) %>%
      filter(onset_age_under_20 == 1
      | is.na(onset_age_under_20)) %>%
      pull(epr_number)
    diabetes_exclusions <- sort(unique(c(
      diabetes_exclusions,
      (pegs_data %>%
        filter(he_c022_diabetes_PARQ == 1
        & .data[["_he_gender_"]] %in% c(1, NA)) %>%
        filter(he_c022a_diabetes_preg_CHILDQ == 1) %>%
        pull(epr_number))
    )))

    pegs_data <- pegs_data %>%
      mutate(Y = replace(
        he_c022_diabetes_PARQ,
        epr_number %in% diabetes_exclusions,
        NA
      ))

    # Put all variables that were used for phenotype preparation into a
    # dataframe
    pheno_data <- pegs_data %>%
      select(all_of(c(
        "epr_number", "he_c022_diabetes_PARQ", "onset_age_under_20",
        "_he_gender_", "he_c022a_diabetes_preg_CHILDQ", "Y"
      )))
  } else {
    stop(paste("Unrecognized phenotype string :"), phenotype)
  }

  cat(
    phenotype,
    "phenotype reclassification (Y) based on inclusion/exclusion criteria:", "\n"
  )
  pheno_data %>%
    select(-epr_number) %>%
    group_by_all() %>%
    count() %>%
    print(n = Inf)

  # Create cases(=1) & controls(=0) dataframe
  pegs_pheno_prepd <- pheno_data %>% select(epr_number, Y)
  cat("Phenotype = ", phenotype, "\n")
  print(table(pegs_pheno_prepd$Y, exclude = NULL))

  return(pegs_pheno_prepd)
}
