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
PEGS_SP_CODES <- qc(
  .M, .S, .N,
  -444444, -555555, -666666, -777777, -888888, -999999
)


### Common functions ###

# Convert columns in the PEGS dataframe to appropriate variable types as
# specified in the metadata file. Special codes are converted to NA.
# Arguments:-
# pegs_df       : The PEGS dataframe to convert
# pegs_df_meta  : The accompanying PEGS metadata dataframe
# Returns:-
# The converted PEGS dataframe
pegs_convert_type_original <- function(pegs_df, pegs_df_meta) {
  # All special codes need to be replaced with NA before type conversion
  pegs_df <- pegs_df %>%
    replace_with_na_all(condition = ~ .x %in% PEGS_SP_CODES)

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
      "date" = as.Date(pegs_df[[col]]),
      pegs_df[[col]]
    )
  }
  return(pegs_df)
}


#' PEGS data converter
#'
#' Convert columns  in the  PEGS data  frame to  appropriate variable  types as
#' specified in the metadata file.  Special codes are converted to NA.
#'
#' @param pegs_df The PEGS dataframe to convert
#' @param pegs_df_meta The accompanying PEGS metadata dataframe
#' @param quiet do not report conversion progress (def=TRUE)
#' @return The converted PEGS dataframe
pegs_convert_type <- function(pegs_df, pegs_df_meta, quiet=TRUE)
{
    ## list for converters, add character="identity"
    CONVERTERS <- c(binary=as.factor, factor=as.factor, numeric=as.numeric,
                    `ordered factor`=as.factor, date=as.Date, character=identity)
    ## set NA
    pegs_df <- as.matrix(pegs_df)
    pegs_df[pegs_df %in% PEGS_SP_CODES] <- NA
    pegs_df <- as.data.frame(pegs_df)
    
    ## convert to R types
    ## look up the meta-data long variable names with table headers, gets true class
    true_class <- with(pegs_df_meta,
                       true_class[match(names(pegs_df), long_variable_name)])
    ## look up the list of converter with true classes, gets converters
    converters <- CONVERTERS[true_class]
    for(i in seq_along(pegs_df))
    {
        pegs_df[[i]] <- converters[[i]](pegs_df[[i]])
        if(!quiet)
            cat(sprintf("%4i %-40s %16s -> %s\n",
                        i, names(pegs_df)[i], true_class[i], class(pegs_df[[i]])))
    }
    return(pegs_df)
}

#' check the equivalence between coverters
#'
#' this works for PEGS Data Freeze V2.
#' All three test must return TRUE.
pegs_convert_type_test <- function()
{
    pgs_dir <- "/ddn/gs1/project/controlled/PEGS/Data_Freezes/freeze_v2"
    ## test health and exposure
    load(file.path(pgs_dir, "Surveys/Health_and_Exposure/healthexposure_26aug21_v2.RData"))
    he1 <- pegs_convert_type_original(epr.he, epr.he.meta)
    he2 <- pegs_convert_type(epr.he, epr.he.meta, quiet=FALSE)
    cat("he1 == he2: ", all.equal(as.data.frame(he1), he2), "\n", sep="")

    ## test exposome A
    load(file.path(pgs_dir, "Surveys/Exposome/exposomea_02jun21_v2.RData"))
    ea1 <- pegs_convert_type_original(epr.ea, epr.ea.meta)
    ea2 <- pegs_convert_type(epr.ea, epr.ea.meta, quiet=FALSE)
    cat("ea1 == ea2: ", all.equal(as.data.frame(ea1), ea2), "\n", sep="")
    
    ## test exposome B
    load(file.path(pgs_dir, "Surveys/Exposome/exposomeb_02jun21_v2.RData"))
    eb1 <- pegs_convert_type_original(epr.eb, epr.eb.meta)
    eb2 <- pegs_convert_type(epr.eb, epr.eb.meta, quiet=FALSE)
    cat("eb1 == eb2: ", all.equal(as.data.frame(eb1), eb2), "\n", sep="")

    return(invisible(NULL))
}

# For backward compatibility
epr_convert_type <- function(pegs_df, pegs_df_meta) pegs_convert_type(pegs_df, pegs_df_meta)


# Convert the specified label string from any PEGS metadata file to a named
# character vector containing key-value pairs.
# Arguments:-
# labels  : character, the label string from a PEGS metadata file
# Returns:-
# A named character vector containing key-value pairs
# Example usage to recode gender values:
# gender_labels <- create_label_vector(
#   epr.bcbb.map.meta %>% filter(variable_name == "gender") %>% pull(label)
# )
# epr.bcbb.map.labeled <- epr.bcbb.map.conv %>%
#   mutate(gender = recode_factor(gender, !!!gender_labels))
create_label_vector <- function(labels) {

  # Separate the multiple labels in the label string (at ';') and then separate
  # each label into key and value pairs (at '=').
  kv_pairs <- strsplit(unlist(str_split(labels, ";")), "=")

  # Remove unwanted characters from the keys and values
  kv_pairs <- lapply(kv_pairs, function(x) {
    x <- trimws(x)
    x <- gsub("(^')|('$)", "", x)
    return(x)
  })

  # Convert to a two-column data frame with the keys in one column and the
  # values in the other column.
  kv.df <- do.call(rbind.data.frame, kv_pairs)
  names(kv.df) <- c("key", "value")

  # Convert to a named character vector containing key-value pairs
  label_vector <- kv.df %>%
    pmap(~ set_names(..2, ..1)) %>%
    unlist()

  return(label_vector)
}


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
  } else {
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
prepare_pegs_phenotype <- function(pegs_data, 
                                   phenotype = c("lower_gi_polyps",'fibroids',
                                                 'boneloss','migraines','IDA',
                                                 'ovariancysts','asthma',"T2D","Allergic_Rhinitis",
                                                 'Stroke','Heart_Attack','Arrhythmia','CAD','CHF',
                                                 'Cholesterol','Hypertension','AVSD','Angina','Angioplasty',
                                                 'PBF','Bclots','Atherogenic')) {
  pheno_data <- NULL

  .lower_gi_phenotype <- function(pegs_data){
    
    cancer_cols <- names(pegs_data)[grep('cancer',names(pegs_data),ignore.case = T)]
    cancer_ex <- c('age|sis|bro|dad|mom')
    cancer_cols <- cancer_cols[-grep(cancer_ex,cancer_cols,ignore.case = T)]
    gi_cols <- c('he_f038_lactose_intolerance','he_f039_crohns','he_f040_ulcerative_colitis')
    exclusion_cols <- c(cancer_cols,gi_cols)
    ex_df <- pegs_data[c(exclusion_cols)] 
    ex_df[] <- lapply(ex_df, function(x) as.numeric(as.character(x))) 
    d <- ex_df %>% 
      as.data.frame(.) %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- pegs_data$epr_number
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- pegs_data %>% dplyr::select(epr_number,he_f041_polyps,he_age_derived)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Lower_GI_Polyps = case_when(he_f041_polyps == 1 ~ 1,
                                         he_f041_polyps == 0 & Flag == 0 ~ 0,
                                         TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Lower_GI_Polyps)
    
    return(pheno_data)
    
  }
  
  .fibroids_phenotype <- function(pegs_data){
    f_he <- pegs_data %>% dplyr::filter(`_he_gender_`==1)
    f_he[] <- lapply(f_he, function(x) as.numeric(as.character(x))) 
    cancer_cols <- names(pegs_data)[grep('cancer',names(pegs_data),ignore.case = T)]
    cancer_ex <- c('age|sis|bro|dad|mom')
    cancer_cols <- cancer_cols[-grep(cancer_ex,cancer_cols,ignore.case = T)]
    exclusion_cols <- c('he_m092_ovarian_cysts','he_m089_endometriosis','he_m090_uterine_polyps')
    f_ex <- f_he[c(exclusion_cols)] %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- f_he$epr_number
    ex_df <- cbind(epr_number,f_ex) %>% dplyr::select(epr_number,sum)
    he_vars <- f_he %>% dplyr::select(epr_number,he_m091_uterine_tumors,he_age_derived)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Fibroids = case_when(he_m091_uterine_tumors == 1 ~ 1,
                                  he_m091_uterine_tumors == 0 & Flag == 0 ~ 0,
                                  TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Fibroids)
    return(pheno_data)
    
  }
  
  .boneloss_phenotype <- function(pegs_data){
    pegs_data[] <- lapply(pegs_data, function(x) as.numeric(as.character(x))) 
    bone_cancer <- pegs_data %>% dplyr::select(epr_number,he_o106_cancer_bone_PARQ_CHILDQ)
    bone_cancer[is.na(bone_cancer)] <- 0
    ex_df <- pegs_data %>% 
      dplyr::filter(!is.na(he_j062_bone_loss) & !is.na(he_j063_osteoporosis) ) %>% 
      dplyr::select(epr_number,he_j062_bone_loss,he_j063_osteoporosis) %>% 
      column_to_rownames('epr_number') %>% 
      dplyr::mutate(sum = rowSums(.,na.rm = T),
                    Bone_Loss = ifelse(sum >0,1,0)) %>% 
      rownames_to_column('epr_number') %>% 
      mutate(epr_number = as.numeric(as.character(epr_number))) %>% 
      left_join(.,bone_cancer, by = "epr_number") %>% 
      mutate(Bone_Loss = ifelse(he_o106_cancer_bone_PARQ_CHILDQ == 1 & Bone_Loss == 0,NA,Bone_Loss) ) %>% 
      dplyr::select(epr_number,Bone_Loss)
    pheno_data <- ex_df 
    return(pheno_data)
    
    
  }
  
  .migraines_phenotype <- function(pegs_data){
    pegs_data[] <- lapply(pegs_data, function(x) as.numeric(as.character(x))) 
    cancer <- pegs_data %>% dplyr::select(epr_number,he_o107_cancer_brain_PARQ_CHILDQ)
    cancer[is.na(cancer)] <- 0
    ex_df <- pegs_data %>% 
      dplyr::select(epr_number,he_e032_migraine) %>% 
      dplyr::mutate(Migraines = ifelse(he_e032_migraine == 1,1,0)) %>% 
      mutate(epr_number = as.numeric(as.character(epr_number))) %>% 
      left_join(.,cancer, by = "epr_number") %>% 
      mutate(Migraines = ifelse(he_o107_cancer_brain_PARQ_CHILDQ == 1 & Migraines == 0,NA,Migraines) ) %>% 
      dplyr::select(epr_number,Migraines)
    pheno_data <- ex_df
    return(pheno_data)
    
  }
  
  .ida_phenotype <- function(pegs_data){
    exclusion_cols <- c('he_i060_pernicious_anemia','he_i061_sickle_cell','he_o115_cancer_leukemia_PARQ_CHILDQ')
    pegs_data[] <- lapply(pegs_data, function(x) as.numeric(as.character(x))) 
    d <- pegs_data[c(exclusion_cols)] %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- pegs_data$epr_number
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- pegs_data %>% dplyr::select(epr_number,he_i059_iron_anemia)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Iron_Def_Anemia = case_when(he_i059_iron_anemia == 1 ~ 1,
                                         he_i059_iron_anemia == 0 & Flag == 0 ~ 0,
                                         TRUE ~ NA_real_)) %>%  
      dplyr::select(epr_number,Iron_Def_Anemia)
    return(pheno_data)
    
  }
  
  .ovarian_cysts_phenotype <- function(pegs_data){ 
    f_he <- pegs_data %>% dplyr::filter(`_he_gender_`==1)
    f_he[] <- lapply(f_he, function(x) as.numeric(as.character(x))) 
    cancer_cols <- names(pegs_data)[grep('cancer',names(pegs_data),ignore.case = T)]
    cancer_ex <- c('age|sis|bro|dad|mom')
    cancer_cols <- cancer_cols[-grep(cancer_ex,cancer_cols,ignore.case = T)]
    exclusion_cols <- c('he_m091_uterine_tumors','he_m089_endometriosis','he_m090_uterine_polyps')
    f_ex <- f_he[c(exclusion_cols)] %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- f_he$epr_number
    ex_df <- cbind(epr_number,f_ex) %>% dplyr::select(epr_number,sum)
    he_vars <- f_he %>% dplyr::select(epr_number,he_m092_ovarian_cysts)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Ovarian_Cysts = case_when(he_m092_ovarian_cysts == 1 ~ 1,
                                       he_m092_ovarian_cysts == 0 & Flag == 0 ~ 0,
                                       TRUE ~ NA_real_)) %>%  
      dplyr::select(epr_number,Ovarian_Cysts)
    return(pheno_data)
    
    
  }
  
  .asthma_phenotype <- function(pegs_data){
    exclusion_cols <- c('he_d025_copd','he_d026_ipf','he_d027_tb_PARQ',"he_d028_cough_breathlessness","he_d029_chest_wheeze")
    ex_df <- pegs_data[c(exclusion_cols)] 
    ex_df[] <- lapply(ex_df, function(x) as.numeric(as.character(x))) 
    d <- ex_df %>% 
      as.data.frame(.) %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- pegs_data$epr_number
    
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- pegs_data %>% dplyr::select(epr_number,he_d030_asthma_PARQ)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Asthma = case_when(he_d030_asthma_PARQ == 1 ~ 1,
                                he_d030_asthma_PARQ == 0 & Flag == 0 ~ 0,
                                TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Asthma)
    return(pheno_data)
    
    
  }
  
  .t2d_phenotype <- function(pegs_data){
    exclusion_cols <- c('he_c022a_diabetes_preg_CHILDQ')
    ex_df <- pegs_data[c(exclusion_cols)] 
    ex_df[] <- lapply(ex_df, function(x) as.numeric(as.character(x))) 
    d <- ex_df %>% 
      as.data.frame(.) %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- pegs_data$epr_number
    
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- pegs_data %>% dplyr::select(epr_number,he_c022_diabetes_PARQ)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Type_2_Diabetes = case_when(he_c022_diabetes_PARQ == 1 ~ 1,
                                         he_c022_diabetes_PARQ == 0 & Flag == 0 ~ 0,
                                         TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Type_2_Diabetes)
    return(pheno_data)
    
    
  }
  
  .ar_phenotype <- function(pegs_data){
    pheno_data <- pegs_data %>% 
      mutate("Allergic_Rhinitis" = ifelse(he_d024_allergies==1,1,0)) %>% 
      dplyr::select(epr_number,Allergic_Rhinitis)
    return(pheno_data)
    
  }
  
  .cvd_phenotype <- function(pegs_data){
    tokeep1 <- c("epr_number","he_b007_hypertension_PARQ","he_b008_high_cholesterol",
                 "he_b009_atherosclerosis", "he_b010_cardiac_arrhythmia","he_b011_angina",
                 "he_b012_heart_attack","he_b013_coronary_artery","he_b014_congestive_heart_failure",
                 "he_b015_poor_blood_flow","he_b017_blood_clots",
                 "he_b018_angioplasty","he_b020_stroke_PARQ", "he_t203_income", "he_s179_100_cigarettes_PARQ")
    CVD.patients <- pegs_data %>% dplyr::select(all_of(tokeep1))
    CVD_Outcomes <- CVD.patients %>% dplyr::mutate(Stroke = ifelse(he_b020_stroke_PARQ== 1,1,0),
                                                   Heart_Attack = ifelse(he_b012_heart_attack== 1,1,0),
                                                   Arrhythmia = ifelse(he_b010_cardiac_arrhythmia== 1,1,0),
                                                   CAD = ifelse(he_b013_coronary_artery== 1,1,0),
                                                   CHF = ifelse(he_b014_congestive_heart_failure== 1,1,0),
                                                   Cholesterol = ifelse(he_b008_high_cholesterol== 1,1,0),
                                                   Hypertension = ifelse(he_b007_hypertension_PARQ== 1,1,0),
                                                   AVSD = ifelse(he_b009_atherosclerosis== 1,1,0),
                                                   Angina = ifelse(he_b011_angina== 1,1,0),
                                                   Angioplasty = ifelse(he_b018_angioplasty== 1,1,0),
                                                   PBF = ifelse(he_b015_poor_blood_flow== 1,1,0),
                                                   Bclots = ifelse(he_b017_blood_clots== 1,1,0)
    ) %>% 
      dplyr::select(epr_number,Stroke,Heart_Attack,Arrhythmia,CAD,CHF,Cholesterol,Hypertension,AVSD,Angina,Angioplasty,PBF,Bclots)
    CVD_Outcomes$Atherogenic <- NA
    CVD_Outcomes$Atherogenic[CVD.patients$he_b020_stroke_PARQ == 1 |
                               CVD.patients$he_b012_heart_attack == 1 |
                               CVD.patients$he_b013_coronary_artery == 1 |
                               CVD.patients$he_b011_angina == 1 |
                               CVD.patients$he_b018_angioplasty == 1 |
                               CVD.patients$he_b009_atherosclerosis == 1] <- 1
    pheno_data <- CVD_Outcomes %>% mutate(Atherogenic = ifelse(is.na(Atherogenic),0,Atherogenic))
    
    return(pheno_data)
    
  }
  
  if(phenotype == "lower_gi_polyps"){
    
    pheno_data <- .lower_gi_phenotype(pegs_data = pegs_data)
    
  }
  
  if(phenotype == "fibroids"){
    pheno_data <- .fibroids_phenotype(pegs_data = pegs_data)
  }
  
  if(phenotype == "boneloss"){
    pheno_data <- .boneloss_phenotype(pegs_data = pegs_data)
    
  }
  
  if(phenotype == "migraines"){
    pheno_data <- .migraines_phenotype(pegs_data = pegs_data)
    
  }
  
  if(phenotype == "IDA"){
    pheno_data <- .ida_phenotype(pegs_data = pegs_data)
    
  }
  
  if(phenotype == "ovariancysts"){
    pheno_data <- .ovarian_cysts_phenotype(pegs_data = pegs_data)
  }
  
  if(phenotype == "asthma"){
    pheno_data <- .asthma_phenotype(pegs_data = pegs_data)
    
    
  }
  
  if(phenotype == "T2D"){
    
    pheno_data <- .t2d_phenotype(pegs_data = pegs_data)
    
  }
  
  if(phenotype == "Allergic_Rhinitis"){
    
    pheno_data <- .ar_phenotype(pegs_data = pegs_data)
    
  }
  
  if(phenotype == "Stroke"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Stroke)
    
  }
  
  if(phenotype == "Heart_Attack"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Heart_Attack)
    
  }
  
  if(phenotype == "Arrhythmia"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Arrhythmia)
    
  }
  
  if(phenotype == "CAD"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,CAD)
    
  }
  
  if(phenotype == "CHF"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,CHF)
    
  }
  if(phenotype == "Cholesterol"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Cholesterol)
    
  }
  if(phenotype == "Hypertension"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Hypertension)
    
  }
  if(phenotype == "AVSD"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,AVSD)
    
  }
  if(phenotype == "Angina"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Angina)
    
  }
  if(phenotype == "Angioplasty"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Angioplasty)
    
  }
  if(phenotype == "PBF"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,PBF)
    
  }
  if(phenotype == "Bclots"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Bclots)
    
  }
  if(phenotype == "Atherogenic"){
    
    pheno_data <- .cvd_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Atherogenic)
    
  }
  
  pheno_data <- sapply( pheno_data, as.numeric )
  pheno_data <- as.data.frame(pheno_data)
  #names(pheno_data) <- c('epr_number','Y')
  
  pheno_data %>%
    select(-epr_number) %>%
    group_by_all() %>%
    count() %>%
    print(n = Inf)

  # Create cases(=1) & controls(=0) dataframe
  #pegs_pheno_prepd <- pheno_data %>% select(epr_number, Y)
  cat("Phenotype = ", phenotype, "\n")
  print(table(pheno_data[2], exclude = NULL))

  return(pheno_data)
}



