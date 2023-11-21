#ExWAS Phenotype Creation

#################
#Fucntions
##############

.lower_gi_phenotype <- function(pegs_data){
  
  cancer_qs <- names(pegs_data)[grep('cancer',names(pegs_data),ignore.case = T)]
  cancer_ex <- c('age|sis|bro|dad|mom|_child_')
  cancer_qs <- cancer_qs[-grep(cancer_ex,cancer_qs,ignore.case = T)]
  cancer_qs <- setdiff(cancer_qs, c("he_o127_cancer_skin_nonmel_PARQ_CHILDQ", "he_o103_cancer_PARQ"))
  
  # Create an indicator variable for any cancer except non-melanoma skin cancer
  sel_cancer_df <- pegs_data %>% select(all_of(cancer_qs))
  sel_cancer_df[] <- lapply(sel_cancer_df, function(x) as.numeric(as.character(x))) 
  sel_cancer_df <- sel_cancer_df %>% 
    mutate(csum = rowSums(.[cancer_qs], na.rm = TRUE)) %>% 
    mutate(non_nmc = ifelse(csum > 0, 1, 0))
  
  polyp_excls <- c(cancer_qs, "he_f038_lactose_intolerance", "he_f039_crohns", "he_f040_ulcerative_colitis")
  polyp_excl_df <- pegs_data %>% select(all_of(polyp_excls))
  polyp_excl_df[] <- lapply(polyp_excl_df, function(x) as.numeric(as.character(x))) 
  polyp_excl_df <- polyp_excl_df %>% 
    mutate(exsum = rowSums(.[polyp_excls], na.rm = TRUE)) %>% 
    mutate(exclude_flag = ifelse(exsum > 0, 1, 0))
  
  epr_number <- pegs_data$epr_number  
  exl_df <- bind_cols("epr_number"=epr_number,polyp_excl_df)
  
  pheno <- pegs_data %>% dplyr::select(epr_number,he_f041_polyps)
  
  pheno_data <- left_join(pheno,exl_df, by = 'epr_number') %>% 
    bind_cols(sel_cancer_df["non_nmc"]) %>% 
    mutate(Lower_GI_Polyps = case_when(he_f041_polyps == 1 ~ 1,
                                       he_f041_polyps == 0 & exclude_flag == 0 & non_nmc == 0 ~ 0,
                                       TRUE ~ NA_real_))
  
  
  return(pheno_data)
  
}
.fibroids_phenotype <- function(pegs_data){
  
  f_he <- pegs_data %>% dplyr::filter(`_he_gender_`==1)
  f_he[] <- lapply(f_he, function(x) as.numeric(as.character(x))) 
  exclusion_cols <- c('he_m092_ovarian_cysts','he_m089_endometriosis','he_m090_uterine_polyps')
  f_ex <- f_he[c(exclusion_cols)] %>% 
    mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
  epr_number <- f_he$epr_number
  ex_df <- cbind(epr_number,f_ex) %>% dplyr::select(epr_number,sum)
  he_vars <- f_he %>% dplyr::select(epr_number,he_m091_uterine_tumors,all_of(exclusion_cols))
  he_vars[] <- lapply(he_vars, function(x) as.numeric(as.character(x))) 
  
  pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
    mutate(Flag = ifelse(sum == 0,0,1)) %>% 
    mutate(Fibroids = case_when(he_m091_uterine_tumors == 1 ~ 1,
                                he_m091_uterine_tumors == 0 & Flag == 0 ~ 0,
                                TRUE ~ NA_real_)) 
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
    dplyr::select(-sum)
  pheno_data <- ex_df 
  return(pheno_data)
  
  
}

.migraines_phenotype <- function(pegs_data){
  pegs_data[] <- lapply(pegs_data, function(x) as.numeric(as.character(x))) 
  cancer <- pegs_data %>% dplyr::select(epr_number,he_o107_cancer_brain_PARQ_CHILDQ)
  cancer[is.na(cancer)] <- 0
  ex_df <- pegs_data %>% 
    dplyr::select(epr_number,he_e032_migraine) %>% 
    #dplyr::mutate(Migraines = ifelse(he_e032_migraine == 1,1,0)) %>% 
    mutate(epr_number = as.numeric(as.character(epr_number))) %>% 
    left_join(.,cancer, by = "epr_number") %>% 
    #mutate(Migraines = ifelse(he_o107_cancer_brain_PARQ_CHILDQ == 1 & Migraines == 0,NA,Migraines) ) %>% 
    mutate(Migraines = case_when(he_e032_migraine == 1 ~ 1,
                                 he_e032_migraine == 0 & he_o107_cancer_brain_PARQ_CHILDQ == 0 ~ 0,
                                 TRUE ~ NA_real_)) 
  pheno_data <- ex_df
  return(pheno_data)
  
}

.ida_phenotype <- function(pegs_data){
  exclusion_cols <- c('he_i060_pernicious_anemia','he_i061_sickle_cell','he_o115_cancer_leukemia_PARQ_CHILDQ')
  pegs_data[] <- lapply(pegs_data, function(x) as.numeric(as.character(x))) 
  d <- pegs_data[c(exclusion_cols)] %>% 
    mutate(sum = rowSums(.,na.rm = T)) 
  epr_number <- pegs_data$epr_number
  ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum,all_of(exclusion_cols))
  he_vars <- pegs_data %>% dplyr::select(epr_number,he_i059_iron_anemia)
  he_vars[] <- lapply(he_vars, function(x) as.numeric(as.character(x))) 
  
  pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
    mutate(Flag = ifelse(sum == 0,0,1)) %>% 
    mutate(Iron_Def_Anemia = case_when(he_i059_iron_anemia == 1 ~ 1,
                                       he_i059_iron_anemia == 0 & Flag == 0 ~ 0,
                                       TRUE ~ NA_real_)) 
  return(pheno_data)
  
}

.ovarian_cysts_phenotype <- function(pegs_data){ 
  f_he <- pegs_data %>% dplyr::filter(`_he_gender_`==1)
  f_he[] <- lapply(f_he, function(x) as.numeric(as.character(x))) 
  exclusion_cols <- c('he_m091_uterine_tumors','he_m089_endometriosis','he_m090_uterine_polyps')
  f_ex <- f_he[c(exclusion_cols)] %>% 
    mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
  epr_number <- f_he$epr_number
  ex_df <- cbind(epr_number,f_ex) %>% dplyr::select(epr_number,sum, all_of(exclusion_cols))
  he_vars <- f_he %>% dplyr::select(epr_number,he_m092_ovarian_cysts)
  he_vars[] <- lapply(he_vars, function(x) as.numeric(as.character(x))) 
  pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
    mutate(Flag = ifelse(sum == 0,0,1)) %>% 
    mutate(Ovarian_Cysts = case_when(he_m092_ovarian_cysts == 1 ~ 1,
                                     he_m092_ovarian_cysts == 0 & Flag == 0 ~ 0,
                                     TRUE ~ NA_real_)) 
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
  
  ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum,all_of(exclusion_cols))
  he_vars <- pegs_data %>% dplyr::select(epr_number,he_d030_asthma_PARQ)
  he_vars[] <- lapply(he_vars, function(x) as.numeric(as.character(x))) 
  pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
    mutate(Flag = ifelse(sum == 0,0,1)) %>% 
    mutate(Asthma = case_when(he_d030_asthma_PARQ == 1 ~ 1,
                              he_d030_asthma_PARQ == 0 & Flag == 0 ~ 0,
                              TRUE ~ NA_real_)) 
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
  
  ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum,all_of(exclusion_cols))
  he_vars <- pegs_data %>% dplyr::select(epr_number,he_c022_diabetes_PARQ,he_c022e_diabetes_age_CHILDQ)
  he_vars[] <- lapply(he_vars, function(x) as.numeric(as.character(x))) 
  pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
    mutate(Flag = ifelse(sum == 0,0,1)) %>% 
    dplyr::mutate(Type_1_Diabetes_Flag = ifelse(he_c022e_diabetes_age_CHILDQ < 20,1,0)) %>% 
    mutate(Type_2_Diabetes = case_when(he_c022_diabetes_PARQ == 1  & Type_1_Diabetes_Flag == 0 & Flag == 0 ~ 1,
                                       he_c022_diabetes_PARQ == 0  ~ 0,
                                       TRUE ~ NA_real_)) 
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
               "he_b018_angioplasty","he_b020_stroke_PARQ", "he_t203_income", "he_s179_100_cigarettes_PARQ","he_b007a_hypertension_preg_CHILDQ","he_b019_stroke_mini")
  CVD.patients <- pegs_data %>% dplyr::select(all_of(tokeep1))
  CVD.patients[] <- lapply(CVD.patients, function(x) as.numeric(as.character(x))) 
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
    #dplyr::mutate(Hypertension = ifelse(he_b007_hypertension_PARQ ==1 & he_b007a_hypertension_preg_CHILDQ == 0,1,0)) %>% 
    dplyr::mutate(Hypertension  = case_when(
      he_b007_hypertension_PARQ == 1 & he_b007a_hypertension_preg_CHILDQ == 0 ~ 1,
      he_b007_hypertension_PARQ == 0 & he_b007a_hypertension_preg_CHILDQ == 0  ~ 0,
      he_b007_hypertension_PARQ == 1 & he_b007a_hypertension_preg_CHILDQ == 1  ~ NA_real_,
      he_b007_hypertension_PARQ == 1  & is.na(he_b007a_hypertension_preg_CHILDQ) ~1,
      is.na(he_b007_hypertension_PARQ) ~ NA_real_,
      TRUE ~ 0)) %>% 
    dplyr::mutate(Stroke = case_when(
      he_b020_stroke_PARQ == 1 | he_b019_stroke_mini == 1 ~1,
      he_b020_stroke_PARQ == 0 | he_b019_stroke_mini == 0 ~0,
      is.na(he_b020_stroke_PARQ) & is.na(he_b019_stroke_mini) ~NA_real_,
      
      TRUE ~ 0
      
    )) %>% 
    #dplyr::mutate(Stroke = ifelse(Stroke == 0 & is.na(he_b019_stroke_mini),0,Stroke)) %>% 
    dplyr::select(epr_number,Stroke,Heart_Attack,Arrhythmia,CAD,CHF,Cholesterol,Hypertension,AVSD,Angina,Angioplasty,
                  PBF,Bclots,he_b009_atherosclerosis,he_b007a_hypertension_preg_CHILDQ,he_b019_stroke_mini,he_b007_hypertension_PARQ,he_b020_stroke_PARQ)
  
  
  
  CVD_Outcomes$Atherogenic <- NA
  CVD_Outcomes$Atherogenic[CVD_Outcomes$Stroke == 1 |
                             CVD.patients$he_b012_heart_attack == 1 |
                             CVD.patients$he_b013_coronary_artery == 1 |
                             CVD.patients$he_b011_angina == 1 |
                             CVD.patients$he_b018_angioplasty == 1 |
                             CVD.patients$he_b009_atherosclerosis == 1] <- 1
  
  
  athero <- CVD.patients %>% 
    dplyr::select(epr_number,he_b020_stroke_PARQ,he_b012_heart_attack,
                  he_b013_coronary_artery,he_b011_angina,he_b018_angioplasty,he_b009_atherosclerosis) %>% 
    mutate(na_count = apply(., MARGIN = 1, function(x) sum(is.na(x)))) %>% 
    dplyr::select(epr_number,na_count)
  
  
  
  pheno_data <- CVD_Outcomes %>% 
    left_join(.,athero, by = "epr_number") %>% 
    mutate(Atherogenic = ifelse(is.na(Atherogenic) & na_count < 6,0,Atherogenic)) %>% 
    dplyr::select(-na_count)
  
  return(pheno_data)
  
}



prepare_pegs_phenotype <- function(pegs_data, 
                                   phenotype = c("lower_gi_polyps",'fibroids',
                                                 'boneloss','migraines','IDA',
                                                 'ovariancysts','asthma',"T2D","Allergic_Rhinitis",
                                                 'Stroke','Heart_Attack','Arrhythmia','CAD','CHF',
                                                 'Cholesterol','Hypertension','AVSD','Angina','Angioplasty',
                                                 'PBF','Bclots','Atherogenic')) {
  pheno_data <- NULL
  
  
  
  if(phenotype == "lower_gi_polyps"){
    
    pheno_data <- .lower_gi_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Lower_GI_Polyps)
    
  }
  
  if(phenotype == "fibroids"){
    pheno_data <- .fibroids_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Fibroids)
  }
  
  if(phenotype == "boneloss"){
    pheno_data <- .boneloss_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Bone_Loss)
    
  }
  
  if(phenotype == "migraines"){
    pheno_data <- .migraines_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Migraines)
    
  }
  
  if(phenotype == "IDA"){
    pheno_data <- .ida_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Iron_Def_Anemia)
    
  }
  
  if(phenotype == "ovariancysts"){
    pheno_data <- .ovarian_cysts_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Ovarian_Cysts)
  }
  
  if(phenotype == "asthma"){
    pheno_data <- .asthma_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Asthma)
    
    
  }
  
  if(phenotype == "T2D"){
    
    pheno_data <- .t2d_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Type_2_Diabetes)
    
  }
  
  if(phenotype == "Allergic_Rhinitis"){
    
    pheno_data <- .ar_phenotype(pegs_data = pegs_data) %>% 
      dplyr::select(epr_number,Allergic_Rhinitis)
    
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
  
  
  return(pheno_data)
}

###################
#Example Usage 
####################


# pheno_codes <- c("lower_gi_polyps",'fibroids',
#                  'boneloss','migraines','IDA',
#                  'ovariancysts','asthma',"T2D",
#                  "Allergic_Rhinitis",
#                  'Cholesterol','Hypertension','Atherogenic')
# final_pheno <- data.frame("epr_number" = clean_he_freeze2$epr_number)
# for(i in seq_along(pheno_codes)){
#   curr_pheno_col <- prepare_pegs_phenotype(pegs_data  = clean_he_freeze2, phenotype = pheno_codes[i])
#   final_pheno <- left_join(final_pheno,curr_pheno_col, by = "epr_number")
#   
#   
#   
# }
# 
# 
# All_ExWAS_Phenotypes <- final_pheno
