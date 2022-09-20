library(tidyverse)
load('/Volumes/PEGS/StatGen/lloyddt/Phenotyping/clean_he.Rdata')



## PEGS Phenotyping function. Input data is the health and exposure survey that has been fed through the pegs_convert_type function
## Outputs one phenotype at a time from the accepted phenotype values
## Excludes people from the controls not the cases when they meet the exclusion criteria
## More Phenotypes will be added to the function 

pegs_phenotyping <- function(input_data, 
                             phenotype = c("lower_gi_polyps",'fibroids',
                                           'boneloss','migraines','IDA',
                                           'ovariancysts','asthma','cvd',"T2D","Allergic_Rhinitis")){
  
  if(phenotype == "lower_gi_polyps"){
    
    cancer_cols <- names(input_data)[grep('cancer',names(input_data),ignore.case = T)]
    cancer_ex <- c('age|sis|bro|dad|mom')
    cancer_cols <- cancer_cols[-grep(cancer_ex,cancer_cols,ignore.case = T)]
    gi_cols <- c('he_f038_lactose_intolerance','he_f039_crohns','he_f040_ulcerative_colitis')
    exclusion_cols <- c(cancer_cols,gi_cols)
    ex_df <- input_data[c(exclusion_cols)] 
    ex_df[] <- lapply(ex_df, function(x) as.numeric(as.character(x))) 
    d <- ex_df %>% 
      as.data.frame(.) %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- input_data$epr_number
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- input_data %>% dplyr::select(epr_number,he_f041_polyps,he_age_derived)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Lower_GI_Polyps = case_when(he_f041_polyps == 1 ~ 1,
                                         he_f041_polyps == 0 & Flag == 0 ~ 0,
                                         TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Lower_GI_Polyps)
  }
  
  if(phenotype == "fibroids"){
    f_he <- input_data %>% dplyr::filter(`_he_gender_`==1)
    f_he[] <- lapply(f_he, function(x) as.numeric(as.character(x))) 
    cancer_cols <- names(input_data)[grep('cancer',names(input_data),ignore.case = T)]
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
  }
  
  if(phenotype == "boneloss"){
    input_data[] <- lapply(input_data, function(x) as.numeric(as.character(x))) 
    bone_cancer <- input_data %>% dplyr::select(epr_number,he_o106_cancer_bone_PARQ_CHILDQ)
    bone_cancer[is.na(bone_cancer)] <- 0
    ex_df <- input_data %>% 
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
    
  }
  
  if(phenotype == "migraines"){
    input_data[] <- lapply(input_data, function(x) as.numeric(as.character(x))) 
    cancer <- input_data %>% dplyr::select(epr_number,he_o107_cancer_brain_PARQ_CHILDQ)
    cancer[is.na(cancer)] <- 0
    ex_df <- input_data %>% 
      dplyr::select(epr_number,he_e032_migraine) %>% 
      dplyr::mutate(Migraines = ifelse(he_e032_migraine == 1,1,0)) %>% 
      mutate(epr_number = as.numeric(as.character(epr_number))) %>% 
      left_join(.,cancer, by = "epr_number") %>% 
      mutate(Migraines = ifelse(he_o107_cancer_brain_PARQ_CHILDQ == 1 & Migraines == 0,NA,Migraines) ) %>% 
      dplyr::select(epr_number,Migraines)
    pheno_data <- ex_df
  }
  
  if(phenotype == "IDA"){
    exclusion_cols <- c('he_i060_pernicious_anemia','he_i061_sickle_cell','he_o115_cancer_leukemia_PARQ_CHILDQ')
    input_data[] <- lapply(input_data, function(x) as.numeric(as.character(x))) 
    d <- input_data[c(exclusion_cols)] %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- input_data$epr_number
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- input_data %>% dplyr::select(epr_number,he_i059_iron_anemia)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Iron_Def_Anemia = case_when(he_i059_iron_anemia == 1 ~ 1,
                                         he_i059_iron_anemia == 0 & Flag == 0 ~ 0,
                                         TRUE ~ NA_real_)) %>%  
      dplyr::select(epr_number,Iron_Def_Anemia)
    
  }
  
  if(phenotype == "ovariancysts"){
    f_he <- input_data %>% dplyr::filter(`_he_gender_`==1)
    f_he[] <- lapply(f_he, function(x) as.numeric(as.character(x))) 
    cancer_cols <- names(input_data)[grep('cancer',names(input_data),ignore.case = T)]
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
    
    
  }
  
  if(phenotype == "asthma"){
    
    exclusion_cols <- c('he_d025_copd','he_d026_ipf','he_d027_tb_PARQ',"he_d028_cough_breathlessness","he_d029_chest_wheeze")
    ex_df <- input_data[c(exclusion_cols)] 
    ex_df[] <- lapply(ex_df, function(x) as.numeric(as.character(x))) 
    d <- ex_df %>% 
      as.data.frame(.) %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- input_data$epr_number
    
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- input_data %>% dplyr::select(epr_number,he_d030_asthma_PARQ)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Asthma = case_when(he_d030_asthma_PARQ == 1 ~ 1,
                                he_d030_asthma_PARQ == 0 & Flag == 0 ~ 0,
                                TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Asthma)
    
    
  }
  
  if(phenotype == "cvd"){
    tokeep1 <- c("epr_number","he_b007_hypertension_PARQ","he_b008_high_cholesterol",
                 "he_b009_atherosclerosis", "he_b010_cardiac_arrhythmia","he_b011_angina",
                 "he_b012_heart_attack","he_b013_coronary_artery","he_b014_congestive_heart_failure",
                 "he_b015_poor_blood_flow","he_b017_blood_clots",
                 "he_b018_angioplasty","he_b020_stroke_PARQ", "he_t203_income", "he_s179_100_cigarettes_PARQ")
    CVD.patients <- input_data %>% dplyr::select(all_of(tokeep1))
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
    
    
    
  }
  
  
  if(phenotype == "T2D"){
    
    exclusion_cols <- c('he_c022a_diabetes_preg_CHILDQ')
    ex_df <- input_data[c(exclusion_cols)] 
    ex_df[] <- lapply(ex_df, function(x) as.numeric(as.character(x))) 
    d <- ex_df %>% 
      as.data.frame(.) %>% 
      mutate(sum = rowSums(.,na.rm = T)) %>% dplyr::select(sum,everything())
    epr_number <- input_data$epr_number
    
    ex_df <- cbind(epr_number,d) %>% dplyr::select(epr_number,sum)
    he_vars <- input_data %>% dplyr::select(epr_number,he_c022_diabetes_PARQ)
    pheno_data <- left_join(he_vars,ex_df, by = 'epr_number') %>% 
      mutate(Flag = ifelse(sum == 0,0,1)) %>% 
      mutate(Type_2_Diabetes = case_when(he_c022_diabetes_PARQ == 1 ~ 1,
                                he_c022_diabetes_PARQ == 0 & Flag == 0 ~ 0,
                                TRUE ~ NA_real_)) %>% 
      dplyr::select(epr_number,Type_2_Diabetes)
    
    
    
    
  }
  
  if(phenotype == "Allergic_Rhinitis"){
    
    pheno_data <- input_data %>% 
      mutate("Allergic_Rhinitis" = ifelse(he_d024_allergies==1,1,0)) %>% 
      dplyr::select(epr_number,Allergic_Rhinitis)
    
  }
  
  pheno_data <- sapply( pheno_data, as.numeric )
  pheno_data <- as.data.frame(pheno_data)
  return(pheno_data)
  
}
