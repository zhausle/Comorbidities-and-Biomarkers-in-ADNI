library(ADNIMERGE)
library(tidyverse)
library(transplantr)

adni_diabetes <-
  adni_condition_scrape("diabetes", c("diabet"), c("pre", "borderline", "related"))
adni_dyslipidemia <-
  adni_condition_scrape("dyslipidemia", c("cholester", "lipid"))
adni_ckd <-
  adni_condition_scrape("ckd", c("chronic kidney", "ckd"), "stones")
adni_hypertension <-
  adni_condition_scrape(
    "hypertension",
    c("hypertension", "htn"),
    c("orthostatic", "white coat",
      "borderline", "ocular")
  )
adni_stroke<-
  adni_condition_scrape("stroke",c("stroke","CVA","cerebrovascular accident","brain attack"))
adni_heart_attack<-
  adni_condition_scrape("heart_attack",c("heart attack","myocardial infarction"))

list_of_conditions<-list(adni_diabetes,adni_dyslipidemia,adni_ckd,adni_hypertension,adni_stroke,adni_heart_attack)
condition_is_chronic<-c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE)

csf_all_conditions<-all_conditions_merger(upenn_merged_csf_biomarkers,list_of_conditions,condition_is_chronic)
csf_all_conditions<-other_df_prep(csf_all_conditions)

plasma_all_conditions<-all_conditions_merger(plasma_merged,list_of_conditions,condition_is_chronic)
plasma_all_conditions<-other_df_prep(plasma_all_conditions)

pet_all_conditions<-all_conditions_merger(amyloid_pet,list_of_conditions,condition_is_chronic)
pet_all_conditions<-other_df_prep(pet_all_conditions)

plasma_ptau_all_conditions<-all_conditions_merger(plasma_ptau,list_of_conditions,condition_is_chronic)
plasma_ptau_all_conditions<-other_df_prep(plasma_ptau_all_conditions)

plasma_nfl_all_conditions<-all_conditions_merger(plasma_nfl,list_of_conditions,condition_is_chronic)
plasma_nfl_all_conditions<-other_df_prep(plasma_nfl_all_conditions)
