library(ADNIMERGE)
library(tidyverse)
library(transplantr)

adni_condition_scrape <-
  function(condition_name,
           include_words,
           exclude_words,
           excluded_records) {
    find.package("ADNIMERGE")
    find.package("tidyverse")
    
    init_health <- ADNIMERGE::inithealth
    rec_hist <- ADNIMERGE::recmhist
    
    init_health$est_onset_date <-
      as.Date((ifelse(
        is.na(stringr::str_detect(init_health$IHDTONSET, "xx-xx")),
        init_health$USERDATE,
        ifelse((stringr::str_detect(init_health$IHDTONSET, "xx")) ==
                 FALSE,
               as.Date(strptime(init_health$IHDTONSET,
                                format = "%Y-%m-%d")),
               ifelse((stringr::str_detect(init_health$IHDTONSET, "xx-xx")) ==
                        TRUE,
                      lubridate::ceiling_date((as.Date(
                        strptime(init_health$IHDTONSET, format = "%Y")
                      )), unit = "years"),
                      lubridate::ceiling_date(lubridate::ym(init_health$IHDTONSET), unit =
                                                "months")
               )
        )
      )),
      origin = "1970-01-01")
    
    init_health$est_cease_date <-
      as.Date((ifelse(
        is.na(stringr::str_detect(init_health$IHCEASE, "xx-xx")),
        as.Date("2100-12-31"),
        ifelse((stringr::str_detect(init_health$IHCEASE, "xx")) ==
                 FALSE,
               as.Date(strptime(init_health$IHCEASE,
                                format = "%Y-%m-%d")),
               ifelse((stringr::str_detect(init_health$IHCEASE, "xx-xx")) ==
                        TRUE,
                      lubridate::floor_date((as.Date(
                        strptime(init_health$IHCEASE, format = "%Y")
                      )), unit = "years"),
                      lubridate::ceiling_date(lubridate::ym(init_health$IHCEASE), unit =
                                                "months")
               )
        )
      )),
      origin = "1970-01-01")
    
    ## this bit of code is meant to work around an error thrown by stringr::ceiling_date for recent medical history
    
    rec_hist$ceilinged_onset <- lubridate::ymd(lubridate::year(rec_hist$MHDTONSET),
                                               truncated = 2L) + lubridate::years(1)
    
    rec_hist$est_onset_date <-
      as.Date(ifelse(
        is.na(rec_hist$MHDTONSET),
        rec_hist$USERDATE,
        rec_hist$ceilinged_onset
      ),
      origin = "1970-01-01")
    
    include_words_case_insensitive <-
      paste(paste0("(?i)", include_words), collapse = "|")
    init_health$condition <-
      ifelse(
        stringr::str_detect(
          stringr::str_to_lower(init_health$IHDESC),
          include_words_case_insensitive
        ),
        1,
        0
      )
    rec_hist$condition <-
      ifelse(
        stringr::str_detect(
          stringr::str_to_lower(rec_hist$MHDESC),
          include_words_case_insensitive
        ),
        1,
        0
      )
    
    if (!missing(exclude_words)) {
      exclude_words_case_insensitive <-
        paste(paste0("(?i)", exclude_words), collapse = "|")
      init_health$condition <-
        ifelse(
          stringr::str_detect(inithealth$IHDESC, exclude_words_case_insensitive),
          0,
          init_health$condition
        )
      rec_hist$condition <-
        ifelse(
          stringr::str_detect(rec_hist$MHDESC, exclude_words_case_insensitive),
          0,
          rec_hist$condition
        )
    }
    
    init_health_reduced <- init_health %>%
      dplyr::filter(condition == 1) %>%
      dplyr::select(RID, VISCODE, est_onset_date, est_cease_date, IHDESC)
    
    rec_hist$est_cease_date <- as.Date("2100-12-31")
    rec_hist_reduced <- rec_hist %>%
      dplyr::filter(MHCUR == "Yes" & condition == 1) %>%
      dplyr::select(RID, VISCODE, est_onset_date, est_cease_date, MHDESC)
    
    colnames(init_health_reduced) <-
      c("RID", "viscode", "onset", "cease", "desc")
    colnames(rec_hist_reduced) <-
      c("RID", "viscode", "onset", "cease", "desc")
    
    merged_recs <-
      rbind.data.frame(init_health_reduced, rec_hist_reduced)
    merged_dates <-
      merged_recs %>% dplyr::group_by(RID) %>% dplyr::mutate(cease = min(cease))
    merged_dates <-
      merged_dates %>% dplyr::arrange(onset) %>% dplyr::distinct_at(vars(RID), .keep_all =
                                                                      TRUE)
    
    colnames(merged_dates) <-
      c(
        "RID",
        "viscode",
        paste0(condition_name, "_onset"),
        paste0(condition_name, "_cease"),
        paste0(condition_name, "_desc")
      )
    
    return(merged_dates)
  }


adni_condition_merge <-
  function(biomarker_data,
           condition_data,
           condition_name) {
    merged_data <-
      dplyr::left_join(biomarker_data, condition_data, by = "RID")
    merged_data$condition <- as.numeric(
      merged_data %>% dplyr::select(EXAMDATE) >= merged_data %>% dplyr::select(paste0(condition_name, "_onset"))
      &
        merged_data %>% dplyr::select(EXAMDATE) < merged_data %>% dplyr::select(paste0(condition_name, "_cease"))
    )
    merged_data$condition <-
      ifelse(is.na(merged_data$condition), 0, merged_data$condition)
    colnames(merged_data)[which(colnames(merged_data) == "condition")] <-
      condition_name
    return(merged_data)
  }

adni_merge <- ADNIMERGE::adnimerge
adni_merge_demog <- ADNIMERGE::ptdemog %>% dplyr::filter(VISCODE != "f")
adni_lab_data <- ADNIMERGE::labdata %>% dplyr::filter(VISCODE != "f")
bateman_abeta_ratio <- ADNIMERGE::batemanlab
washu_abeta_ratio <- ADNIMERGE::plasma_abeta_project_wash_u
bateman_abeta_ratio <-
  readr::read_delim("~/Downloads/batemanlab_20190621.csv")
adni_web_demog <-
  readr::read_delim("~/Downloads/PTDEMOG.csv") %>% dplyr::filter(VISCODE !=
                                                                   "f")

upenn_mk_12 <-readr::read_delim("~/Downloads/UPENNBIOMK12_01_04_21.csv")
upenn_mk_10 <-readr::read_delim("~/Downloads/UPENNBIOMK10_07_29_19.csv") %>% dplyr::rename(EXAMDATE=DRAWDATE,ABETA=ABETA42)
upenn_mk_9 <-readr::read_delim("~/Downloads/UPENNBIOMK9_04_19_17.csv")


upenn_mk_9 <- upenn_mk_9 %>% dplyr::mutate(ABETA=as.numeric(ABETA),TAU=as.numeric(TAU),PTAU=as.numeric(PTAU))
upenn_mk_9$ABETA <- ifelse(is.na(upenn_mk_9$ABETA),readr::parse_number(upenn_mk_9$COMMENT),upenn_mk_9$ABETA)


adni_mod_hach <- ADNIMERGE::modhach %>% dplyr::filter(VISCODE != "f")
adni_vitals <- ADNIMERGE::vitals %>% dplyr::filter(VISCODE != "f")


## only has ADNI data
upenn_adni_dian <-
  readr::read_delim("~/Downloads/UPENNBIOMKADNIDIAN2017.csv")

upenn_merged_csf_biomarkers <- dplyr::bind_rows(upenn_mk_12,upenn_mk_10,upenn_mk_9)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::select(-VISCODE) %>% dplyr::rename(VISCODE=VISCODE2)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::mutate(ABETA=round(ABETA))
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::group_by(RID,EXAMDATE) %>% dplyr::filter(RUNDATE==max(RUNDATE)) %>% dplyr::ungroup()

bateman_abeta_ratio <-
  bateman_abeta_ratio %>% dplyr::rename(abeta_ratio = RATIO_ABETA42_40_BY_ISTD_TOUSE)
washu_abeta_ratio <-
  washu_abeta_ratio %>% dplyr::rename(abeta_ratio = STANDARDIZED_PLASMAAB4240)
washu_abeta_ratio <-
  washu_abeta_ratio %>% dplyr::rename(abeta_ratio_non_std = PLASMAAB4240)

washu_abeta_ratio$origin <- "WashU"
bateman_abeta_ratio$origin <- "Bateman"

abeta_ratio_merged <-
  rbind.data.frame(
    washu_abeta_ratio %>% dplyr::filter(QC_STATUS == "Passed") %>%
      dplyr::select(RID, VISCODE, EXAMDATE, abeta_ratio, origin),
    bateman_abeta_ratio %>% dplyr::filter(QC_STATUS ==
                                            "Passed", INSTRUMENT == "Lumos", INJECTION == "a") %>%
      dplyr::select(RID, VISCODE, EXAMDATE, abeta_ratio, origin)
  ) %>% dplyr::filter(RID != 999999)

av45 = ADNIMERGE::ucberkeleyav45
fbb = ADNIMERGE::ucberkeleyfbb
av45$AV45.DATE = av45$EXAMDATE
av45 = av45[, c(
  'RID',
  'AV45.DATE',
  'VISCODE',
  "SUMMARYSUVR_WHOLECEREBNORM",
  "SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF"
)]
fbb = ADNIMERGE::ucberkeleyfbb
fbb$FBB.DATE = fbb$EXAMDATE
fbb = fbb[, c(
  'RID',
  'FBB.DATE',
  'VISCODE',
  "SUMMARYSUVR_WHOLECEREBNORM",
  "SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF"
)]
d = merge(av45, fbb, by = c('RID', 'VISCODE'), all = T)
d$Centiloid = 188.22 * d$SUMMARYSUVR_WHOLECEREBNORM.x - 189.16
d$Centiloid[(is.na(d$Centiloid))] = 157.15 * d$SUMMARYSUVR_WHOLECEREBNORM.y[is.na(d$Centiloid)] -
  151.87
d$Cent.DATE = d$AV45.DATE
d$Cent.DATE[is.na(d$Cent.DATE)] = d$FBB.DATE[is.na(d$Cent.DATE)]
d$AmyloidPos = d$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF
d$AmyloidPos[is.na(d$AmyloidPos)] = d$SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF[is.na(d$AmyloidPos)]
d = d[!is.na(d$Centiloid), c('RID', 'Cent.DATE', 'VISCODE', 'Centiloid', 'AmyloidPos')]
d <- d[order(d$RID, d$Cent.DATE), ]
d <- d[!duplicated(d[c("RID", "Cent.DATE")]), ]
d$Cent.DATE <- as.Date(d$Cent.DATE)
d$EXAMDATE <- as.Date(d$Cent.DATE)
amyloid_pet <- d

adni_merge_demog_uniques <- adni_merge_demog %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(EXAMDATE == min(EXAMDATE))

adni_web_demog_uniques <- adni_web_demog %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(USERDATE == min(USERDATE))

adni_merge_demog_uniques <-
  dplyr::left_join(
    adni_merge_demog_uniques,
    adni_web_demog_uniques %>%
      dplyr::select(PTDOBYY, PTDOBMM, RID),
    by = "RID"
  )

adni_merge_demog_uniques_reduced <-
  adni_merge_demog_uniques %>% dplyr::select(RID,
                                             ORIGPROT,
                                             PTGENDER,
                                             PTETHCAT,
                                             PTRACCAT,
                                             PTDOBYY,
                                             PTDOBMM,
                                             PTEDUCAT)


adni_lab_data_demog_reduced <-
  dplyr::left_join(
    adni_lab_data,
    adni_merge_demog_uniques,
    by = "RID",
    suffix = c(".labs", ".demos")
  )

adni_lab_data_demog_reduced <-
  adni_lab_data_demog_reduced %>% mutate(
    glucose = as.numeric(RCT11),
    creatinine = as.numeric(RCT392),
    cholesterol = as.numeric(RCT20)
  )

adni_lab_data_demog_reduced$age_at_lab <-
  round(as.numeric(lubridate::year(
    adni_lab_data_demog_reduced$EXAMDATE.labs
  )) -
    adni_lab_data_demog_reduced$PTDOBYY +
    ((
      as.numeric(
        lubridate::month(adni_lab_data_demog_reduced$EXAMDATE.labs)
      ) - adni_lab_data_demog_reduced$PTDOBMM
    ) / 12),
  1)

adni_lab_data_demog_reduced$race_for_egfr <-
  ifelse(
    adni_lab_data_demog_reduced$PTRACCAT == "Black or African American",
    "black",
    "non-black"
  )

adni_lab_data_demog_reduced$gender_for_egfr <-
  ifelse(
    adni_lab_data_demog_reduced$PTGENDER == "Male",
    "M",
    "F"
  )

adni_lab_data_demog_reduced$eGFR <-
  transplantr::ckd_epi(
    adni_lab_data_demog_reduced$creatinine,
    adni_lab_data_demog_reduced$age_at_lab,
    adni_lab_data_demog_reduced$gender_for_egfr,
    adni_lab_data_demog_reduced$race_for_egfr,
    units = "US"
  )

adni_lab_data_demog_reduced <- adni_lab_data_demog_reduced %>% dplyr::mutate(eGFR_cat = with(.,case_when(
  (eGFR>90)~"Greater than 90",
  (eGFR>60 & eGFR<=90)~"61-90",
  (eGFR>45 & eGFR <=60)~"45-60",
  (eGFR<=45)~"Less than 45",
  (is.na(eGFR))~"No eGFR Available")))

adni_lab_data_demog_reduced <- adni_lab_data_demog_reduced %>%
  mutate(., diabetes_cat = with(
    .,
    case_when(
      (glucose >= 125) ~ 'Diabetes',
      (glucose >= 100 & glucose < 125) ~ 'Impaired Glycemia',
      (glucose < 100) ~ 'Normoglycemic',
      is.na(glucose) ~ 'Failed Draw',
    )
  )) %>%
  mutate(., ckd_cat = with(
    .,
    case_when(
      (eGFR <= 60) ~ 'Suspected CKD',
      (eGFR > 60) ~ 'No Suspected CKD',
      is.na(eGFR) ~ 'Could Not Calculate eGFR',
    )
  )) %>%
  mutate(., dyslipidemia_cat = with(
    .,
    case_when(
      (cholesterol >= 200) ~ 'Hyperlipidemia',
      (cholesterol < 200) ~ 'Normal Cholesterol',
      is.na(cholesterol) ~ 'Failed Draw',
    )
  )) %>%
  dplyr::select(
    RID,
    age_at_lab,
    glucose,
    diabetes_cat,
    creatinine,
    eGFR,
    ckd_cat,
    cholesterol,
    dyslipidemia_cat,
    race_for_egfr,
    gender_for_egfr,
    eGFR_cat
  )

adni_lab_data_demog_reduced <- adni_lab_data_demog_reduced %>%
  dplyr::distinct_at(vars(RID), .keep_all = TRUE)

adni_mod_hach$htn_hach <-
  ifelse(adni_mod_hach$HMHYPERT == "Present - 1 point", 1, 0)

adni_vitals <-
  adni_vitals %>% dplyr::mutate(
    diastolic_bp = as.numeric(VSBPDIA),
    systolic_bp = as.numeric(VSBPSYS),
    vitals_date = USERDATE
  )
adni_vitals$htn_vitals <-
  ifelse(adni_vitals$diastolic_bp >= 90 |
           adni_vitals$systolic_bp >= 140,
         1,
         0)

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

csf_all_conditions <-
  adni_condition_merge(upenn_merged_csf_biomarkers, adni_diabetes, "diabetes")
csf_all_conditions <-
  adni_condition_merge(csf_all_conditions, adni_ckd, "ckd")
csf_all_conditions <-
  adni_condition_merge(csf_all_conditions, adni_dyslipidemia, "dyslipidemia")
csf_all_conditions <-
  adni_condition_merge(csf_all_conditions, adni_hypertension, "hypertension")
csf_all_conditions<-
  adni_condition_merge(csf_all_conditions, adni_stroke, "stroke",is_chronic=FALSE)
csf_all_conditions<-
  adni_condition_merge(csf_all_conditions, adni_heart_attack, "heart_attack",is_chronic=FALSE)
csf_all_conditions <-
  dplyr::left_join(csf_all_conditions, adni_lab_data_demog_reduced, by = "RID",suffix=c(".CSF",".labs"))
csf_all_conditions <-
  dplyr::left_join(csf_all_conditions,
                   adni_mod_hach %>% dplyr::select(RID, htn_hach),
                   by = "RID",suffix=c(".CSF_labs",".hach"))
csf_all_conditions <-
  dplyr::left_join(
    csf_all_conditions,
    adni_vitals %>% dplyr::select(RID, diastolic_bp, systolic_bp, htn_vitals, vitals_date),
    by = "RID", suffix=c(".CSF_labs_hach",".vitals")
  )
csf_all_conditions$date_diff <-
  ifelse(
    csf_all_conditions$EXAMDATE - csf_all_conditions$vitals_date < 0,
    999999,
    csf_all_conditions$EXAMDATE - csf_all_conditions$vitals_date
  )
csf_all_conditions <- csf_all_conditions %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(date_diff == min(date_diff))
csf_all_conditions$diabetes_joined <-
  ifelse(csf_all_conditions$diabetes == 1 |
           csf_all_conditions$diabetes_cat == "Diabetes",
         1,
         0)
csf_all_conditions$ckd_joined <-
  ifelse(csf_all_conditions$ckd == 1 |
           csf_all_conditions$ckd_cat == "Suspected CKD",
         1,
         0)
csf_all_conditions$dyslipidemia_joined <-
  ifelse(
    csf_all_conditions$dyslipidemia == 1 |
      csf_all_conditions$dyslipidemia_cat == "Hyperlipidemia",
    1,
    0
  )
csf_all_conditions$hypertension_joined <-
  ifelse(
    csf_all_conditions$hypertension == 1 |
      (
        csf_all_conditions$htn_vitals == 1 &
          csf_all_conditions$date_diff != 9999
      ) | csf_all_conditions$htn_hach == 1,
    1,
    0
  )
csf_all_conditions <-
  dplyr::left_join(csf_all_conditions, adni_merge_demog_uniques_reduced, by =
                     "RID")
csf_all_conditions$age_at_exam <-
  round(as.numeric(lubridate::year(csf_all_conditions$EXAMDATE)) -
          csf_all_conditions$PTDOBYY +
          ((
            as.numeric(lubridate::month(csf_all_conditions$EXAMDATE)) - csf_all_conditions$PTDOBMM
          ) / 12),
        1)

abeta_all_conditions <-
  adni_condition_merge(abeta_ratio_merged, adni_diabetes, "diabetes")
abeta_all_conditions <-
  adni_condition_merge(abeta_all_conditions, adni_ckd, "ckd")
abeta_all_conditions <-
  adni_condition_merge(abeta_all_conditions, adni_dyslipidemia, "dyslipidemia")
abeta_all_conditions <-
  adni_condition_merge(abeta_all_conditions, adni_hypertension, "hypertension")
abeta_all_conditions <-
  dplyr::left_join(abeta_all_conditions, adni_lab_data_demog_reduced, by =
                     "RID")
abeta_all_conditions <-
  dplyr::left_join(abeta_all_conditions,
                   adni_mod_hach %>% dplyr::select(RID, htn_hach),
                   by = "RID")
abeta_all_conditions <-
  dplyr::left_join(
    abeta_all_conditions,
    adni_vitals %>% dplyr::select(RID, diastolic_bp, systolic_bp, htn_vitals, vitals_date),
    by = "RID"
  )
abeta_all_conditions$date_diff <-
  ifelse(
    abeta_all_conditions$EXAMDATE - abeta_all_conditions$vitals_date < 0,
    999999,
    abeta_all_conditions$EXAMDATE - abeta_all_conditions$vitals_date
  )
abeta_all_conditions <- abeta_all_conditions %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(date_diff == min(date_diff))
abeta_all_conditions$diabetes_joined <-
  ifelse(
    abeta_all_conditions$diabetes == 1 |
      abeta_all_conditions$diabetes_cat == "Diabetes",
    1,
    0
  )
abeta_all_conditions$ckd_joined <-
  ifelse(abeta_all_conditions$ckd == 1 |
           abeta_all_conditions$ckd_cat == "Suspected CKD",
         1,
         0)
abeta_all_conditions$dyslipidemia_joined <-
  ifelse(
    abeta_all_conditions$dyslipidemia == 1 |
      abeta_all_conditions$dyslipidemia_cat == "Hyperlipidemia",
    1,
    0
  )
abeta_all_conditions$hypertension_joined <-
  ifelse(
    abeta_all_conditions$hypertension == 1 |
      (
        abeta_all_conditions$htn_vitals == 1 &
          abeta_all_conditions$date_diff != 9999
      ) | abeta_all_conditions$htn_hach == 1,
    1,
    0
  )
abeta_all_conditions <-
  dplyr::left_join(abeta_all_conditions, adni_merge_demog_uniques_reduced, by =
                     "RID")
abeta_all_conditions$age_at_exam <-
  round(as.numeric(lubridate::year(abeta_all_conditions$EXAMDATE)) -
          abeta_all_conditions$PTDOBYY +
          ((
            as.numeric(lubridate::month(abeta_all_conditions$EXAMDATE)) - abeta_all_conditions$PTDOBMM
          ) / 12),
        1)

pet_all_conditions <-
  adni_condition_merge(amyloid_pet, adni_diabetes, "diabetes")
pet_all_conditions <-
  adni_condition_merge(pet_all_conditions, adni_ckd, "ckd")
pet_all_conditions <-
  adni_condition_merge(pet_all_conditions, adni_dyslipidemia, "dyslipidemia")
pet_all_conditions <-
  adni_condition_merge(pet_all_conditions, adni_hypertension, "hypertension")
pet_all_conditions <-
  dplyr::left_join(pet_all_conditions, adni_lab_data_demog_reduced, by = "RID")
pet_all_conditions <-
  dplyr::left_join(pet_all_conditions,
                   adni_mod_hach %>% dplyr::select(RID, htn_hach),
                   by = "RID")
pet_all_conditions <-
  dplyr::left_join(
    pet_all_conditions,
    adni_vitals %>% dplyr::select(RID, diastolic_bp, systolic_bp, htn_vitals, vitals_date),
    by = "RID"
  )
pet_all_conditions$date_diff <-
  ifelse(
    pet_all_conditions$EXAMDATE - pet_all_conditions$vitals_date < 0,
    999999,
    pet_all_conditions$EXAMDATE - pet_all_conditions$vitals_date
  )
pet_all_conditions <- pet_all_conditions %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(date_diff == min(date_diff))
pet_all_conditions$diabetes_joined <-
  ifelse(pet_all_conditions$diabetes == 1 |
           pet_all_conditions$diabetes_cat == "Diabetes",
         1,
         0)
pet_all_conditions$ckd_joined <-
  ifelse(pet_all_conditions$ckd == 1 |
           pet_all_conditions$ckd_cat == "Suspected CKD",
         1,
         0)
pet_all_conditions$dyslipidemia_joined <-
  ifelse(
    pet_all_conditions$dyslipidemia == 1 |
      pet_all_conditions$dyslipidemia_cat == "Hyperlipidemia",
    1,
    0
  )
pet_all_conditions$hypertension_joined <-
  ifelse(
    pet_all_conditions$hypertension == 1 |
      (
        pet_all_conditions$htn_vitals == 1 &
          pet_all_conditions$date_diff != 9999
      ) | pet_all_conditions$htn_hach == 1,
    1,
    0
  )
pet_all_conditions <-
  dplyr::left_join(pet_all_conditions, adni_merge_demog_uniques_reduced, by =
                     "RID")
pet_all_conditions$age_at_exam <-
  round(as.numeric(lubridate::year(pet_all_conditions$EXAMDATE)) -
          pet_all_conditions$PTDOBYY +
          ((
            as.numeric(lubridate::month(pet_all_conditions$EXAMDATE)) - pet_all_conditions$PTDOBMM
          ) / 12),
        1)
