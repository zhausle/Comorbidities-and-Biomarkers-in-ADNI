library(ADNIMERGE)
library(tidyverse)
library(transplantr)

## would recommend downloading fresh copies of below files from ADNI 
## as a rule, I try to use .csv instead of ADNImerge because ADNImerge data is labelled and consequently doesn't play well with many functions in R

init_health<-readr::read_delim("~/Projects/Comorbidities Project/INITHEALTH.csv")
rec_hist<- readr::read_delim("~/Projects/Comorbidities Project/RECMHIST.csv")

## creates dataframes with cleaned date formats for initial health and recent medical history

init_health <- init_health %>% dplyr::mutate(
  onset=case_when(
    (is.na(IHDTONSET)) ~ USERDATE,
    (lubridate::month(IHDTONSET)==7 & (lubridate::day(IHDTONSET)==2|lubridate::day(IHDTONSET)==1)) ~ lubridate::ceiling_date(init_health$IHDTONSET,unit="years"), ## Initial health incorrectly has most records as having onset/cease dates of 7/1 or 7/2
    (!is.na(IHDTONSET)) ~ lubridate::ymd(IHDTONSET)
                  ),
  cease=case_when(
      (lubridate::month(IHCEASE)==7 & (lubridate::day(IHCEASE)==2|lubridate::day(IHCEASE)==1)) ~ lubridate::floor_date(init_health$IHCEASE,unit="years"), ## Initial health incorrectly has most records as having onset/cease dates of 7/1 or 7/2
      (!is.na(IHCEASE)) ~ lubridate::ymd(IHCEASE)
)
)

rec_hist$MHDTONSET<-ifelse(rec_hist$MHDTONSET=="x1",NA,rec_hist$MHDTONSET)
rec_hist$MHDTONSET<-ifelse(rec_hist$MHDTONSET=="x4",NA,rec_hist$MHDTONSET)

rec_hist$onset_year<-lubridate::year(strptime(stringr::str_sub(start=-4,rec_hist$MHDTONSET),format="%Y"))
rec_hist$onset_month<-lubridate::month(strptime(paste(stringr::str_sub(end=2,rec_hist$MHDTONSET),"01",stringr::str_sub(start=-4,rec_hist$MHDTONSET),sep="-"),format="%m-%d-%Y"))
rec_hist$onset_day<-lubridate::day(strptime(stringr::str_sub(start=4,end=5,rec_hist$MHDTONSET),format="%d"))

rec_hist <- rec_hist %>% dplyr::mutate(onset=case_when(
    (is.na(onset_year)) ~ USERDATE,
    (is.na(onset_month)) ~ lubridate::ceiling_date(lubridate::ymd(paste(onset_year,"01","01",sep="-")),unit="years"),
    (is.na(onset_day)) ~ lubridate::ceiling_date(lubridate::ymd(paste(onset_year,onset_month,"01",sep="-")),unit="months"),
    (!is.na(onset_day)) ~ lubridate::mdy(MHDTONSET)
    ))

rec_hist_scraper<- rec_hist %>%
  dplyr::filter(MHCUR == 1) %>%
  dplyr::select(RID, VISCODE2, onset, MHDESC,USERDATE) %>%
  dplyr::rename(viscode=VISCODE2,desc=MHDESC)

init_health_scraper <- init_health %>%
  dplyr::select(RID, VISCODE2, onset, cease, IHDESC,USERDATE) %>%
  dplyr::rename(viscode=VISCODE2,desc=IHDESC)

merged_health_scraper<- dplyr::bind_rows(rec_hist_scraper,init_health_scraper)

adni_condition_scrape <-
  function(condition_name,
           include_words,
           exclude_words=NULL) {
    find.package("tidyverse")
    
    include_words_case_insensitive <-
      paste(paste0("(?i)", include_words), collapse = "|")

    merged_health_scraper$condition <-
      ifelse(
        stringr::str_detect(
          stringr::str_to_lower(merged_health_scraper$desc),
          include_words_case_insensitive
        ),
        1,
        0
      )
    
    if (!is.null(exclude_words)) {
      exclude_words_case_insensitive <-
        paste(paste0("(?i)", exclude_words), collapse = "|")
      merged_health_scraper$condition <-
        ifelse(
          stringr::str_detect(merged_health_scraper$desc, exclude_words_case_insensitive),
          0,
          merged_health_scraper$condition
        )
    }
    merged_dates <- merged_health_scraper %>% dplyr::filter(condition==1)
    
merged_dates<-list(df=merged_dates,condition=condition_name)
    return(merged_dates)
  }


adni_condition_merge <-
  function(biomarker_data,
           condition_list) {
    matched_records <-
      dplyr::left_join(biomarker_data, condition_list$df, by = "RID")
    condition_name<-condition_list$condition
    condition_name_enquo<-quo_name(condition_list$condition)
    
    matched_records <-date_comparison(matched_records,EXAMDATE,USERDATE,vars(RID,EXAMDATE),comp_method="before") %>% dplyr::filter(retain_flag==1|alt_flag==1) %>% dplyr::select(RID,EXAMDATE,onset,cease,viscode,desc)

    merged_data <- dplyr::left_join(biomarker_data,matched_records,by=c("RID","EXAMDATE"))
    
    merged_data <- merged_data %>% dplyr::mutate(condition=case_when(
                (is.na(onset)) ~ 0,
                (is.na(cease) & EXAMDATE>=onset) ~ 1,
                (!is.na(cease) & EXAMDATE>=onset & EXAMDATE<cease) ~ 1,
                TRUE ~ 0
    )) %>% dplyr::rename_with(.cols=c(onset,cease,viscode,desc),.fn = ~ paste0(condition_name,"_",.x))  %>% dplyr::rename((!!condition_name_enquo):=condition)
    return(merged_data)
  }

adni_merge <- ADNIMERGE::adnimerge
adni_merge_demog <- ADNIMERGE::ptdemog %>% dplyr::filter(VISCODE != "f")
adni_lab_data <- ADNIMERGE::labdata %>% dplyr::filter(VISCODE != "f")
bateman_plasma <- ADNIMERGE::batemanlab
washu_plasma <- ADNIMERGE::plasma_abeta_project_wash_u

adni_mod_hach <- ADNIMERGE::modhach %>% dplyr::filter(VISCODE != "f")
adni_vitals <- ADNIMERGE::vitals %>% dplyr::filter(VISCODE != "f")

adni_web_demog <-
  readr::read_delim("~/Downloads/PTDEMOG.csv") %>% dplyr::filter(VISCODE !=
                                                                   "f")

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
  adni_condition_merge(upenn_merged_csf_biomarkers, adni_diabetes)
csf_all_conditions <-
  adni_condition_merge(csf_all_conditions, adni_ckd)
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

plasma_all_conditions <-
  adni_condition_merge(plasma_ratio_merged, adni_diabetes, "diabetes")
plasma_all_conditions <-
  adni_condition_merge(abeta_all_conditions, adni_ckd, "ckd")
plasma_all_conditions <-
  adni_condition_merge(plasma_all_conditions, adni_dyslipidemia, "dyslipidemia")
plasma_all_conditions <-
  adni_condition_merge(plasma_all_conditions, adni_hypertension, "hypertension")
plasma_all_conditions <-
  dplyr::left_join(plasma_all_conditions, adni_lab_data_demog_reduced, by =
                     "RID")
plasma_all_conditions <-
  dplyr::left_join(plasma_all_conditions,
                   adni_mod_hach %>% dplyr::select(RID, htn_hach),
                   by = "RID")
plasma_all_conditions <-
  dplyr::left_join(
    plasma_all_conditions,
    adni_vitals %>% dplyr::select(RID, diastolic_bp, systolic_bp, htn_vitals, vitals_date),
    by = "RID"
  )
plasma_all_conditions$date_diff <-
  ifelse(
    plasma_all_conditions$EXAMDATE - plasma_all_conditions$vitals_date < 0,
    999999,
    plasma_all_conditions$EXAMDATE - plasma_all_conditions$vitals_date
  )
plasma_all_conditions <- plasma_all_conditions %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(date_diff == min(date_diff))
plasma_all_conditions$diabetes_joined <-
  ifelse(
    plasma_all_conditions$diabetes == 1 |
      plasma_all_conditions$diabetes_cat == "Diabetes",
    1,
    0
  )
plasma_all_conditions$ckd_joined <-
  ifelse(plasma_all_conditions$ckd == 1 |
           plasma_all_conditions$ckd_cat == "Suspected CKD",
         1,
         0)
plasma_all_conditions$dyslipidemia_joined <-
  ifelse(
    plasma_all_conditions$dyslipidemia == 1 |
      plasma_all_conditions$dyslipidemia_cat == "Hyperlipidemia",
    1,
    0
  )
plasma_all_conditions$hypertension_joined <-
  ifelse(
    plasma_all_conditions$hypertension == 1 |
      (
        plasma_all_conditions$htn_vitals == 1 &
          plasma_all_conditions$date_diff != 9999
      ) | plasma_all_conditions$htn_hach == 1,
    1,
    0
  )
plasma_all_conditions <-
  dplyr::left_join(plasma_all_conditions, adni_merge_demog_uniques_reduced, by =
                     "RID")
plasma_all_conditions$age_at_exam <-
  round(as.numeric(lubridate::year(plasma_all_conditions$EXAMDATE)) -
          plasma_all_conditions$PTDOBYY +
          ((
            as.numeric(lubridate::month(plasma_all_conditions$EXAMDATE)) - plasma_all_conditions$PTDOBMM
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
