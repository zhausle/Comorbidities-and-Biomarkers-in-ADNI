library(ADNIMERGE)
library(tidyverse)
library(transplantr)
## would recommend downloading fresh copies of any listed CSVs

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

## Plasma data load-in
bateman_plasma <-
  readr::read_delim("~/Projects/Comorbidities Project/batemanlab_20190621.csv")
washu_plasma <-
  readr::read_delim("~/Projects/Comorbidities Project/PLASMA_ABETA_PROJECT_WASH_U_11_05_21.csv")

bateman_plasma <-
  bateman_plasma %>% dplyr::rename(abeta_ratio = RATIO_ABETA42_40_BY_ISTD_TOUSE) %>% dplyr::mutate(origin="Bateman")
washu_plasma <-
  washu_plasma %>% dplyr::rename(abeta_ratio = STANDARDIZED_PLASMAAB4240,abeta_ratio_non_std = PLASMAAB4240) %>% dplyr::mutate(origin="WashU")

plasma_merged <- dplyr::bind_rows(
  washu_plasma %>% dplyr::filter(QC_STATUS == "Passed") %>% dplyr::mutate(MS_RUN_DATE=as.Date(MS_RUN_DATE,origin="1970-01-01")),
  bateman_plasma %>% dplyr::filter(QC_STATUS == "Passed", INSTRUMENT == "Lumos", INJECTION == "a")) %>% dplyr::filter(RID != 999999)

## Plasma ptau-181 load-in

plasma_ptau <- 
  readr::read_delim("~/Projects/Comorbidities Project/UGOTPTAU181_06_18_20.csv") %>% dplyr::select(RID,VISCODE2,EXAMDATE,PLASMAPTAU181) %>% dplyr::rename(VISCODE=VISCODE2,ptau_181=PLASMAPTAU181)

## Plasma NfL load-in
plasma_nfl <-
  readr::read_delim("~/Projects/Comorbidities Project/ADNI_BLENNOWPLASMANFLLONG_10_03_18.csv") %>% dplyr::select(RID,VISCODE2,EXAMDATE,PLASMA_NFL) %>% dplyr::rename(VISCODE=VISCODE2,nfl=PLASMA_NFL)

## PET data load-in
av45 <- readr::read_delim("~/Projects/Comorbidities Project/UCBERKELEYAV45_04_26_22.csv") %>% 
  dplyr::select(RID,EXAMDATE,VISCODE2,SUMMARYSUVR_WHOLECEREBNORM,SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF) %>% 
  dplyr::rename(suvr_summary=SUMMARYSUVR_WHOLECEREBNORM,AmyloidPos=SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF) %>%
  dplyr::mutate(Centiloid = (188.22 * suvr_summary) - 189.16,tracer="av45") ## reflects most recent values for AV45 from "Instructions for converting ADNI processing results to Centiloids (PDF)" on LONI
fbb <- readr::read_delim("~/Projects/Comorbidities Project/UCBERKELEYFBB_04_26_22.csv") %>% 
  dplyr::select(RID,EXAMDATE,VISCODE2,SUMMARYSUVR_WHOLECEREBNORM,SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF) %>% 
  dplyr::rename(suvr_summary=SUMMARYSUVR_WHOLECEREBNORM,AmyloidPos=SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF) %>%
  dplyr::mutate(Centiloid = (157.15 * suvr_summary) - 151.87,tracer="fbb") ## reflects most recent values for FBB from "Instructions for converting ADNI processing results to Centiloids (PDF)" on LONI
amyloid_pet<-dplyr::bind_rows(av45,fbb) %>% dplyr::mutate(tracer=as.factor(tracer)) %>% dplyr::rename(VISCODE=VISCODE2)

## CSF data load-in
upenn_mk_12 <-readr::read_delim("~/Projects/Comorbidities Project/UPENNBIOMK12_01_04_21.csv")
upenn_mk_10 <-readr::read_delim("~/Downloads/UPENNBIOMK10_07_29_19.csv") %>% dplyr::rename(EXAMDATE=DRAWDATE,ABETA=ABETA42)

upenn_mk_9 <-readr::read_delim("~/Projects/Comorbidities Project/UPENNBIOMK9_04_19_17.csv")
upenn_mk_9 <- upenn_mk_9 %>% dplyr::mutate(ABETA=as.numeric(ABETA),TAU=as.numeric(TAU),PTAU=as.numeric(PTAU))
upenn_mk_9$ABETA <- ifelse(is.na(upenn_mk_9$ABETA),readr::parse_number(upenn_mk_9$COMMENT),upenn_mk_9$ABETA)

upenn_merged_csf_biomarkers <- dplyr::bind_rows(upenn_mk_12,upenn_mk_10,upenn_mk_9)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::select(-VISCODE) %>% dplyr::rename(VISCODE=VISCODE2)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::mutate(ABETA=round(ABETA),ptau_pos=PTAU>24)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::group_by(RID,EXAMDATE) %>% dplyr::filter(RUNDATE==max(RUNDATE)) %>% dplyr::ungroup()

adni_mod_hach <- ADNIMERGE::modhach %>% dplyr::filter(VISCODE != "f")
adni_vitals <- ADNIMERGE::vitals %>% dplyr::filter(VISCODE != "f")

adni_web_demog <-
  readr::read_delim("~/Projects/Comorbidities Project/PTDEMOG.csv") %>% dplyr::filter(VISCODE !=
                                                                   "f")

## Demographic data
## I use ADNIMERGE here because the demographic categories are well-labeled in ADNIMERGE but not LONI
adni_lab_data <- ADNIMERGE::labdata %>% dplyr::filter(VISCODE != "f")
adni_merge_demog <- ADNIMERGE::ptdemog %>% dplyr::filter(VISCODE != "f")
adni_web_demog <-
  readr::read_delim("~/Projects/Comorbidities Project/PTDEMOG.csv") %>% dplyr::filter(VISCODE !=
                                                                   "f")

adni_merge_demog_uniques <- adni_merge_demog %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adni_web_demog_uniques <- adni_web_demog %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adni_joined_demog_uniques <-
  dplyr::left_join(
    adni_merge_demog_uniques,
    adni_web_demog_uniques %>%
      dplyr::select(PTDOBYY, PTDOBMM, RID),
    by = "RID"
  )

adni_joined_demog_uniques_reduced <-
  adni_joined_demog_uniques %>% dplyr::select(RID,PTGENDER,PTETHCAT,PTRACCAT,PTDOBYY,PTDOBMM,PTEDUCAT,ORIGPROT)

## Lab data load-in
## This is somewhat involved because I needed to calculate eGFR for CKD, which depends on age, gender, and race

adni_lab_data_demog_reduced <-
  dplyr::left_join(adni_lab_data,adni_joined_demog_uniques_reduced,by = "RID",suffix = c(".labs", ".demos"))

adni_lab_data_demog_reduced <-
  adni_lab_data_demog_reduced %>% mutate(
    glucose = as.numeric(RCT11),
    creatinine = as.numeric(RCT392),
    cholesterol = as.numeric(RCT20)
  )

adni_lab_data_demog_reduced$age_at_lab <-round(as.numeric(lubridate::year(adni_lab_data_demog_reduced$EXAMDATE)) -
    adni_lab_data_demog_reduced$PTDOBYY +
    ((as.numeric(lubridate::month(adni_lab_data_demog_reduced$EXAMDATE)) - 
        adni_lab_data_demog_reduced$PTDOBMM) / 12), digits=1)

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
  dplyr::select(RID,age_at_lab,glucose,diabetes_cat,creatinine,eGFR,ckd_cat,cholesterol,dyslipidemia_cat,race_for_egfr,gender_for_egfr,eGFR_cat
  )

adni_lab_data_use<- adni_lab_data_demog_reduced %>%
  dplyr::distinct_at(vars(RID), .keep_all = TRUE)

## Hachinski & Vitals load-in

adni_mod_hach <- ADNIMERGE::modhach %>% dplyr::filter(VISCODE != "f")
adni_mod_hach$htn_hach <-
  ifelse(adni_mod_hach$HMHYPERT == "Present - 1 point", 1, 0)

adni_vitals <- ADNIMERGE::vitals %>% dplyr::filter(VISCODE != "f")

adni_vitals <-
  adni_vitals %>% dplyr::mutate(
    diastolic_bp = as.numeric(VSBPDIA),
    systolic_bp = as.numeric(VSBPSYS),
    vitals_date = USERDATE) 

adni_vitals$htn_vitals <-
  ifelse(adni_vitals$diastolic_bp >= 90 |
           adni_vitals$systolic_bp >= 140,
         1,
         0)

## Diagnostic info load-in

adni_diagnoses<-ADNIMERGE::adnimerge %>% 
  dplyr::mutate(RID=as.numeric(RID)) %>%
  dplyr::select(RID,DX,DX.bl,EXAMDATE) %>% dplyr::rename(dx_date=EXAMDATE) %>%
  dplyr::mutate(DX.baseline = case_when(
    (DX.bl %in% c("EMCI","LMCI")) ~ "MCI",
    (DX.bl %in% c("SMC","CN")) ~ "CN",
    (DX.bl == "AD" ~ "AD")
  )) %>%
  dplyr::filter(!is.na(DX))

## APOE info load-in
adni_apoe<-ADNIMERGE::adnimerge %>% 
  dplyr::mutate(RID=as.numeric(RID)) %>%
  dplyr::select(RID,APOE4) %>% dplyr::mutate(apoe_status=(APOE4>0)) %>%
  dplyr::filter(!is.na(APOE4)) %>%
  dplyr::distinct_at(vars(RID),.keep_all = TRUE)
