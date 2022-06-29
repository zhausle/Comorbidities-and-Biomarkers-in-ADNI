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

## PET data load-in

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
d$Centiloid = 188.22 * d$SUMMARYSUVR_WHOLECEREBNORM.x - 189.16 ## reflects most recent values for AV45 from "Instructions for converting ADNI processing results to Centiloids (PDF)" on LONI
d$Centiloid[(is.na(d$Centiloid))] = 157.15 * d$SUMMARYSUVR_WHOLECEREBNORM.y[is.na(d$Centiloid)] -
  151.87 ## reflects most recent values for FBB from "Instructions for converting ADNI processing results to Centiloids (PDF)" on LONI
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

## CSF data load-in
upenn_mk_12 <-readr::read_delim("~/Downloads/UPENNBIOMK12_01_04_21.csv")
upenn_mk_10 <-readr::read_delim("~/Downloads/UPENNBIOMK10_07_29_19.csv") %>% dplyr::rename(EXAMDATE=DRAWDATE,ABETA=ABETA42)

upenn_mk_9 <-readr::read_delim("~/Downloads/UPENNBIOMK9_04_19_17.csv")
upenn_mk_9 <- upenn_mk_9 %>% dplyr::mutate(ABETA=as.numeric(ABETA),TAU=as.numeric(TAU),PTAU=as.numeric(PTAU))
upenn_mk_9$ABETA <- ifelse(is.na(upenn_mk_9$ABETA),readr::parse_number(upenn_mk_9$COMMENT),upenn_mk_9$ABETA)

upenn_merged_csf_biomarkers <- dplyr::bind_rows(upenn_mk_12,upenn_mk_10,upenn_mk_9)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::select(-VISCODE) %>% dplyr::rename(VISCODE=VISCODE2)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::mutate(ABETA=round(ABETA))
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::group_by(RID,EXAMDATE) %>% dplyr::filter(RUNDATE==max(RUNDATE)) %>% dplyr::ungroup()

