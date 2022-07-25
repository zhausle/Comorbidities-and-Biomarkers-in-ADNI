library(lmerTest)
library(optimx)
library(flextable)

plasma_abeta_pet_use<-biomarker_response_join(plasma_abeta_all_conditions,amyloid_pet,abeta_ratio,AmyloidPos)
plasma_abeta_csf_use<-biomarker_response_join(plasma_abeta_all_conditions,csf_all_conditions,abeta_ratio,ptau_pos)
plasma_ptau_pet_use<-biomarker_response_join(plasma_ptau_all_conditions,amyloid_pet,ptau_181,AmyloidPos)
plasma_ptau_csf_use<-biomarker_response_join(plasma_ptau_all_conditions,csf_all_conditions,ptau_181,ptau_pos)
plasma_nfl_pet_use<-biomarker_response_join(plasma_nfl_all_conditions,pet_all_conditions,nfl,AmyloidPos)
plasma_nfl_csf_use<-biomarker_response_join(plasma_nfl_all_conditions,csf_all_conditions,nfl,ptau_pos)

plasma_ptau_pet_use<- plasma_ptau_pet_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))
plasma_abeta_pet_use<- plasma_abeta_pet_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))
plasma_nfl_pet_use<- plasma_nfl_pet_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))
plasma_ptau_csf_use<- plasma_ptau_csf_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))
plasma_nfl_csf_use<- plasma_nfl_csf_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))
plasma_abeta_csf_use<- plasma_abeta_csf_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))

plasma_abeta_pet_first_obs<- plasma_abeta_pet_use %>% dplyr::filter(!is.na(AmyloidPos)) %>% dplyr::arrange(EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::filter(EXAMDATE==EXAMDATE[1L])
plasma_ptau_pet_first_obs<- plasma_ptau_pet_use %>% dplyr::filter(!is.na(AmyloidPos)) %>% dplyr::arrange(EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::filter(EXAMDATE==EXAMDATE[1L])
plasma_nfl_pet_first_obs<- plasma_nfl_pet_use %>% dplyr::filter(!is.na(AmyloidPos)) %>% dplyr::arrange(EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::filter(EXAMDATE==EXAMDATE[1L])
plasma_abeta_csf_first_obs<- plasma_abeta_csf_use %>% dplyr::filter(!is.na(ptau_pos)) %>% dplyr::arrange(EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::filter(EXAMDATE==EXAMDATE[1L])
plasma_ptau_csf_first_obs<- plasma_ptau_csf_use %>% dplyr::filter(!is.na(ptau_pos)) %>% dplyr::arrange(EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::filter(EXAMDATE==EXAMDATE[1L])
plasma_nfl_csf_first_obs<- plasma_nfl_csf_use %>% dplyr::filter(!is.na(ptau_pos)) %>% dplyr::arrange(EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::filter(EXAMDATE==EXAMDATE[1L])


plasma_abeta_pet_roc_and_measure<-mielke_roc_mega_wrapper(plasma_abeta_pet_use,"biomarker_z_scaled","AmyloidPos","Amyloid PET ~ Plasma AB42/40")
plasma_ptau_pet_roc_and_measure<-mielke_roc_mega_wrapper(plasma_ptau_pet_use,"biomarker_z_scaled","AmyloidPos","Amyloid PET ~ Plasma ptau-181",longitudinal=TRUE)
plasma_nfl_pet_roc_and_measure<-mielke_roc_mega_wrapper(plasma_nfl_pet_use,"biomarker_z_scaled","AmyloidPos","Amyloid PET ~ Plasma NfL",longitudinal=TRUE)
plasma_abeta_csf_roc_and_measure<-mielke_roc_mega_wrapper(plasma_abeta_csf_use,"biomarker_z_scaled","ptau_pos","CSF Tau ~ Plasma AB42/40")
plasma_ptau_csf_roc_and_measure<-mielke_roc_mega_wrapper(plasma_ptau_csf_use,"biomarker_z_scaled","ptau_pos","CSF Tau ~ Plasma ptau-181",longitudinal=TRUE)
plasma_nfl_csf_roc_and_measure<-mielke_roc_mega_wrapper(plasma_nfl_csf_use,"biomarker_z_scaled","ptau_pos","CSF Tau ~ Plasma NfL",longitudinal=TRUE)


plasma_abeta_pet_roc_and_measure_first_obs<-mielke_roc_mega_wrapper(plasma_abeta_pet_first_obs,"biomarker_z_scaled","AmyloidPos","Amyloid PET ~ Plasma AB42/40 (Earliest Obs. Only)")
plasma_ptau_pet_roc_and_measure_first_obs<-mielke_roc_mega_wrapper(plasma_ptau_pet_first_obs,"biomarker_z_scaled","AmyloidPos","Amyloid PET ~ Plasma ptau-181 (Earliest Obs. Only)")
plasma_nfl_pet_roc_and_measure_first_obs<-mielke_roc_mega_wrapper(plasma_nfl_pet_first_obs,"biomarker_z_scaled","AmyloidPos","Amyloid PET ~ Plasma NfL (Earliest Obs. Only)")
plasma_abeta_csf_roc_and_measure_first_obs<-mielke_roc_mega_wrapper(plasma_abeta_csf_first_obs,"biomarker_z_scaled","ptau_pos","CSF Tau ~ Plasma AB42/40 (Earliest Obs. Only)")
plasma_ptau_csf_roc_and_measure_first_obs<-mielke_roc_mega_wrapper(plasma_ptau_csf_first_obs,"biomarker_z_scaled","ptau_pos","CSF Tau ~ Plasma ptau-181 (Earliest Obs. Only)")
plasma_nfl_csf_roc_and_measure_first_obs<-mielke_roc_mega_wrapper(plasma_nfl_csf_first_obs,"biomarker_z_scaled","ptau_pos","CSF Tau ~ Plasma NfL (Earliest Obs. Only)")


plasma_abeta_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_abeta_pet_use,"biomarker_z_scaled","AmyloidPos")
mielke_joined_roc_plots(plasma_abeta_roc_and_measure,"Plasma AB42/40 Logistic Models for Amyloid PET")
mielke_joined_logistic_measures(plasma_abeta_roc_and_measure,"Plasma AB42/40 Logistic Models for Amyloid PET")
mielke_coef_table_fn(plasma_abeta_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]],"Amyloid PET ~ Plasma AB42/40 Full Model")
mielke_coef_table_fn(plasma_abeta_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]],"Amyloid PET ~ Plasma AB42/40 Full Model, Lab Comorbidities Only")  

plasma_abeta_roc_and_measure_comparisons<-mielke_analysis_roc_and_measure_comparisons(plasma_abeta_pet_use,"biomarker_z_scaled","AmyloidPos")
mielke_joined_roc_plots_comparisons(plasma_abeta_roc_and_measure_comparisons,"Plasma AB42/40 Logistic Models for Amyloid PET")
mielke_joined_logistic_measures_comparisons(plasma_abeta_roc_and_measure_comparisons,"Plasma AB42/40 Logistic Models for Amyloid PET")

test_cutoff<-plasma_abeta_roc_and_measure_comparisons$`All Observations`[[4]][which.max(plasma_abeta_roc_and_measure_comparisons$`All Observations`[[4]]$Youden),] ## this should be correct but unsure


test<-ROCit::rocit(score=-plasma_abeta_pet_first_obs$biomarker_z_scaled,class=plasma_abeta_pet_first_obs$AmyloidPos)
test<-ROCit::rocit(score=plasma_ptau_pet_first_obs$biomarker_z_scaled,class=plasma_ptau_pet_first_obs$AmyloidPos)
test<-ROCit::rocit(score=plasma_ptau_pet_first_obs$biomarker_z_scaled,class=plasma_ptau_pet_first_obs$AmyloidPos)



flextable::flextable(as.data.frame(summary(plasma_abeta_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]])$coefficients) %>% tibble::rownames_to_column(var="Predictor")) %>% flextable::autofit() %>% flextable::colformat_double(digits=3) %>% flextable::add_header_lines(value=title)
summary(plasma_abeta_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]])
mielke_joined_logistic_measures(plasma_abeta_roc_and_measure,"Plasma AB42/40 Logistic Models for Amyloid PET")


plasma_abeta_age_bin<-mielke_analysis_age_bin(plasma_abeta_pet_use,abeta_ratio,AmyloidPos)
plasma_abeta_cutpoints<-mielke_analysis_cutpoints(plasma_abeta_pet_use,abeta_ratio,AmyloidPos)
plasma_abeta_univariate_regs<-mielke_analysis_univariate_analyses(plasma_abeta_pet_use,"biomarker_z_scaled","AmyloidPos")



plasma_abeta_csf_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_abeta_csf_use,"biomarker_z_scaled","ptau_pos")
mielke_joined_roc_plots(plasma_abeta_csf_roc_and_measure)

plasma_abeta_csf_roc_and_measure$Unadjusted[[2]]$AUC
plasma_abeta_csf_roc_and_measure$`Adjusted for Age & Sex`[[2]]$AUC
plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, and APOE`[[2]]$AUC
plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$AUC
plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$AUC
plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$AUC
plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, and Comorbidities`[[2]]$AUC

plasma_abeta_csf_roc_and_measure_comparisons<-mielke_analysis_roc_and_measure_comparisons(plasma_abeta_csf_first_obs,"biomarker_z_scaled","ptau_pos")
mielke_joined_roc_plots_comparisons(plasma_abeta_roc_and_measure_comparisons,"Plasma AB42/40 Logistic Models for CSF Tau")
mielke_joined_logistic_measures_comparisons(plasma_abeta_csf_roc_and_measure_comparisons,"Plasma AB42/40 Logistic Models for CSF Tau")


summary(plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]])
summary(plasma_abeta_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]])

plasma_ptau_age_bin<-mielke_analysis_age_bin(plasma_ptau_pet_use,ptau_181,AmyloidPos)
plasma_ptau_cutpoints<-mielke_analysis_cutpoints(plasma_ptau_pet_use,ptau_181,AmyloidPos)
plasma_ptau_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_pet_use,"biomarker_z_scaled","AmyloidPos",longitudinal=TRUE)
mielke_joined_roc_plots(plasma_ptau_roc_and_measure)
title(main="Amyloid PET Status ~ Plasma ptau-181 + Covariates")
plasma_ptau_univariate_regs<-mielke_analysis_univariate_analyses(plasma_ptau_pet_use,"biomarker_z_scaled","AmyloidPos")

plasma_ptau_roc_and_measure$Unadjusted[[2]]$AUC
plasma_ptau_roc_and_measure$`Adjusted for Age & Sex`[[2]]$AUC
plasma_ptau_roc_and_measure$`Adjusted for Age, Sex, and APOE`[[2]]$AUC
plasma_ptau_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$AUC
plasma_ptau_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$AUC
plasma_ptau_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$AUC
plasma_ptau_roc_and_measure$`Adjusted for Age, Sex, and Comorbidities`[[2]]$AUC

summary(plasma_ptau_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]])
coef(summary(plasma_ptau_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]]))[,"Estimate"]
summary(plasma_ptau_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]])

plasma_ptau_csf_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_csf_use,"biomarker_z_scaled","ptau_pos",longitudinal=TRUE)
mielke_joined_roc_plots(plasma_ptau_csf_roc_and_measure)
title(main="Tau CSF Status ~ Plasma ptau-181 + Covariates")

plasma_ptau_pet_use %>% dplyr::group_by(DX=="CN",AmyloidPos) %>% dplyr::summarize(mean=mean(ptau_181,na.rm=TRUE),sd=sd(ptau_181,na.rm=TRUE),count=n())

plasma_nfl_age_bin<-mielke_analysis_age_bin(plasma_nfl_pet_use,nfl,AmyloidPos)
plasma_nfl_cutpoints<-mielke_analysis_cutpoints(plasma_nfl_pet_use,ptau_181,AmyloidPos)
plasma_nfl_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_nfl_pet_use,"biomarker_z_scaled","AmyloidPos",longitudinal=TRUE)
plasma_nfl_univariate_regs<-mielke_analysis_univariate_analyses(plasma_nfl_pet_use,"biomarker_z_scaled","AmyloidPos")

mielke_joined_roc_plots(plasma_nfl_roc_and_measure)
title(main="Amyloid PET Status ~ Plasma NfL + Covariates")

summary(plasma_nfl_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]])
summary(plasma_nfl_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]])

plasma_ptau_csf_age_bin<-mielke_analysis_age_bin(plasma_ptau_csf_use,ptau_181,ptau_pos)
plasma_ptau_csf_cutpoints<-mielke_analysis_cutpoints(plasma_ptau_csf_use,ptau_181,ptau_pos)
plasma_ptau_csf_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_csf_use,"biomarker_z_scaled","ptau_pos",longitudinal=TRUE)
plasma_ptau_csf_univariate_regs<-mielke_analysis_univariate_analyses(plasma_ptau_csf_use,"ptau_181","ptau_pos")

summary(plasma_ptau_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]])
summary(plasma_ptau_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]])

plasma_nfl_csf_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_nfl_csf_use,"biomarker_z_scaled","ptau_pos",longitudinal=TRUE)

plasma_nfl_csf_roc_and_measure$Unadjusted[[2]]$AUC
plasma_nfl_csf_roc_and_measure$`Adjusted for Age & Sex`[[2]]$AUC
plasma_nfl_csf_roc_and_measure$`Adjusted for Age, Sex, and APOE`[[2]]$AUC
plasma_nfl_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$AUC
plasma_nfl_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$AUC
plasma_nfl_csf_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$AUC
plasma_nfl_csf_roc_and_measure$`Adjusted for Age, Sex, and Comorbidities`[[2]]$AUC

plasma_ptau_pet_use <- plasma_ptau_pet_use %>% dplyr::mutate(dyslipidemia_cat=relevel(as.factor(dyslipidemia_cat),ref="Normal Cholesterol"),diabetes_cat=relevel(as.factor(diabetes_cat),ref="Normoglycemic"),ckd_cat=relevel(as.factor(ckd_cat),ref="No Suspected CKD"))
plasma_ptau_pet_multivariate<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_ptau_pet_use)
plasma_ptau_pet_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)
plasma_ptau_pet_multivariate_glmm_change_from_bl<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)

plasma_abeta_pet_use <- plasma_abeta_pet_use %>% dplyr::mutate(dyslipidemia_cat=relevel(as.factor(dyslipidemia_cat),ref="Normal Cholesterol"),diabetes_cat=relevel(as.factor(diabetes_cat),ref="Normoglycemic"),ckd_cat=relevel(as.factor(ckd_cat),ref="No Suspected CKD"),eGFR_cat=relevel(as.factor(eGFR_cat),ref="Less than 45"))

plasma_abeta_pet_multivariate<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_abeta_pet_use)
sjPlot::plot_model(plasma_abeta_pet_multivariate) + ggtitle("LM Coefficients for Plasma AB42/40 (PET Data)")
mielke_coef_table_fn(plasma_abeta_pet_multivariate,"LM Coefficients for Plasma AB42/40 (PET Data)")

plasma_abeta_csf_multivariate<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos,data=plasma_abeta_csf_use)
sjPlot::plot_model(plasma_abeta_csf_multivariate) + ggtitle("LM Coefficients for Plasma AB42/40 (CSF Data)")
mielke_coef_table_fn(plasma_abeta_csf_multivariate,"LM Coefficients for Plasma AB42/40 (CSF Data)")

plasma_ptau_pet_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)
sjPlot::plot_model(plasma_ptau_pet_multivariate_glmm) + ggtitle("LMM Coefficients for ptau-181 (PET data)")
mielke_coef_table_fn(plasma_ptau_pet_multivariate_glmm,"LMM Coefficients for ptau-181 (PET data)")

plasma_ptau_csf_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_ptau_csf_use)
sjPlot::plot_model(plasma_ptau_csf_multivariate_glmm) + ggtitle("LMM Coefficients for ptau-181 (CSF data)")
mielke_coef_table_fn(plasma_ptau_csf_multivariate_glmm,"LMM Coefficients for ptau-181 (CSF data)")

plasma_ptau_pet_multivariate_glmm_change_from_bl<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)
sjPlot::plot_model(plasma_ptau_pet_multivariate_glmm_change_from_bl) + ggtitle("LMM Coefficients for ptau-181 Change from Baseline (PET data)")
mielke_coef_table_fn(plasma_ptau_pet_multivariate_glmm_change_from_bl,"LMM Coefficients for ptau-181 Change from Baseline (PET data)")

plasma_ptau_csf_multivariate_glmm_change_from_bl<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_ptau_csf_use)
sjPlot::plot_model(plasma_ptau_csf_multivariate_glmm_change_from_bl) + ggtitle("LMM Coefficients for ptau-181 Change from Baseline (CSF data)")
mielke_coef_table_fn(plasma_ptau_csf_multivariate_glmm_change_from_bl,"LMM Coefficients for ptau-181 Change from Baseline (CSF data)")

plasma_nfl_pet_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_nfl_pet_use)
sjPlot::plot_model(plasma_nfl_pet_multivariate_glmm) + ggtitle("LMM Coefficients for NfL (PET data)")
mielke_coef_table_fn(plasma_nfl_pet_multivariate_glmm,"LMM Coefficients for NfL (PET data)")

plasma_nfl_csf_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_nfl_csf_use)
sjPlot::plot_model(plasma_nfl_csf_multivariate_glmm) + ggtitle("LMM Coefficients for NfL (CSF data)")
mielke_coef_table_fn(plasma_nfl_csf_multivariate_glmm,"LMM Coefficients for NfL (CSF data)")

plasma_nfl_pet_multivariate_glmm_change_from_bl<-lmerTest::lmer(nfl_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_nfl_pet_use)
sjPlot::plot_model(plasma_nfl_pet_multivariate_glmm_change_from_bl) + ggtitle("LMM Coefficients for NfL Change from Baseline (PET data)")
mielke_coef_table_fn(plasma_nfl_pet_multivariate_glmm_change_from_bl,"LMM Coefficients for NfL Change from Baseline (PET data)")

plasma_nfl_csf_multivariate_glmm_change_from_bl<-lmerTest::lmer(nfl_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_nfl_csf_use)
sjPlot::plot_model(plasma_nfl_csf_multivariate_glmm_change_from_bl) + ggtitle("LMM Coefficients for NfL Change from Baseline (CSF data)")
mielke_coef_table_fn(plasma_nfl_csf_multivariate_glmm_change_from_bl,"LMM Coefficients for NfL Change from Baseline (CSF data)")

plasma_abeta_pet_multivariate_first_obs<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_abeta_pet_first_obs)
sjPlot::plot_model(plasma_abeta_pet_multivariate_first_obs) + ggtitle("LM Coefficients for Plasma AB42/40 (PET Data, Earliest Obs. Only)")
mielke_coef_table_fn(plasma_abeta_pet_multivariate_first_obs,"LM Coefficients for Plasma AB42/40 (PET Data, Earliest Obs. Only)")

plasma_abeta_csf_multivariate_first_obs<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos,data=plasma_abeta_csf_first_obs)
sjPlot::plot_model(plasma_abeta_csf_multivariate_first_obs) + ggtitle("LM Coefficients for Plasma AB42/40 (CSF Data, Earliest Obs. Only)")
mielke_coef_table_fn(plasma_abeta_csf_multivariate_first_obs,"LM Coefficients for Plasma AB42/40 (CSF Data, Earliest Obs. Only)")

plasma_ptau_pet_multivariate_lm_first_obs<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_ptau_pet_first_obs)
sjPlot::plot_model(plasma_ptau_pet_multivariate_lm_first_obs) + ggtitle("LM Coefficients for ptau-181 (PET data, Earliest Obs. Only)")
mielke_coef_table_fn(plasma_ptau_pet_multivariate_lm_first_obs,"LM Coefficients for ptau-181 (PET data, Earliest Obs. Only)")

plasma_ptau_csf_multivariate_lm_first_obs<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos,data=plasma_ptau_csf_first_obs)
sjPlot::plot_model(plasma_ptau_csf_multivariate_lm_first_obs) + ggtitle("LM Coefficients for ptau-181 (CSF data, Earliest Obs. Only)")
mielke_coef_table_fn(plasma_ptau_csf_multivariate_lm_first_obs,"LM Coefficients for ptau-181 (CSF data, Earliest Obs. Only)")

plasma_nfl_pet_multivariate_lm_first_obs<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_nfl_pet_first_obs)
sjPlot::plot_model(plasma_nfl_pet_multivariate_lm_first_obs) + ggtitle("LM Coefficients for nfl (PET data, Earliest Obs. Only)")
mielke_coef_table_fn(plasma_nfl_pet_multivariate_lm_first_obs,"LM Coefficients for nfl (PET data, Earliest Obs. Only)")

plasma_nfl_csf_multivariate_lm_first_obs<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos,data=plasma_nfl_csf_first_obs)
sjPlot::plot_model(plasma_nfl_csf_multivariate_lm_first_obs) + ggtitle("LM Coefficients for nfl (CSF data, Earliest Obs. Only)")
mielke_coef_table_fn(plasma_nfl_csf_multivariate_lm_first_obs,"LM Coefficients for nfl (CSF data, Earliest Obs. Only)")


plasma_nfl_csf_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_ptau_csf_use)

plasma_ptau_pet_multivariate_glmm_change_from_bl<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)
plasma_nfl_pet_multivariate_glmm_change_from_bl<-lmerTest::lmer(nfl_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_nfl_pet_use)
plasma_ptau_csf_multivariate_glmm_change_from_bl<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_ptau_csf_use)
plasma_nfl_csf_multivariate_glmm_change_from_bl<-lmerTest::lmer(nfl_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+ptau_pos+(1|RID),data=plasma_nfl_csf_use)

plasma_ptau_pet_multivariate_glmm_change_from_bl_eGFR<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+eGFR_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)

plasma_nfl_pet_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_nfl_pet_use)
sjPlot::plot_model(plasma_nfl_pet_multivariate_glmm) + ggtitle("GLM Coefficients for NfL (PET data)")
summary(plasma_nfl_pet_multivariate_glmm)
plasma_nfl_pet_multivariate_glmm_change_from_bl_eGFR<-lmerTest::lmer(nfl_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+eGFR_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_nfl_pet_use)

plasma_ptau_csf_use <- plasma_ptau_csf_use %>% dplyr::mutate(dyslipidemia_cat=relevel(as.factor(dyslipidemia_cat),ref="Normal Cholesterol"),diabetes_cat=relevel(as.factor(diabetes_cat),ref="Normoglycemic"),ckd_cat=relevel(as.factor(ckd_cat),ref="No Suspected CKD"),eGFR_cat=relevel(as.factor(eGFR_cat),ref="Less than 45"))
plasma_ptau_csf_multivariate_glmm_change_from_bl_eGFR<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+eGFR_cat+apoe_status+ptau_pos+(1|RID),data=plasma_ptau_csf_use)

plasma_nfl_csf_use <- plasma_nfl_csf_use %>% dplyr::mutate(dyslipidemia_cat=relevel(as.factor(dyslipidemia_cat),ref="Normal Cholesterol"),diabetes_cat=relevel(as.factor(diabetes_cat),ref="Normoglycemic"),ckd_cat=relevel(as.factor(ckd_cat),ref="No Suspected CKD"),eGFR_cat=relevel(as.factor(eGFR_cat),ref="Less than 45"))
plasma_nfl_csf_multivariate_glmm_change_from_bl_eGFR<-lmerTest::lmer(nfl_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+eGFR_cat+apoe_status+ptau_pos+(1|RID),data=plasma_nfl_csf_use)

eGFR_test<-lm(eGFR~age_exam_scaled+PTGENDER,data=plasma_ptau_pet_use)
ckd_cat_test<-glm(as.factor(ckd_cat)~age_exam_scaled+PTGENDER,family="binomial",data=plasma_ptau_pet_use)
ckd_cat_roc<-ROCit::rocit(score=qlogis(ckd_cat_test$fitted.values),class=ckd_cat_test$y)

table(adni_lab_data$VISCODE)
plasma_ptau_pet_earliest_obs<- plasma_ptau_pet_use %>% dplyr::arrange(RID,EXAMDATE) %>% dplyr::distinct_at(vars(RID),.keep_all = TRUE)
plasma_ptau_earliest_obs_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_pet_earliest_obs,"ptau_181","AmyloidPos")

plasma_ptau_earliest_obs_roc_and_measure$Unadjusted[[2]]$AUC
plasma_ptau_earliest_obs_roc_and_measure$`Adjusted for Age & Sex`[[2]]$AUC
plasma_ptau_earliest_obs_roc_and_measure$`Adjusted for Age, Sex, and APOE`[[2]]$AUC
plasma_ptau_earliest_obs_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$AUC
plasma_ptau_earliest_obs_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$AUC
plasma_ptau_earliest_obs_roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$AUC

plasma_ptau_pet_earliest_obs_multivariate<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_ptau_pet_earliest_obs)

plasma_ptau_pet_use<- plasma_ptau_pet_use %>% dplyr::mutate(age_exam_scaled_for_glmm=scale(age_at_exam))

test_glm<-glm(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","PTGENDER","age_exam_scaled_for_glmm",
                                                  "dyslipidemia_cat","htn_vitals","diabetes_cat","as.factor(apoe_status)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use)

test_roc_glm<-  ROCit::rocit(score=qlogis(test_glm$fitted.values),class=test_glm$y)
test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","(1|RID)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","PTGENDER","age_exam_scaled_for_glmm","apoe_status","(1|RID)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","PTGENDER","age_exam_scaled_for_glmm",
                                                  "dyslipidemia_joined","hypertension_joined","diabetes_joined","ckd_joined","apoe_status","(1|RID)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","PTGENDER","age_exam_scaled_for_glmm",
                                                  "dyslipidemia_cat","htn_vitals","diabetes_cat","apoe_status","creatinine","(1|RID)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","PTGENDER","age_exam_scaled_for_glmm",
                                                  "dyslipidemia","hypertension","diabetes","ckd","apoe_status","(1|RID)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("biomarker_z_scaled","PTGENDER","age_exam_scaled_for_glmm",
                                                  "dyslipidemia_joined","hypertension_joined","diabetes_joined","ckd_joined","apoe_status","(1|RID)",sep=" + ")),family=binomial,data=plasma_nfl_pet_use,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
test<-lme4::glmer(AmyloidPos ~ (1|RID),family=binomial,data=plasma_ptau_pet_use,control=glmerControl(optimizer="bobyqa"))

test<-lme4::glmer(paste0("AmyloidPos"," ~ ",paste("(1|RID)","(age_exam_scaled_for_glmm|RID)",sep=" + ")),family=binomial,data=plasma_ptau_pet_use,control=glmerControl(optimizer="bobyqa"))
test_roc<-  ROCit::rocit(score=predict(test,re.form=NA),class=test@resp$`.->y`)
                         
predict(test,re.form=NA)

plasma_abeta_pet_rf<-randomForest::randomForest(as.factor(plasma_abeta_pet_first_obs$AmyloidPos) ~ .,importance=TRUE,data=plasma_abeta_pet_first_obs %>% dplyr::select(AmyloidPos,biomarker_z_scaled,PTGENDER,PTRACCAT,age_at_exam,
                                                                                                                   ckd_joined,dyslipidemia_joined,hypertension_joined,
                                                                                                                                            diabetes_joined,apoe_status,stroke,heart_attack,PTETHCAT,PTEDUCAT,glucose,creatinine,eGFR,systolic_bp,diastolic_bp,cholesterol),na.action=na.omit)
randomForest::varImpPlot(plasma_abeta_pet_rf)
plasma_abeta_rf_roc<-ROCit::rocit(plasma_abeta_pet_rf$votes[,2],plasma_abeta_pet_rf$y)


plasma_ptau_pet_rf<-randomForest::randomForest(as.factor(plasma_ptau_pet_first_obs$AmyloidPos) ~ .,importance=TRUE,data=plasma_ptau_pet_first_obs %>% dplyr::select(AmyloidPos,biomarker_z_scaled,PTGENDER,PTRACCAT,age_at_exam,
                                                                                                                                                                       ckd_joined,dyslipidemia_joined,hypertension_joined,
                                                                                                                                                                       diabetes_joined,apoe_status,stroke,heart_attack,PTETHCAT,PTEDUCAT,glucose,creatinine,eGFR,systolic_bp,diastolic_bp,cholesterol),na.action=na.omit)

randomForest::varImpPlot(plasma_ptau_pet_rf)
randomForest::importance(plasma_ptau_pet_rf)
plasma_ptau_rf_roc<-ROCit::rocit(plasma_ptau_pet_rf$votes[,2],plasma_ptau_pet_rf$y)

plasma_nfl_pet_rf<-randomForest::randomForest(as.factor(plasma_nfl_pet_first_obs$AmyloidPos) ~ .,importance=TRUE,data=plasma_nfl_pet_first_obs %>% dplyr::select(AmyloidPos,biomarker_z_scaled,PTGENDER,PTRACCAT,age_at_exam,
                                                                                                                                                                    ckd_joined,dyslipidemia_joined,hypertension_joined,
                                                                                                                                                                    diabetes_joined,apoe_status,stroke,heart_attack,PTETHCAT,PTEDUCAT,glucose,creatinine,eGFR,systolic_bp,diastolic_bp,cholesterol),na.action=na.omit)

randomForest::varImpPlot(plasma_nfl_pet_rf)
randomForest::importance(plasma_nfl_pet_rf)
plasma_nfl_rf_roc<-ROCit::rocit(plasma_nfl_pet_rf$votes[,2],plasma_nfl_pet_rf$y)

