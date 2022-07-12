library(lmerTest)

plasma_abeta_pet_use<-biomarker_response_join(plasma_abeta_all_conditions,amyloid_pet,abeta_ratio,AmyloidPos)
plasma_ptau_pet_use<-biomarker_response_join(plasma_ptau_all_conditions,amyloid_pet,ptau_181,AmyloidPos)
plasma_ptau_csf_use<-biomarker_response_join(plasma_ptau_all_conditions,csf_all_conditions,ptau_181,ptau_pos)

plasma_abeta_age_bin<-mielke_analysis_age_bin(plasma_abeta_pet_use,abeta_ratio,AmyloidPos)
plasma_abeta_cutpoints<-mielke_analysis_cutpoints(plasma_abeta_pet_use,abeta_ratio,AmyloidPos)
plasma_abeta_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_abeta_pet_use,"abeta_ratio","AmyloidPos")
plasma_abeta_univariate_regs<-mielke_analysis_univariate_analyses(plasma_abeta_pet_use,"abeta_ratio_z_scaled","AmyloidPos")

plasma_ptau_age_bin<-mielke_analysis_age_bin(plasma_ptau_pet_use,ptau_181,AmyloidPos)
plasma_ptau_cutpoints<-mielke_analysis_cutpoints(plasma_ptau_pet_use,ptau_181,AmyloidPos)
plasma_ptau_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_pet_use,"ptau_181","AmyloidPos")
plasma_ptau_univariate_regs<-mielke_analysis_univariate_analyses(plasma_ptau_pet_use,"biomarker_z_scaled","AmyloidPos")

plasma_ptau_pet_use %>% dplyr::group_by(DX=="CN",AmyloidPos) %>% dplyr::summarize(mean=mean(ptau_181,na.rm=TRUE),sd=sd(ptau_181,na.rm=TRUE),count=n())

plasma_ptau_csf_age_bin<-mielke_analysis_age_bin(plasma_ptau_csf_use,ptau_181,ptau_pos)
plasma_ptau_csf_cutpoints<-mielke_analysis_cutpoints(plasma_ptau_csf_use,ptau_181,ptau_pos)
plasma_ptau_csf_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_csf_use,"ptau_181","ptau_pos")
plasma_ptau_csf_univariate_regs<-mielke_analysis_univariate_analyses(plasma_ptau_csf_use,"ptau_181","ptau_pos")

plasma_ptau_pet_use <- plasma_ptau_pet_use %>% dplyr::mutate(dyslipidemia_cat=relevel(as.factor(dyslipidemia_cat),ref="Normal Cholesterol"),diabetes_cat=relevel(as.factor(diabetes_cat),ref="Normoglycemic"),ckd_cat=relevel(as.factor(ckd_cat),ref="No Suspected CKD"))
plasma_ptau_pet_multivariate<-lm(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos,data=plasma_ptau_pet_use)
plasma_ptau_pet_multivariate_glmm<-lmerTest::lmer(biomarker_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)
plasma_ptau_pet_multivariate_glmm_change_from_bl<-lmerTest::lmer(ptau_change_from_bl~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+ckd_cat+apoe_status+AmyloidPos+(1|RID),data=plasma_ptau_pet_use)

eGFR_test<-lm(eGFR~age_exam_scaled+PTGENDER,data=plasma_ptau_pet_use)
ckd_cat_test<-glm(as.factor(ckd_cat)~age_exam_scaled+PTGENDER,family="binomial",data=plasma_ptau_pet_use)
ckd_cat_roc<-ROCit::rocit(score=qlogis(ckd_cat_test$fitted.values),class=ckd_cat_test$y)
