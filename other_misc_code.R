plasma_abeta_pet_use<-dplyr::left_join(plasma_all_conditions,amyloid_pet %>% dplyr::select(RID,EXAMDATE,AmyloidPos) %>% dplyr::rename(pet_date=EXAMDATE),by="RID")
plasma_abeta_pet_use<-date_comparison(pet_plasma_use,EXAMDATE,pet_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60) %>% dplyr::filter(retain_flag==1|alt_flag==1)
plasma_abeta_pet_use<-pet_plasma_use %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE) 
plasma_abeta_pet_use<- pet_plasma_use %>% dplyr::ungroup() %>% dplyr::mutate(abeta_ratio_z_scaled=(abeta_ratio-mean(abeta_ratio,na.rm=TRUE))/sd(abeta_ratio,na.rm=TRUE))

plasma_abeta_age_bin<-mielke_analysis_age_bin(plasma_abeta_pet_use,abeta_ratio,AmyloidPos)
plasma_abeta_cutpoints<-mielke_analysis_cutpoints(plasma_abeta_pet_use,abeta_ratio,AmyloidPos)
plasma_abeta_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_abeta_pet_use,"abeta_ratio","AmyloidPos")
plasma_abeta_univariate_regs<-mielke_analysis_univariate_analyses(plasma_abeta_pet_use,"abeta_ratio_z_scaled","AmyloidPos")

plasma_ptau_pet_use<-dplyr::left_join(plasma_ptau_all_conditions,amyloid_pet %>% dplyr::select(RID,EXAMDATE,AmyloidPos) %>% rename(pet_date=EXAMDATE),by="RID")
plasma_ptau_pet_use<-date_comparison(plasma_ptau_pet_use,EXAMDATE,pet_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60) %>% dplyr::filter(retain_flag==1|alt_flag==1)
plasma_ptau_pet_use<-plasma_ptau_pet_use %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE) 
plasma_ptau_pet_use<- plasma_ptau_pet_use %>% dplyr::ungroup() %>% dplyr::mutate(ptau_z_scaled=(ptau_181-mean(ptau_181,na.rm=TRUE))/sd(ptau_181,na.rm=TRUE))

biomarker_response_join<-function(biomarker_data,response_data,biomarker,response){
  
biomarker_enquo=enquo(biomarker)
response_enquo=enquo(response)

joined_data<-dplyr::left_join(biomarker_data, response_data %>% dplyr::select(RID,EXAMDATE,!!response_enquo) %>% rename(response_date=EXAMDATE),by="RID") 
joined_data<-date_comparison(joined_data,EXAMDATE,response_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60) %>% dplyr::filter(retain_flag==1|alt_flag==1)
joined_data<-joined_data %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE) 
return(joined_data)
}

test<-biomarker_response_join(plasma_ptau_all_conditions,pet_all_conditions,ptau_181,AmyloidPos)

plasma_ptau_age_bin<-mielke_analysis_age_bin(plasma_ptau_pet_use,ptau_181,AmyloidPos)
plasma_ptau_cutpoints<-mielke_analysis_cutpoints(plasma_ptau_pet_use,ptau_181,AmyloidPos)
plasma_ptau_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_pet_use,"ptau_181","AmyloidPos")
plasma_ptau_univariate_regs<-mielke_analysis_univariate_analyses(plasma_ptau_pet_use,"ptau_z_scaled","AmyloidPos")

plasma_ptau_pet_use %>% dplyr::group_by(DX=="CN",AmyloidPos) %>% dplyr::summarize(mean=mean(ptau_181,na.rm=TRUE),sd=sd(ptau_181,na.rm=TRUE))

plasma_ptau_csf_use<-dplyr::left_join(plasma_ptau_all_conditions,csf_all_conditions %>% dplyr::select(RID,EXAMDATE,ptau_pos) %>% rename(csf_date=EXAMDATE),by="RID")
plasma_ptau_csf_use<-date_comparison(plasma_ptau_csf_use,EXAMDATE,csf_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60) %>% dplyr::filter(retain_flag==1|alt_flag==1)
plasma_ptau_csf_use<- plasma_ptau_csf_use %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE) 
plasma_ptau_csf_use<- plasma_ptau_csf_use %>% dplyr::ungroup() %>% dplyr::mutate(ptau_z_scaled=(ptau_181-mean(ptau_181,na.rm=TRUE))/sd(ptau_181,na.rm=TRUE))

plasma_ptau_csf_age_bin<-mielke_analysis_age_bin(plasma_ptau_csf_use,ptau_181,ptau_pos)
plasma_ptau_csf_cutpoints<-mielke_analysis_cutpoints(plasma_ptau_csf_use,ptau_181,ptau_pos)
plasma_ptau_csf_roc_and_measure<-mielke_analysis_roc_and_measure_outer(plasma_ptau_csf_use,"ptau_181","ptau_pos")
plasma_ptau_csf_univariate_regs<-mielke_analysis_univariate_analyses(plasma_ptau_csf_use,"ptau_181","ptau_pos")

plasma_ptau_pet_multivariate<-lm(ptau_z_scaled~age_exam_scaled+htn_vitals+apoe_status+dyslipidemia_cat+PTGENDER+diabetes_cat+eGFR_scaled+apoe_status+AmyloidPos,data=plasma_ptau_pet_use)
eGFR_test<-lm(eGFR~age_exam_scaled+PTGENDER,data=plasma_ptau_pet_use)
ckd_cat_test<-glm(as.factor(ckd_cat)~age_exam_scaled+PTGENDER,family="binomial",data=plasma_ptau_pet_use)
ckd_cat_roc<-ROCit::rocit(score=qlogis(ckd_cat_test$fitted.values),class=ckd_cat_test$y)
