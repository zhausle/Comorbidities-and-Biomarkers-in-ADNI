biomarker_response_join<-function(biomarker_data,response_data,biomarker,response){
  
  biomarker_enquo=enquo(biomarker)
  response_enquo=enquo(response)
  
  joined_data<-dplyr::left_join(biomarker_data, response_data %>% dplyr::select(RID,EXAMDATE,!!response_enquo) %>% rename(response_date=EXAMDATE),by="RID") 
  joined_data<-date_comparison(joined_data,EXAMDATE,response_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60) %>% dplyr::filter(retain_flag==1|alt_flag==1)
  joined_data<-joined_data %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE)
  joined_data<-joined_data %>% dplyr::ungroup() %>% dplyr::mutate(biomarker_z_scaled=((!!biomarker_enquo-mean(!!biomarker_enquo,na.rm=TRUE)))/sd(!!biomarker_enquo,na.rm=TRUE))
  return(joined_data)
}

mielke_analysis_age_bin <- function(data,biomarker,response){
  biomarker_enquo<-enquo(biomarker)
  response_enquo<-enquo(response)
  
  biomarker_by_age_bin_all <- data %>% dplyr::group_by(age_exam_cat) %>% dplyr::filter((!is.na(age_exam_cat)) & age_at_exam>65) %>% dplyr::summarize(mean=mean(!!biomarker_enquo,na.rm=TRUE),
                                                                                                                                                     se=sd(!!biomarker_enquo,na.rm=TRUE)/sqrt(n())) %>% dplyr::mutate(upper=mean+se,lower=mean-se)
  
  biomarker_by_age_bin_am_pos <- data %>% dplyr::group_by(age_exam_cat,!!response_enquo) %>% dplyr::filter((!is.na(age_exam_cat)) & age_at_exam>65) %>% dplyr::summarize(mean=mean(!!biomarker_enquo,na.rm=TRUE),
                                                                                                                                                                         se=sd(!!biomarker_enquo,na.rm=TRUE)/sqrt(n())) %>% dplyr::mutate(upper=mean+se,lower=mean-se)
  
  print(ggplot2::ggplot(data=biomarker_by_age_bin_all,
                        aes(x=age_exam_cat,y=mean)) + 
          geom_point(position=position_dodge(width=1)) +
          geom_errorbar(aes(ymin=mean-se,ymax=mean+se),position=position_dodge(width=1)))
  
  print(ggplot2::ggplot(data=biomarker_by_age_bin_am_pos,
                        aes(x=age_exam_cat,y=mean,color=as.factor(!!response_enquo))) + 
          geom_point(position=position_dodge(width=1)) +
          geom_errorbar(aes(ymin=mean-se,ymax=mean+se),position=position_dodge(width=1)))
  
  results <- list(biomarker_by_age_bin_all,biomarker_by_age_bin_am_pos)
  
  return(results)
}


mielke_analysis_cutpoints <- function(data,biomarker,response){
  biomarker_enquo<-enquo(biomarker)
  response_enquo<-enquo(response)
  
  cutpoint_no_comorbidities <- data %>% dplyr::ungroup() %>% dplyr::filter(DX=="CN",!!response_enquo==0,ckd_joined==0) %>%
    dplyr::summarize(cutpoint=mean(!!biomarker_enquo)-(1.96*sd(!!biomarker_enquo)),n_obs=n())
  
  cutpoint_inc_comorbidities <- data %>% dplyr::ungroup() %>% dplyr::filter(DX=="CN",!!response_enquo==0) %>%
    dplyr::summarize(cutpoint=mean(!!biomarker_enquo)-(1.96*sd(!!biomarker_enquo)),n_obs=n())
  
  cutpoint_df <- data %>% dplyr::mutate(pred_response_no_comorbidities=(!!response_enquo<=cutpoint_no_comorbidities),
                                        pred_response_inc_comorbidities=(!!response_enquo<=cutpoint_inc_comorbidities)) %>%
    dplyr::select(pred_response_no_comorbidities,pred_response_inc_comorbidities,!!response_enquo)
  
  results<-list(cutpoint_no_comorbidities,cutpoint_inc_comorbidities,cutpoint_df)
  
  return(results)
}

mielke_analysis_roc_and_measure_inner <- function(formula,data){
  model<-glm(formula,family=binomial,data=data)
  roc<-ROCit::rocit(score=qlogis(model$fitted.values),class=model$y)
  measure<-ROCit::measureit(score=qlogis(model$fitted.values),class=model$y,
                            measure=c("ACC","SENS","SPEC","PPV","NPV"))
  measure_table<-cbind.data.frame(measure$Cutoff,measure$ACC,
                                  measure$SENS, measure$SPEC,
                                  measure$PPV,measure$NPV) %>% dplyr::mutate(youden=measure$SENS+measure$SPEC-1)
  colnames(measure_table)<-c("cutoff","ACC","SENS","SPEC","PPV","NPV","Youden")
  results<-list(model,roc,measure,measure_table)
  plot(roc)
  return(results)
}

mielke_analysis_roc_and_measure_outer <- function(data,biomarker,response){
  unadj_formula<-paste0(response," ~ ",biomarker)
  adj_formula<-paste0(response," ~ ",paste(biomarker,"age_at_exam","PTGENDER","apoe_status",sep=" + "))
  adj_formula_w_comorbidities<-paste0(response," ~ ",paste(biomarker,"age_at_exam","PTGENDER",
                                                           "dyslipidemia_joined","ckd_joined","hypertension_joined","diabetes_joined","apoe_status",sep=" + "))
  adj_formula_w_comorbidities_labs<-paste0(response," ~ ",paste(biomarker,"age_at_exam","PTGENDER",
                                                           "dyslipidemia_cat","ckd_cat","htn_vitals","diabetes_cat","apoe_status",sep=" + "))
  unadj<-mielke_analysis_roc_and_measure_inner(unadj_formula,data)
  adj<-mielke_analysis_roc_and_measure_inner(adj_formula,data)
  adj_comorbidities<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities,data)
  adj_comorbidities_labs<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities_labs,data)
  results<-list(unadj,adj,adj_comorbidities,adj_comorbidities_labs)
  results<-set_names(results,c("Unadjusted","Adjusted for Age & Sex","Adjusted for Age, Sex, and Comorbidities","Adjusted for Age, Sex, and Comorbidities (Labs Only)"))
  return(results)}

mielke_analysis_univariate_analyses <- function(data,biomarker,response){
  
  lm_vars<-c("diabetes_joined","hypertension_joined","ckd_joined","apoe_status","PTEDUCAT","PTGENDER","age_exam_scaled","dyslipidemia_joined","eGFR_scaled",response)
  lm_var_names<-c("Diabetes","HTN","CKD","Has APOE4 Allele","Education (Yrs)","Gender","Age at Plasma Exam (10 Yrs.)","Dyslipidemia","eGFR (10 u.)","Gold Standard+")
  
  lm_vars_adj<-c("diabetes_joined","hypertension_joined","ckd_joined","apoe_status","PTEDUCAT","dyslipidemia_joined","eGFR_scaled",response)
  lm_var_names_adj<-c("Diabetes","HTN","CKD","Has APOE4 Allele","Education (Yrs)","Dyslipidemia","eGFR (10 u.)","Gold Standard+")
  
  lm_vars_adj_obj_status<-c("diabetes_joined","hypertension_joined","ckd_joined","apoe_status","PTEDUCAT","dyslipidemia_joined","eGFR")
  lm_var_names_adj_obj_status<-c("Diabetes","HTN","CKD","Has APOE4 Allele","Education (Yrs)","Dyslipidemia", "eGFR (10 u.)")
  
  unadj<-univariateRegLooper(lm_vars,lm_var_names,biomarker,df=data,model_name="Unadjusted")
  adj<-univariateRegLooper(lm_vars_adj,lm_var_names_adj,biomarker,covariates=c("age_at_exam","PTGENDER"),
                           df=data,model_name="Adjusted for Age and Sex")
  adj_pet<-univariateRegLooper(lm_vars_adj_obj_status,lm_var_names_adj_obj_status,biomarker,covariates=c("age_at_exam","PTGENDER",response),
                               df=data,model_name="Adjusted for Age, Sex, and Gold Standard+")
  
  forest_data <- dplyr::bind_rows(unadj,adj,adj_pet)
  
  ## Forest plot for covariate comparison
  forest_plot<-ggplot(data=forest_data, 
                      aes(y = var_names, x = coef, color=model)) +
    geom_point(position=position_dodge(width=1)) +
    geom_errorbar(aes(xmin=lower,xmax=upper),position=position_dodge(width=1)) +
    geom_vline(xintercept=0,linetype="dashed") +
    geom_hline(yintercept=seq(1.5, length(unique(forest_data$var_names))-0.5, 1), 
               lwd=0.2, colour="black")
  
  
  print(forest_plot)
  return(forest_data)
}

## Used in replication of Mielke paper. Utility function to produce individual univariate regressions predicting biomarker levels
univariateRegLooper<- function(vars,var_names,response,covariates=NULL,df,model_name){
  if(!missing(covariates)){  
    init_formula <- paste0(paste0(response," ~ "),paste0(covariates,collapse=" + ")," + ")
  } else {
    init_formula <- paste0(paste0(response," ~ "))
  }
  models<-lapply(paste0(init_formula,vars),
                 function(frm) lm(as.formula(frm),data=df))
  num_coef<-length(covariates)+2
  coef<-unlist(purrr::map(models,
                          function(x) coef(x)[num_coef]))
  upper<-unlist(purrr::map(models,
                           function(x) confint(x)[num_coef,2]))
  lower<-unlist(purrr::map(models,
                           function(x) confint(x)[num_coef,1]))
  tibble<-dplyr::tibble(var_names,coef,lower,upper) %>% dplyr::mutate(model=model_name)
  return(tibble)
}
