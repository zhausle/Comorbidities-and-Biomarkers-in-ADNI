library(lme4)
library(lmerTest)

biomarker_response_join<-function(biomarker_data,response_data,biomarker,response){
  
  biomarker_enquo=enquo(biomarker)
  response_enquo=enquo(response)
  
  joined_data<-dplyr::left_join(biomarker_data, response_data %>% dplyr::select(RID,EXAMDATE,!!response_enquo) %>% rename(response_date=EXAMDATE),by="RID") 
  joined_data<-date_comparison(joined_data,EXAMDATE,response_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60) %>% dplyr::filter(retain_flag==1|alt_flag==1)
  joined_data<-joined_data %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE)
  joined_data<-joined_data %>% dplyr::ungroup() %>% dplyr::mutate(biomarker_z_scaled=((!!biomarker_enquo-mean(!!biomarker_enquo,na.rm=TRUE)))/sd(!!biomarker_enquo,na.rm=TRUE))
  joined_data <- joined_data %>% dplyr::mutate(dyslipidemia_cat=relevel(as.factor(dyslipidemia_cat),ref="Normal Cholesterol"),
                                               diabetes_cat=relevel(as.factor(diabetes_cat),ref="Normoglycemic"),
                                               ckd_cat=relevel(as.factor(ckd_cat),ref="No Suspected CKD"),
                                               eGFR_cat=relevel(as.factor(eGFR_cat),ref="Less than 45"))
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

mielke_analysis_roc_and_measure_inner <- function(formula,data,longitudinal=FALSE){
  if(longitudinal){
  model<-lme4::glmer(paste0(formula,"+ (1|RID)"),family=binomial,data=data,glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
  roc<-ROCit::rocit(score=predict(model,re.form=NA),class=model@resp$`.->y`)
  measure<-ROCit::measureit(score=predict(model,re.form=NA),class=model@resp$`.->y`,
                            measure=c("ACC","SENS","SPEC","PPV","NPV"))
  } else {
  model<-glm(formula,family=binomial,data=data)
  roc<-ROCit::rocit(score=qlogis(model$fitted.values),class=model$y)
  measure<-ROCit::measureit(score=predict(model,re.form=NA),class=model$y,
                            measure=c("ACC","SENS","SPEC","PPV","NPV"))}
  measure_table<-cbind.data.frame(measure$Cutoff,measure$ACC,
                                  measure$SENS, measure$SPEC,
                                  measure$PPV,measure$NPV) %>% dplyr::mutate(youden=(measure$SENS+measure$SPEC-1))
  colnames(measure_table)<-c("cutoff","ACC","SENS","SPEC","PPV","NPV","Youden")
  results<-list(model,roc,measure,measure_table)
  plot(roc)
  return(results)
}

mielke_analysis_roc_and_measure_outer <- function(data,biomarker,response,longitudinal=FALSE){
  unadj_formula<-paste0(response," ~ ",biomarker)
  adj_formula<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",sep=" + "))
  adj_formula_apoe<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER","apoe_status",sep=" + "))
  adj_formula_w_comorbidities<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                           "dyslipidemia_joined","ckd_joined","hypertension_joined","diabetes_joined","apoe_status",sep=" + "))
  adj_formula_w_comorbidities_ethnicity<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                           "dyslipidemia_joined","ckd_joined","hypertension_joined","diabetes_joined","apoe_status","PTETHCAT",sep=" + "))
  if(longitudinal){
  adj_formula_w_comorbidities_labs<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                                "dyslipidemia_cat","creatinine","htn_vitals","diabetes_cat","apoe_status",sep=" + "))
  } else {
  adj_formula_w_comorbidities_labs<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                                "dyslipidemia_cat","ckd_cat","htn_vitals","diabetes_cat","apoe_status",sep=" + "))}
  adj_formula_w_comorbidities_hist<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                                "dyslipidemia","ckd","hypertension","diabetes","apoe_status",sep=" + "))
  adj_formula_w_lab_levels<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                                "cholesterol","systolic_bp","diastolic_bp","glucose","creatinine","apoe_status",sep=" + "))
  adj_formula_w_comorbidities_no_apoe<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER",
                                                           "dyslipidemia_joined","ckd_joined","hypertension_joined","diabetes_joined",sep=" + "))
  unadj<-mielke_analysis_roc_and_measure_inner(unadj_formula,data,longitudinal)
  adj<-mielke_analysis_roc_and_measure_inner(adj_formula,data,longitudinal)
  adj_apoe<-mielke_analysis_roc_and_measure_inner(adj_formula_apoe,data,longitudinal)
  adj_comorbidities<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities,data,longitudinal)
  adj_comorbidities_labs<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities_labs,data,longitudinal)
  adj_comorbidities_hist<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities_hist,data,longitudinal)
  adj_comorbidities_ethnicity<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities_ethnicity,data,longitudinal)
  adj_labs<-mielke_analysis_roc_and_measure_inner(adj_formula_w_lab_levels,data,longitudinal)
  adj_comorbidities_no_apoe<-mielke_analysis_roc_and_measure_inner(adj_formula_w_comorbidities_no_apoe,data,longitudinal)
  results<-list(unadj,adj,adj_apoe,adj_comorbidities,adj_comorbidities_labs,adj_comorbidities_ethnicity,adj_comorbidities_hist,adj_labs,adj_comorbidities_no_apoe)
  results<-set_names(results,c("Unadjusted","Adjusted for Age & Sex","Adjusted for Age, Sex, and APOE",
                               "Adjusted for Age, Sex, APOE, and Comorbidities","Adjusted for Age, Sex, APOE, Comorbidities, and Ethnicity",
                               "Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)","Adjusted for Age, Sex, APOE, and Comorbidities (History Only)","Adjusted for Age, Sex, APOE, and Lab Levels", "Adjusted for Age, Sex, and Comorbidities"))
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

mielke_joined_roc_plots<-function(data,title){
  plot(data$Unadjusted[[2]], col = c(1,"gray50"), 
       legend = FALSE, YIndex = FALSE)
  lines(data$`Adjusted for Age & Sex`[[2]]$TPR ~ data$`Adjusted for Age & Sex`[[2]]$FPR, 
        col = 2, lwd = 2)
  lines(data$`Adjusted for Age, Sex, and APOE`[[2]]$TPR ~ data$`Adjusted for Age, Sex, and APOE`[[2]]$FPR, 
        col = 3, lwd = 2)
  lines(data$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$TPR ~ data$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$FPR, 
        col = 4, lwd = 2)
  lines(data$`Adjusted for Age, Sex, APOE, Comorbidities, and Ethnicity`[[2]]$TPR ~ data$`Adjusted for Age, Sex, APOE, Comorbidities, and Ethnicity`[[2]]$FPR, 
        col = 5, lwd = 2)
  lines(data$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$TPR ~ data$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$FPR, 
        col = 6, lwd = 2)
  lines(data$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$TPR ~ data$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$FPR, 
        col = 7, lwd = 2)
  lines(data$`Adjusted for Age, Sex, APOE, and Lab Levels`[[2]]$TPR ~ data$`Adjusted for Age, Sex, APOE, and Lab Levels`[[2]]$FPR, 
        col = 8, lwd = 2)
  lines(data$`Adjusted for Age, Sex, and Comorbidities`[[2]]$TPR ~ data$`Adjusted for Age, Sex, and Comorbidities`[[2]]$FPR, 
        col = 9, lwd = 2)
  legend("bottomright",col = c(1,2,3,4,5,6,7,8,9),
         c(paste0("Unadjusted: AUC = ", round(data$Unadjusted[[2]]$AUC,digits=3)),
           paste0("Age & Sex: AUC = ", round(data$`Adjusted for Age & Sex`[[2]]$AUC,digits=3)),
           paste0("Age, Sex & APOE: AUC = ", round(data$`Adjusted for Age, Sex, and APOE`[[2]]$AUC,digits=3)),
           paste0("Age, Sex, APOE & Comorbidities: AUC = ", round(data$`Adjusted for Age, Sex, APOE, and Comorbidities`[[2]]$AUC,digits=3)), 
           paste0("Age, Sex, APOE, Comorbidities & Ethnicity: AUC = ", round(data$`Adjusted for Age, Sex, APOE, Comorbidities, and Ethnicity`[[2]]$AUC,digits=3)), 
           paste0("Age, Sex, APOE & Comorbidities (Labs Only): AUC = ", round(data$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[2]]$AUC,digits=3)),
           paste0("Age, Sex, APOE & Comorbidities (History Only): AUC = ", round(data$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[2]]$AUC,digits=3)), 
           paste0("Age, Sex, APOE & Lab Levels: AUC = ", round(data$`Adjusted for Age, Sex, APOE, and Lab Levels`[[2]]$AUC,digits=3)), 
           paste0("Age, Sex & Comorbidities: AUC = ", round(data$`Adjusted for Age, Sex, and Comorbidities`[[2]]$AUC,digits=3))), 
         lwd = 2,cex=0.7)
  title(main = title)
}

mielke_joined_logistic_measures <- function(data,title) {
  results<-dplyr::bind_rows(data$Unadjusted[[4]][which.max(data$Unadjusted[[4]]$Youden),] %>% dplyr::mutate(Model = "Unadjusted"),
                            data$`Adjusted for Age & Sex`[[4]][which.max(data$`Adjusted for Age & Sex`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age & Sex"),
                            data$`Adjusted for Age, Sex, and APOE`[[4]][which.max(data$`Adjusted for Age, Sex, and APOE`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex, and APOE"),
                            data$`Adjusted for Age, Sex, APOE, and Comorbidities`[[4]][which.max(data$`Adjusted for Age, Sex, APOE, and Comorbidities`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex, APOE, and Comorbidities"),
                            data$`Adjusted for Age, Sex, APOE, Comorbidities, and Ethnicity`[[4]][which.max(data$`Adjusted for Age, Sex, APOE, Comorbidities, and Ethnicity`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex, APOE, Comorbidities, and Ethnicity"),
                            data$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[4]][which.max(data$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex, APOE, and Comorbidities (Labs Only)"),
                            data$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[4]][which.max(data$`Adjusted for Age, Sex, APOE, and Comorbidities (History Only)`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex, APOE, and Comorbidities (History Only)"),
                            data$`Adjusted for Age, Sex, APOE, and Lab Levels`[[4]][which.max(data$`Adjusted for Age, Sex, APOE, and Lab Levels`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex, APOE, and Lab Levels"),
                            data$`Adjusted for Age, Sex, and Comorbidities`[[4]][which.max(data$`Adjusted for Age, Sex, and Comorbidities`[[4]]$Youden),] %>% dplyr::mutate(Model = "Age, Sex and Comorbidities"),
  ) %>%
    dplyr::relocate(Model)
  flextable::flextable(results) %>% flextable::autofit() %>% flextable::add_header_lines(value=title) %>% vline(j = 1) %>% colformat_double(digits=3)
}

mielke_coef_table_fn<-function(model,title){
  flextable::flextable(as.data.frame(summary(model)$coefficients) %>% tibble::rownames_to_column(var="Predictor")) %>% flextable::autofit() %>% flextable::colformat_double(digits=3) %>% flextable::add_header_lines(value=title)
}

mielke_roc_mega_wrapper<-function(data,biomarker,response,base_title,longitudinal=FALSE){
  roc_and_measure<-mielke_analysis_roc_and_measure_outer(data,biomarker,response,longitudinal)
  mielke_joined_roc_plots(roc_and_measure,base_title)
  print(mielke_joined_logistic_measures(roc_and_measure,base_title))
  print(mielke_coef_table_fn(roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities`[[1]],paste0(base_title," Full Model")))
  print(mielke_coef_table_fn(roc_and_measure$`Adjusted for Age, Sex, APOE, and Comorbidities (Labs Only)`[[1]],paste0(base_title,", Full Model, Lab Comorbidities Only")))
  return(roc_and_measure)
}  

mielke_analysis_roc_and_measure_comparisons <- function(data,biomarker,response,longitudinal=FALSE){
  adj_formula_apoe<-paste0(response," ~ ",paste(biomarker,"age_exam_scaled_for_glmm","PTGENDER","apoe_status",sep=" + "))
  adj_apoe_all_obs<-mielke_analysis_roc_and_measure_inner(adj_formula_apoe,data,longitudinal)
  adj_apoe_no_ckd<-mielke_analysis_roc_and_measure_inner(adj_formula_apoe,data %>% dplyr::filter(ckd_joined==0),longitudinal)
  adj_apoe_ckd_only<-mielke_analysis_roc_and_measure_inner(adj_formula_apoe,data %>% dplyr::filter(ckd_joined==1),longitudinal)
  adj_apoe_all_comorbidities<-mielke_analysis_roc_and_measure_inner(adj_formula_apoe,data %>% dplyr::filter(ckd_joined==1|dyslipidemia_joined==1|hypertension_joined==1|diabetes_joined==1),longitudinal)
  results<-list(adj_apoe_all_obs,adj_apoe_no_ckd,adj_apoe_ckd_only,adj_apoe_all_comorbidities)
  results<-set_names(results,c("All Observations","No CKD","CKD Only","All Comorbidities"))
  return(results)}

mielke_joined_roc_plots_comparisons<-function(data,title){
  plot(data$`All Observations`[[2]], col = c(1,"gray50"), 
       legend = FALSE, YIndex = FALSE)
  lines(data$`No CKD`[[2]]$TPR ~ data$`No CKD`[[2]]$FPR, 
        col = 2, lwd = 2)
  lines(data$`CKD Only`[[2]]$TPR ~ data$`CKD Only`[[2]]$FPR, 
        col = 3, lwd = 2)
  lines(data$`All Comorbidities`[[2]]$TPR ~ data$`All Comorbidities`[[2]]$FPR, 
        col = 4, lwd = 2)
  legend("bottomright",col = c(1,2,3,4),
         c(paste0("All Observations: AUC = ", round(data$`All Observations`[[2]]$AUC,digits=3)),
           paste0("No CKD: AUC = ", round(data$`No CKD`[[2]]$AUC,digits=3)),
           paste0("CKD Only: AUC = ", round(data$`CKD Only`[[2]]$AUC,digits=3)),
           paste0("All Comorbidities: AUC = ", round(data$`All Comorbidities`[[2]]$AUC,digits=3))), 
         lwd = 2,cex=0.7)
  title(main = title)
}

mielke_joined_logistic_measures_comparisons <- function(data,title) {
  results<-dplyr::bind_rows(
    data$`All Observations`[[4]][which.max(data$`All Observations`[[4]]$Youden),] %>% dplyr::mutate(Dataset = "All Observations"),
    data$`No CKD`[[4]][which.max(data$`No CKD`[[4]]$Youden),] %>% dplyr::mutate(Dataset = "No CKD"),
    data$`CKD Only`[[4]][which.max(data$`CKD Only`[[4]]$Youden),] %>% dplyr::mutate(Dataset = "CKD Only"),
    data$`All Comorbidities`[[4]][which.max(data$`All Comorbidities`[[4]]$Youden),] %>% dplyr::mutate(Dataset = "All Comorbidities")
  )
#    dplyr::relocate(Dataset)
  flextable::flextable(results) %>% flextable::autofit() %>% flextable::add_header_lines(value=title) %>% vline(j = 1) %>% colformat_double(digits=3)
}
