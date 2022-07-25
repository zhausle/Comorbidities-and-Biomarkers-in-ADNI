## the function flags all tied records in the presence of ties; decide on ties outside of the context of this function
## function always throws a warning because giving an input with no in-range values for max returns "-Inf"; this is handled
date_comparison <- function(df,ref_date_col,comp_date_col,grouping_vars,comp_method="nearest",diff_tol_abs=NULL,diff_tol_before=NULL){
  ref_date_quo<-enquo(ref_date_col)
  comp_date_quo<-enquo(comp_date_col)
  df <- df %>% dplyr::mutate(date_diff=!!comp_date_quo-!!ref_date_quo,alt_flag=0,retain_flag=0)
  df <- df %>% dplyr::group_by(!!!grouping_vars) %>%
    dplyr::mutate(max_non_pos_diff=max(date_diff[which(date_diff<=0)]),min_abs_diff=min(abs(date_diff)),min_diff=min(date_diff))
  df$max_non_pos_diff<-ifelse(df$max_non_pos_diff=="-Inf",NA,df$max_non_pos_diff)
  if(comp_method=="nearest"){
    if(is.null(diff_tol_abs)){
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (is.na(date_diff)) ~ 1,
        (abs(date_diff)==min_abs_diff) ~ 1,
        (abs(date_diff)!=min_abs_diff) ~ 0))
    } else {
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (is.na(date_diff)) ~ 1,
        (abs(date_diff)==min_abs_diff & min_abs_diff<=diff_tol_abs) ~ 1,
        (abs(date_diff)!=min_abs_diff|min_abs_diff>diff_tol_abs) ~ 0
      ))}
  }
  if(comp_method=="before"){
    if(is.null(diff_tol_before)){
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (is.na(date_diff)) ~ 1,
        (is.na(max_non_pos_diff)) ~ 0,
        (date_diff==max_non_pos_diff) ~ 1,
        (date_diff!=max_non_pos_diff) ~ 0))
    } else {
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (is.na(date_diff)) ~ 1,
        (is.na(max_non_pos_diff)) ~ 0,
        (date_diff==min_abs_diff & max_non_pos_diff>=diff_tol_before) ~ 1,
        (date_diff!=max_non_pos_diff|date_diff<diff_tol_before) ~ 0
      ))}
    if(!is.null(diff_tol_abs)){
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(alt_flag=case_when(
        (retain_flag==0 & date_diff==min_abs_diff & date_diff<=diff_tol_abs) ~ 1,
        TRUE ~ 0
      )
      )
    }
  }
  return(df %>% select(-c(date_diff,max_non_pos_diff,min_abs_diff,min_diff)))
}

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
           condition_list,
           is_chronic=TRUE) {
    matched_records <-
      dplyr::left_join(biomarker_data, condition_list$df, by = "RID") ## copy=TRUE must be included because of how the condition looper used later processes dfs
    condition_name<-condition_list$condition
    condition_name_enquo<-quo_name(condition_list$condition)
    
    matched_records <-date_comparison(matched_records,EXAMDATE,USERDATE,vars(RID,EXAMDATE),comp_method="before") %>% dplyr::filter(retain_flag==1|alt_flag==1) %>% dplyr::select(RID,EXAMDATE,onset,cease,viscode,desc)
    
    merged_data <- dplyr::left_join(biomarker_data,matched_records,by=c("RID","EXAMDATE"))
    
    merged_data <- merged_data %>% dplyr::mutate(condition=case_when(
      (is.na(onset)) ~ 0,
      (is_chronic==FALSE) ~ 1,
      (is.na(cease) & EXAMDATE>=onset) ~ 1,
      (!is.na(cease) & EXAMDATE>=onset & EXAMDATE<cease) ~ 1,
      TRUE ~ 0
    )) %>% dplyr::rename_with(.cols=c(onset,cease,viscode,desc),.fn = ~ paste0(condition_name,"_",.x))  %>% dplyr::rename((!!condition_name_enquo):=condition)
    return(merged_data)
  }

all_conditions_merger <- function(biomarker_df,condition_df_list,chronic_list){
  df<-biomarker_df
  for (i in 1:(length(condition_df_list))){
    temp_cond<-condition_df_list[[i]]
    df<-adni_condition_merge(df,temp_cond,chronic_list[i])
    gc(reset=TRUE)
    rm(temp_cond)
    df<-df %>% distinct_at(vars(RID,VISCODE),.keep_all=TRUE)
    print(i)
  }
  return(df)
}

other_df_prep <- function(df){
  df <-
    dplyr::left_join(df, adni_lab_data_use, by = "RID")
  df <-
    dplyr::left_join(df, adni_mod_hach %>% dplyr::select(RID, htn_hach),by = "RID")
  df <-
    dplyr::left_join(
      df,
      adni_vitals %>% dplyr::select(RID, diastolic_bp, systolic_bp, htn_vitals, vitals_date),
      by = "RID")
  df<-date_comparison(df,EXAMDATE,vitals_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 120) %>% dplyr::filter(retain_flag==1|alt_flag==1|is.na(retain_flag)) %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all = TRUE)
  df$diabetes_joined <-
    ifelse(df$diabetes == 1 |
             df$diabetes_cat == "Diabetes",
           1,
           0)
  df$ckd_joined <-
    ifelse(df$ckd == 1 |
             df$ckd_cat == "Suspected CKD",
           1,
           0)
  df$dyslipidemia_joined <-
    ifelse(
      df$dyslipidemia == 1 |
        df$dyslipidemia_cat == "Hyperlipidemia",
      1,
      0
    )
  df$hypertension_joined <-
    ifelse(
      df$hypertension == 1 |
        df$htn_vitals == 1 | 
        df$htn_hach == 1,
      1,
      0
    )
  
  df <-
    dplyr::left_join(df, adni_joined_demog_uniques_reduced, by =
                       "RID")
  
  df$age_at_exam <-
    round(as.numeric(lubridate::year(df$EXAMDATE)) -
            df$PTDOBYY +
            ((
              as.numeric(lubridate::month(df$EXAMDATE)) - df$PTDOBMM
            ) / 12),
          1)
  
  df <- df %>% dplyr::mutate(age_exam_cat=cut(age_at_exam,breaks=seq(from=50,to=90,by=5)),age_exam_scaled=age_at_exam/10,eGFR_scaled=eGFR/10,age_exam_scaled_for_glmm=scale(age_at_exam))
  
  df<-dplyr::left_join(df,adni_diagnoses,by="RID")
  df<-date_comparison(df,EXAMDATE,dx_date,vars(RID,EXAMDATE),comp_method="before",diff_tol_abs = 60)
  df<- df %>% dplyr::filter(retain_flag==1|alt_flag==1|is.na(retain_flag)) %>% dplyr::distinct_at(vars(RID,EXAMDATE),.keep_all = TRUE)
  
  df<- dplyr::left_join(df,adni_apoe,by="RID")
  
  return(df)
}
