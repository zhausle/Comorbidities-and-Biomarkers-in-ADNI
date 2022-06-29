## the function flags all tied records in the presence of ties; decide on ties outside of the context of this function
date_comparison <- function(df,ref_date_col,comp_date_col,grouping_vars,comp_method="nearest",diff_tol_abs=NULL,diff_tol_before=NULL){
  ref_date_quo<-enquo(ref_date_col)
  comp_date_quo<-enquo(comp_date_col)
  df <- df %>% dplyr::mutate(date_diff=!!ref_date_quo-!!comp_date_quo,alt_flag=0,retain_flag=0)
  df <- df %>% dplyr::group_by(!!!grouping_vars) %>%
    dplyr::mutate(max_non_pos_diff=max(date_diff[which(date_diff<=0)]),min_abs_diff=min(abs(date_diff)),min_diff=min(date_diff))
  df$max_non_pos_diff<-ifelse(df$max_non_pos_diff=="-Inf",NA,df$max_non_pos_diff)
  if(comp_method=="nearest"){
    if(is.null(diff_tol_abs)){
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (abs(date_diff)==min_abs_diff) ~ 1,
        (abs(date_diff)!=min_abs_diff) ~ 0))
    } else {
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (abs(date_diff)==min_abs_diff & min_abs_diff<=diff_tol_abs) ~ 1,
        (abs(date_diff)!=min_abs_diff|min_abs_diff>diff_tol_abs) ~ 0
      ))}
  }
  if(comp_method=="before"){
    if(is.null(diff_tol_before)){
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
        (is.na(max_non_pos_diff)) ~ 0,
        (date_diff==max_non_pos_diff) ~ 1,
        (date_diff!=max_non_pos_diff) ~ 0))
    } else {
      df <- df %>% dplyr::group_by(!!!grouping_vars) %>% dplyr::mutate(retain_flag=case_when(
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