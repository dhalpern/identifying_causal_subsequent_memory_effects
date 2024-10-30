get_t_test_tbl <- function(clf_tbl, metric = "roc_auc", paired = TRUE, 
                           alternative = 'two.sided', mu = 0, ret_zero_if_const = FALSE) {
  if(paired) {
    col2 = 'glmer_metrics'
  } else {
    col2 = NULL
  }
  clf_tbl %>% unnest(c(glmer_metrics, glmnet_metrics), names_sep = "") %>%
    group_by(var, var_group, study_block, model_type, perm_id) %>%
    nest() %>%
    mutate(
      mult_comps = if_else(!(study_block %in% c('mean', 'all')), 'ind', study_block),
      roc_auc_t_test = map(data, mod_pred_t.test, metric, 
                           col1 = 'glmnet_metrics',
                           col2 = col2, paired = paired, 
                           alternative = alternative, mu = mu,
                           ret_zero_if_const = ret_zero_if_const)
    ) %>% unnest(roc_auc_t_test) %>% select(-data)
}

read_clf_file <- function(perm_id, results_dir,  prefix = 'clf_metric_df_', 
                          t_test_tbl = TRUE, use_model_type = 'causal', ret_zero_if_const = FALSE,
                          study_block_hack = TRUE) {
  print(perm_id)
  if (use_model_type == 'causal') {
    paired = TRUE
    mu = 0
  } else {
    paired = FALSE
    mu = 0.5
  }
  alternative = 'greater'
  clf_file_path <- file.path(results_dir, paste0(prefix, perm_id, '.RDS'))
  if (file.exists(clf_file_path)) {
    clf_tbl <- readRDS(clf_file_path) %>% filter(model_type == use_model_type)
    if (study_block_hack) {
      clf_tbl <- clf_tbl %>% mutate(
        study_block = if_else(is.na(study_block), 'all', study_block)
      )
    }
    clf_tbl$perm_id <- perm_id
    if (t_test_tbl) {
      clf_t_test_tbl <- get_t_test_tbl(clf_tbl, "roc_auc", alternative = alternative, 
                                       paired = paired, mu = mu, ret_zero_if_const = ret_zero_if_const)
      return(clf_t_test_tbl)
    } else {
      return(clf_tbl)
    }
  }
  else {
    return(NULL)
  }
}