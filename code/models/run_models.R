library(tidyverse)
library(glmnet)
library(glmnetUtils)
library(tidymodels)
library(lme4)

base_dir <- '/Users/davidhalpern/Documents/GitHub/identifying_causal_subsequent_memory_effects'
if(dir.exists(base_dir)) {
  args <- c(1, 'irt')
  laptop <- TRUE
  # args = commandArgs(trailingOnly=TRUE)
  data_dir <- file.path(base_dir, 'data')
  code_dir <- file.path(base_dir, 'code', 'models')
  output_dir <- file.path(base_dir, 'code', 'models')
} else {
  data_dir <- file.path('/scratch', 'djh389', 'fmri_data')
  code_dir <- file.path('/home', 'djh389', 'identifying_causal_subsequent_memory_effects', 'models')
  output_dir <- file.path('/scratch', 'djh389', 'classifier_results')
  args = commandArgs(trailingOnly=TRUE)
  laptop <- FALSE
}

arr_ind <- as.numeric(args[1])
base_mod <- args[2]

source(file.path(code_dir, 'model_utils.R'))
print('standardizing within subject!')

raw_classifier_df <- read_csv(file.path(data_dir, 'classifier_rep_df.csv'),
                          col_types = cols(roi = col_character(), study_block = col_character()))
mean_classifier_df <- raw_classifier_df %>% 
  group_by(sub_id, roi, lith_word, var, var_group) %>% 
  summarise(value = mean(value)) %>% mutate(study_block = "mean")
classifier_df <- bind_rows(raw_classifier_df, mean_classifier_df)

classifier_df <- classifier_df %>% 
  group_by(study_block, sub_id, roi, var, var_group) %>% mutate(
  value = scale2(value, na.rm = TRUE)
) %>% ungroup()

rs_obj <- readRDS(file.path(data_dir, 'irt_mod_fit_df_test_loc.RDS'))

events_df <- read_csv(file.path(data_dir, 'mri_subs_events_wide.csv'))
trial_word_df <- events_df %>% 
  dplyr::select(lith_word, sub_id, recall_acc) %>% distinct()

jol_df <- events_df %>% 
  dplyr::select(lith_word, sub_id, jol) %>% distinct()

classifier_df_wide <- classifier_df %>%
  pivot_wider(names_from = c(var_group, var, study_block, roi), values_from = value) %>% 
  mutate(
    fake_intercept = 1
    ) %>% left_join(jol_df, by = c("lith_word", "sub_id"))

set.seed(arr_ind)
if(arr_ind > 1) {
  #shuffle words within block
  lith_word_perm_df <- trial_word_df %>% 
    distinct(lith_word) %>% 
    mutate(
      block = rep(1:9, each = 5)
    ) %>% group_by(block) %>%
    mutate(
      lith_word_perm = sample(lith_word)
    )
  
  #switch lith_word to permuted version
  classifier_df_wide <- classifier_df_wide %>% 
    left_join(lith_word_perm_df, by = "lith_word") %>%
    mutate(
      lith_word = lith_word_perm
    )
}

classifier_rs_df <- rs_obj %>% mutate(
  analysis_df = map(analysis_df, left_join, classifier_df_wide, by = c("lith_word", "sub_id")),
  assessment_df = map(assessment_df, left_join, classifier_df_wide, by = c("lith_word", "sub_id"))
)

# Get formulas

clf_feature_df <- classifier_df %>% select(var_group, var, study_block, roi) %>% 
  distinct() %>% mutate(
    feature = paste(var_group, var, study_block, roi, sep = "_")
  ) %>% select(-roi)

clf_ind_study_block_formula_df <- clf_feature_df %>%
  group_by(var_group, var, study_block) %>% nest() %>%
  mutate(
    formulas = map(data, ~ reformulate(.x$feature, response = "recall_acc"))
  ) %>% select(-data)

clf_all_study_block_formula_df <- clf_feature_df %>%
  filter(study_block != "mean") %>%
  group_by(var_group, var) %>% nest() %>%
  mutate(
    formulas = map(data, ~ reformulate(.x$feature, response = "recall_acc")),
    study_block = "all"
  ) %>% select(-data)


#features across var groups
clf_all_feature_df <- clf_feature_df %>%
  mutate(
    study_block = if_else(study_block == "mean", study_block, "all")
    )

clf_all_feature_formula_df <- clf_all_feature_df %>% 
  group_by(study_block) %>%
  nest() %>%
  mutate(
    formulas = map(data, ~ reformulate(.x$feature, response = "recall_acc")),
    var = "all",
    var_group = "all"
  ) %>% select(-data)

clf_all_var_group_ind_study_block_formula_df <- clf_feature_df %>%
  group_by(var_group, study_block) %>% nest() %>%
  mutate(
    formulas = map(data, ~ reformulate(.x$feature, response = "recall_acc")),
    var = "all"
  ) %>% select(-data)

clf_all_var_group_all_study_block_formula_df <- clf_feature_df %>%
  filter(study_block != "mean") %>%
  group_by(var_group) %>% nest() %>%
  mutate(
    formulas = map(data, ~ reformulate(.x$feature, response = "recall_acc")),
    study_block = "all",
    var = "all"
  ) %>% select(-data)

clf_jol_formula_df <- tibble(var = "jol", var_group = "beh", study_block = "all", formulas = list(recall_acc ~ jol))

clf_formula_df <- bind_rows(clf_all_feature_formula_df,
                            clf_jol_formula_df, clf_all_study_block_formula_df, clf_ind_study_block_formula_df)

clf_formula_df <- clf_formula_df %>% expand_grid(model_type = c('causal', 'standard'))

clf_full_res_df <- clf_formula_df %>% 
  mutate(
    mod_fit_df = map2(formulas, model_type, fit_cvglmnet, 
                      rs_obj = classifier_rs_df, alpha = 0, s = "lambda.min", 
                      standardize = FALSE)
  )

clf_res_df <- clf_full_res_df %>% unnest(mod_fit_df) %>% unnest(pred_df) %>% 
  select(-cvglmnet_fit, -heldout_sub, -formulas)

multi_metric <- metric_set(roc_auc, mn_log_loss)
clf_metric_df <- clf_res_df %>%
  group_by(var, var_group, study_block, sub_id, model_type) %>% nest() %>% 
  mutate(
    glmer_metrics = map(data, multi_metric, factor(recall_acc, levels = c(1, 0)), glmer_preds),
    glmnet_metrics = map(data, multi_metric, factor(recall_acc, levels = c(1, 0)), glmnet_preds)
  )

if(arr_ind < 3) {
  saveRDS(clf_full_res_df, file.path(output_dir, paste0('full_standardized_res_df_', base_mod, '_', arr_ind, '.RDS')))
}
saveRDS(clf_res_df, file.path(output_dir, paste0('standardized_res_df_', base_mod, '_',  arr_ind, '.RDS')))
saveRDS(clf_metric_df, file.path(output_dir, paste0('standardized_metric_df_', base_mod, '_',  arr_ind, '.RDS')))




