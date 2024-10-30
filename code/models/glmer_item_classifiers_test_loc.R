library(tidyverse)
library(glmnet)
library(tidymodels)
library(lme4)

base_dir <- '/Users/djhalp/Documents/GitHub/identifying_causal_subsequent_memory_effects'

source(file.path(base_dir, 'code', 'models', 'model_utils.R'))

trial_word_df <- read_csv(file.path(base_dir, 'data', 'mri_subs_events_wide.csv'))
trial_word_df <- trial_word_df %>% 
  dplyr::select(lith_word, sub_id, recall_acc) %>% unique()

trial_word_df <- trial_word_df %>% mutate(
  sub_id_num = str_extract(sub_id, "[0-9]+"),
  test_location_cap = if_else(sub_id_num >= 62, "Online", "Lab"),
  test_location = if_else(sub_id_num >= 62, "online", "lab")
)

trial_word_summary <- trial_word_df %>% group_by(test_location) %>% 
  summarise(n_words = n_distinct(lith_word), n_subs = n_distinct(sub_id), 
            n = n(), mean_recall_acc = mean(recall_acc))
trial_word_summary

trial_word_df <- trial_word_df %>%
  mutate(
    recall_acc = as.factor(recall_acc)
  )

k <- 10
group_col <- 'sub_id'

set.seed(1234)

# Item Models -------------------------------------------------------------

rs_obj <- group_vfold_cv(trial_word_df, group = all_of(group_col))
ff_glmer_causal_test_loc <- recall_acc ~ 0 + lith_word + test_location + (1 | sub_id)
ff_glmer_standard_test_loc <- recall_acc ~ 0 + test_location + (1 | sub_id)
ff_glmer_causal <- recall_acc ~ 0 + lith_word + (1 | sub_id)
ff_glmer_standard <- recall_acc ~ (1 | sub_id)
# leads to singular fit
# ff_glmer <- recall_acc ~ 0 + lith_word + test_location + (1 | sub_id) + (1 | lith_word:test_location)

# gmod <- fit_glmer(ff_glmer,
#           family = binomial, data = trial_word_df)

rs_obj <- rs_obj %>% mutate(
  causal_glmer_fit = map(splits, function(splits) fit_glmer(ff_glmer_causal_test_loc, analysis(splits), family = binomial)),
  standard_glmer_fit = map(splits, function(splits) fit_glmer(ff_glmer_standard_test_loc, analysis(splits), family = binomial)),
  base_mod = 'test_loc'
)

set.seed(1234)
rs_obj <- rs_obj %>% mutate(
  analysis_df = pmap(list(splits, causal_glmer_fit, standard_glmer_fit), 
                     function(splits, causal_m, standard_m) 
                       add_glmer_offset(analysis(splits), causal_m, standard_m)),
  # foldid = map(analysis_df, function(df) df$foldid), #for updating
  # analysis_df = map2(analysis_df, foldid, function(df, foldid_col) df %>% mutate(foldid = foldid_col)), #for updating
  analysis_df = map(analysis_df, add_folds_col, group_col_str = group_col, k = k),
  assessment_df = pmap(list(splits, causal_glmer_fit, standard_glmer_fit), 
                       function(splits, causal_m, standard_m) 
                         add_glmer_offset(assessment(splits), causal_m, standard_m, re.form = ~0)),
  heldout_sub = map_chr(assessment_df, function(x) x$sub_id[[1]])
)


saveRDS(rs_obj, file.path(base_dir, 'data', 'irt_mod_fit_df_test_loc.RDS'))
