library(tidyverse)
library(tidymodels)

base_dir <- '/Users/davidhalpern/Documents/GitHub/identifying_causal_subsequent_memory_effects'
fig_dir <- file.path(base_dir, 'doc/figs')
results_dir <- file.path(base_dir, 'results')
test_loc_str <- 'test_loc_'

source(file.path(base_dir, 'code', 'models', 'model_utils.R'))
source(file.path(base_dir, 'code', 'models', 'load_utils.R'))

clf_perm_t_test_df <- c(2:502) %>%
  map_dfr(read_clf_file, results_dir = results_dir,
          prefix = paste0('standardized_metric_df_irt_', test_loc_str), 
          use_model_type = 'causal', ret_zero_if_const = TRUE)
write_csv(clf_perm_t_test_df, file.path(results_dir, paste0(test_loc_str, 
                                                            'clf_perm_causal_one_tailed_t_test_df.csv')))