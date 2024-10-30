library(tidyverse)
library(glmnet)
library(tidymodels)
library(lme4)

base_dir <- '/Users/djhalp/Documents/GitHub/identifying_causal_subsequent_memory_effects'
fig_dir <- file.path(base_dir, 'doc/figs')
results_dir <- file.path(base_dir, 'results')

source(file.path(base_dir, 'code', 'models', 'model_utils.R'))

test_loc_str <- 'test_loc_'
# plot_model_type <- 'causal'
plot_model_type <- 'standard'

if(plot_model_type == "standard") {
  star_loc_col <- 'star_loc'
  estimate_col <- 'estimate'
} else {
  star_loc_col <- 'star_loc_causal'
  estimate_col <- 'estimate2'
}

clf_metric_df_1 <- readRDS(file.path(results_dir, 
                                     paste0("standardized_metric_df_irt_", test_loc_str, "1.RDS")))
# clf_metric_df_1_old <- readRDS(file.path(results_dir_old, "clf_standardized_metric_df_1.RDS"))

# clf_full_df_1 <- readRDS(file.path(results_dir, "clf_full_standardized_res_df_1.RDS"))
# clf_full_df_1_old <- readRDS(file.path(results_dir_old, "clf_full_standardized_res_df_1.RDS"))
# 
# t1 <- clf_full_df_1 %>% filter(model_type == 'standard')
# t1_old <- clf_full_df_1_old %>% filter(use_offset == FALSE)

classifier_metric_df <- clf_metric_df_1 %>% mutate(
  var_group = recode(var_group, whole_hipp = "hipp_whole"),
  study_block = if_else(is.na(study_block), 'all', study_block) #hack, possibly remove
)
classifier_metric_nest_df <- classifier_metric_df %>% unnest(c(glmer_metrics, glmnet_metrics), names_sep = "") %>%
  group_by(var, var_group, study_block, model_type) %>%
  nest()

classifier_model_type_metric_nest_df <- classifier_metric_nest_df %>%
  filter(model_type == plot_model_type)

mod_pred_t.test_model_type <- function(data, metric, model_type) {
  if (model_type == "standard") {
    mod_pred_t.test(data, metric, 
                    col1 = 'glmnet_metrics', col2 = NULL,
                    paired = FALSE, alternative = 'greater', mu = .5)
  } else {
    mod_pred_t.test(data, metric, col1 = 'glmnet_metrics',
    col2 = 'glmer_metrics', paired = TRUE, alternative = 'greater')
  }
} 

classifier_model_type_ttest_df <- classifier_model_type_metric_nest_df %>% 
  mutate(
    mult_comps = if_else(!(study_block %in% c('mean', 'all')), 'ind', study_block),
    mn_log_loss_t_test = map(data, mod_pred_t.test_model_type, "mn_log_loss", 
                             plot_model_type),
    roc_auc_t_test = map(data, mod_pred_t.test_model_type, "roc_auc", 
                         plot_model_type)
  ) #%>% filter(mult_comps != 'ind')

clf_perm_t_test_df <- read_csv(file.path(results_dir, 
                                         paste0(test_loc_str, 
                                                'clf_perm_', 
                                                plot_model_type,
                                                '_one_tailed_t_test_df.csv'))) 
clf_perm_t_test_df <- clf_perm_t_test_df %>% mutate(
  var_group = recode(var_group, whole_hipp = "hipp_whole")
)

clf_t_test_df_1 <- classifier_model_type_ttest_df %>% unnest(roc_auc_t_test) %>% 
  select(-data, -mn_log_loss_t_test) %>%
  mutate(perm_id = 1)
clf_auc_t_df <- bind_rows(clf_perm_t_test_df, clf_t_test_df_1)
clf_auc_t_perm_df <- clf_auc_t_df %>%
  group_by(var_group, var, study_block, mult_comps) %>%
  mutate(
    one_tailed_t_rank = rank(statistic),
    two_tailed_t_rank = rank(abs(statistic)),
    max_rank = max(one_tailed_t_rank),
    perm_p_value = 1 - (one_tailed_t_rank / (max_rank))
  )  %>% 
  filter(perm_id == 1) %>% ungroup %>% 
  select(var, var_group, study_block, mult_comps, estimate, perm_p_value, one_tailed_t_rank, max_rank) %>%
  mutate(
    perm_p_value = if_else(estimate == 0, 1, perm_p_value) 
    #some models have all permutation estimates = 0, these should have a p-value of 1 but get ranked first
  )

clf_auc_t_perm_df <- clf_auc_t_perm_df %>%
  group_by(mult_comps) %>%
  mutate(
    perm_q = p.adjust(perm_p_value, method = "fdr")
  )

plot_df <- classifier_model_type_metric_nest_df %>% unnest(cols = c(data)) %>%
  filter(glmer_metrics.metric == "roc_auc") %>% 
  dplyr::select(var_group, var, study_block, ends_with(".estimate")) %>% #select(heldout_sub, ends_with(".estimate")) %>%
  pivot_longer(cols = ends_with(".estimate"), 
               names_to = "model", 
               names_pattern = "(.*)_metrics.estimate",
               values_to = "roc_auc")

item_plot_df <- plot_df %>% 
  filter(model == "glmer", var == "ISC_z", var_group == "schaefer", study_block == "all") %>%
  mutate(
    var = "item",
    var_group = "beh"
  )
jol_plot_df <- plot_df %>% filter(var == "jol", model == "glmnet")
beh_plot_df <- item_plot_df %>% bind_rows(jol_plot_df) %>%
  ungroup() %>%
  mutate(
    study_block = NULL,
    var_group = "Behavior",
    var = recode(var,
                 jol = "JOL",
                 item = "IRT"),
    feature = var,
    var_group = factor(var_group, levels = c("Behavior",
                                             "All",
                                             "Hippocampus-Whole",
                                             "Hippocampus-Subparts",
                                             "Hippocampus-Voxels",
                                             "Schaefer",
                                             "Targeted",
                                             "Whole Brain (Schaefer ROIs)"))
  )
mri_plot_df <- plot_df %>% filter(model == "glmnet", var != "jol") %>% 
  ungroup() %>%
  arrange(var_group, var) %>%
  mutate(
    var = recode(var,
                 all = "",
                 pattern_sim_z_Within = "IPS",
                 pattern_sim_z_Cross = "GPS",
                 ISC_z = "ISPC"),
    feature = var,
    var_group = recode(var_group, 
                       hipp_whole = "Hippocampus-Whole",
                       hipp_sub = "Hippocampus-Subparts",
                       hipp_vox = "Hippocampus-Voxels",
                       schaefer = "Schaefer",
                       tr = "Targeted",
                       whole_brain_schaefer = "Whole Brain (Schaefer ROIs)",
                       all = "All"),
    var_group = factor(var_group, levels = c("Behavior", 
                                            "All",
                                            "Hippocampus-Whole",
                                            "Hippocampus-Subparts",
                                            "Hippocampus-Voxels",
                                            "Schaefer",
                                            "Targeted",
                                            "Whole Brain (Schaefer ROIs)"))
  ) 

mri_plot_df_subset <- mri_plot_df %>% filter(study_block %in% c('mean', 'all'))

perm_plot_df <- clf_perm_t_test_df %>% 
  mutate(
    var = recode(var,
                 all = "",
                 pattern_sim_z_Within = "IPS",
                 pattern_sim_z_Cross = "GPS",
                 jol = "JOL",
                 item = "IRT",
                 ISC_z = "ISPC"),
    feature = var,
    var_group = recode(var_group, 
                       beh = "Behavior",
                       hipp_whole = "Hippocampus-Whole",
                       hipp_sub = "Hippocampus-Subparts",
                       hipp_vox = "Hippocampus-Voxels",
                       schaefer = "Schaefer",
                       tr = "Targeted",
                       whole_brain_schaefer = "Whole Brain (Schaefer ROIs)",
                       all = "All"),
    var_group = factor(var_group, levels = c("Behavior", 
                                             "All",
                                             "Hippocampus-Whole",
                                             "Hippocampus-Subparts",
                                             "Hippocampus-Voxels",
                                             "Schaefer",
                                             "Targeted",
                                             "Whole Brain (Schaefer ROIs)")),
    estimate2 = estimate + .72
  )

perm_plot_sig_df <- clf_auc_t_perm_df %>%
  mutate(
    var = recode(var,
                 all = "",
                 pattern_sim_z_Within = "IPS",
                 pattern_sim_z_Cross = "GPS",
                 jol = "JOL",
                 item = "IRT",
                 ISC_z = "ISPC"),
    feature = var,
    var_group = recode(var_group,
                       beh = "Behavior",
                       hipp_whole = "Hippocampus-Whole",
                       hipp_sub = "Hippocampus-Subparts",
                       hipp_vox = "Hippocampus-Voxels",
                       schaefer = "Schaefer",
                       tr = "Targeted",
                       whole_brain_schaefer = "Whole Brain (Schaefer ROIs)", 
                       all = "All"),
    var_group = factor(var_group, levels = c("Behavior", 
                                             "All",
                                             "Hippocampus-Whole",
                                             "Hippocampus-Subparts",
                                             "Hippocampus-Voxels",
                                             "Schaefer",
                                             "Targeted",
                                             "Whole Brain (Schaefer ROIs)")),
    sig = if_else(perm_q < 0.05, "**", if_else(perm_p_value < 0.05, "*", "")),
    star_loc = .66,
    estimate2 = estimate + .72,
    star_loc_causal = estimate2 + .01,
  )

beh_perm_plot_df <- perm_plot_df %>% filter(var_group == "Behavior", mult_comps != "ind") %>% mutate(study_block = NULL)
mri_perm_plot_df <- perm_plot_df %>% filter(var_group != "Behavior", mult_comps != "ind")
beh_perm_plot_sig_df <- perm_plot_sig_df %>% filter(var_group == "Behavior", mult_comps != "ind") %>% mutate(study_block = NULL)
mri_perm_plot_sig_df <- perm_plot_sig_df %>% filter(var_group != "Behavior", mult_comps != "ind")

# permutation plot
p <- ggplot(mri_plot_df_subset, aes(y = feature, x = roc_auc, color = var_group)) + 
  geom_violin(data = mri_perm_plot_df,
              aes(y = feature, x = .data[[estimate_col]]), color = "grey") +
  stat_summary() + 
  geom_text(data = mri_perm_plot_sig_df,
            aes(y = feature, x = .data[[star_loc_col]], label = sig)) +
  facet_grid(cols = vars(study_block), rows = vars(var_group), 
             scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_blank(), # Make facet label background white.
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.text.y.left = element_text(angle = 0, hjust = 1)
  ) + xlab("ROC AUC") +
  scale_x_continuous(n.breaks = 3) +
  scale_colour_viridis_d()

if(plot_model_type == "standard") {
  p <- p + geom_vline(xintercept = .5)
} else {
  p <- p + geom_violin(data = beh_perm_plot_df,
              aes(y = feature, x = .data[[estimate_col]]), color = "grey") +
    stat_summary(data = beh_plot_df) +
    geom_text(data = beh_perm_plot_sig_df,
              aes(y = feature, x = .data[[star_loc_col]], label = sig))
}
p
ggsave(file.path(fig_dir, paste0(test_loc_str, plot_model_type,
                                 '_clf_plot.pdf')), width = 6, height = 6, plot = p)


# Study-block plot --------------------------------------------------------

beh_plot_df2 <- beh_plot_df %>% mutate(feature2 = NULL, var = '')
mri_plot_df_subset2 <- mri_plot_df %>% filter(!(study_block %in% c('mean', 'all'))) %>% 
  mutate(feature = factor(study_block), feature2 = var_group)

beh_perm_plot_df2 <- perm_plot_df %>% filter(var_group == "Behavior", mult_comps == "ind") %>% mutate(study_block = NULL)
mri_perm_plot_df2 <- perm_plot_df %>% filter(var_group != "Behavior", mult_comps == "ind") %>% mutate(feature2 = var_group)
beh_perm_plot_sig_df2 <- perm_plot_sig_df %>% filter(var_group == "Behavior", mult_comps == "ind") %>% mutate(study_block = NULL)
mri_perm_plot_sig_df2 <- perm_plot_sig_df %>% filter(var_group != "Behavior", mult_comps == "ind") %>% mutate(feature2 = var_group)

p <- ggplot(mri_plot_df_subset2, aes(y = feature, x = roc_auc, color = var_group)) + 
  stat_summary() + 
  geom_text(data = mri_perm_plot_sig_df2,
            aes(y = study_block, x = .data[[star_loc_col]], label = sig)) +
  facet_grid(cols = vars(feature2), rows = vars(var),
             scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() + 
  theme(strip.placement = "outside",                      # Place facet labels outside x axis labels.
        strip.background = element_blank(),# Make facet label background white.
        # strip.text.x = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        legend.position = "none",
        #axis.title.y = element_blank()
  ) + xlab("ROC AUC") + ylab("Study block/Study block pair (by feature)") +
  scale_x_continuous(n.breaks = 3) +
  scale_colour_viridis_d()

if(plot_model_type == "standard") {
  p <- p + geom_vline(xintercept = .5)
} else {
  p <- p + stat_summary(data = beh_plot_df2)
}
p

ggsave(file.path(fig_dir, paste0(test_loc_str, plot_model_type, '_study_block_clf_plot.pdf')), scale = 1.5, plot = p)

