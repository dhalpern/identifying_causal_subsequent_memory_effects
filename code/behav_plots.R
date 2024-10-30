library(tidyverse)
library(lme4)
library(cowplot)
library(rstanarm)
options(mc.cores = parallel::detectCores())

base_dir <- '/Users/david/Documents/GitHub/identifying_causal_subsequent_memory_effects'
fig_dir <- file.path(base_dir, 'doc', 'figs')

trial_word_df <- read_csv(file.path(base_dir, 'data', 'mri_subs_events_wide.csv'))
trial_word_df <- trial_word_df %>% 
  dplyr::select(lith_word, sub_id, recall_acc) %>% unique()
 
memorability_df <- trial_word_df %>% group_by(lith_word) %>% 
  summarise(memorability = mean(recall_acc))
sub_ability_df <- trial_word_df %>% group_by(sub_id) %>% 
  summarise(sub_ability = mean(recall_acc))

memorability_df %>% summarise(mean(memorability), median(memorability), sd(memorability))
sub_ability_df %>% summarise(mean(sub_ability), median(sub_ability), sd(sub_ability))

m_plot <- memorability_df %>% 
  ggplot(aes(x = memorability)) + 
  geom_histogram(binwidth = .1, color = "darkgoldenrod1", fill = "darkgoldenrod1", alpha = .4) + 
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  xlab('Word Pair Prop. Correct') +
  ylab('Number of Pairs') + 
  scale_y_continuous(breaks = c(5, 10), expand = c(0, 0)) + 
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8), limits = c(0, 1), expand = c(0, 0))
s_plot <- sub_ability_df %>% 
  ggplot(aes(x = sub_ability)) + 
  geom_histogram(binwidth = .1, color = "darkred", fill = "darkred", alpha = .4) + 
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  xlab("Subject Prop. Correct") +
  ylab('Number of Subjects') +
  scale_y_continuous(breaks = c(5, 10), expand = c(0, 0)) + 
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8), limits = c(0, 1), expand = c(0, 0))

save_plot(file.path(fig_dir, 'm_plot.pdf'), m_plot)
save_plot(file.path(fig_dir, 's_plot.pdf'), s_plot)
pg <- plot_grid(m_plot, s_plot, labels = "auto", label_size = 20)
save_plot(file.path(fig_dir, 'behav_plots.pdf'), pg, ncol = 2)

pg <- plot_grid(m_plot, s_plot, labels = c("B", "C"), label_size = 20)
plot_grid(p3, bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)


# online vs.  lab comparison ----------------------------------------------


trial_word_df <- trial_word_df %>% mutate(
  sub_id_num = str_extract(sub_id, "[0-9]+"),
  test_location_cap = if_else(sub_id_num >= 62, "Online", "Lab"),
  test_location = if_else(sub_id_num >= 62, "online", "lab")
)

memorability_df <- trial_word_df %>% group_by(lith_word, test_location) %>% 
  summarise(memorability = mean(recall_acc), 
            n_subs = n(),
            n_correct = sum(recall_acc)
            ) %>% mutate(
              se = sqrt(memorability * (1 - memorability) * (1 / n_subs)),
              ci = 1.96 * se,
              memorability_ac = (n_correct + 2) / (n_subs + 4),
              se_ac = sqrt(memorability_ac * (1 - memorability_ac) * (1 / (n_subs + 4))),
              ci_ac = 1.96 * se_ac
)

sub_ability_df <- trial_word_df %>% group_by(sub_id, test_location_cap) %>% 
  summarise(sub_ability = mean(recall_acc))

s_plot <- sub_ability_df %>% 
  ggplot(aes(x = sub_ability, fill = test_location_cap)) + 
  geom_histogram(binwidth = .1) + 
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  labs(fill = "Test Location") +
  xlab("Subject Prop. Correct") +
  ylab('Number of Subjects') +
  scale_y_continuous(breaks = c(5, 10), expand = c(0, 0)) + 
  scale_x_continuous(breaks = c(0, .2, .4, .6, .8), limits = c(0, 1), expand = c(0, 0))
s_plot

memorability_plot_df <- memorability_df %>% pivot_wider(id_cols = c(lith_word), 
                                names_from = test_location, 
                                values_from = c(memorability, ci)) %>% mutate(
                                  online_lo = memorability_online - ci_online,
                                  online_hi = memorability_online + ci_online,
                                  lab_lo = memorability_lab - ci_lab,
                                  lab_hi = memorability_lab + ci_lab
                                )

memorability_ac_plot_df <- memorability_df %>% pivot_wider(id_cols = c(lith_word), 
                                names_from = test_location, 
                                values_from = c(memorability, n_subs, memorability_ac, ci_ac)) %>% mutate(
                                  online_ac_lo = memorability_ac_online - ci_ac_online,
                                  online_ac_hi = memorability_ac_online + ci_ac_online,
                                  lab_ac_lo = memorability_ac_lab - ci_ac_lab,
                                  lab_ac_hi = memorability_ac_lab + ci_ac_lab,
                                  # Truncate following Gelman, Hill and Vehtari (2020, pg. 52)
                                  online_ac_hi = if_else(online_ac_hi > 1, 1, online_ac_hi),
                                  online_ac_lo = if_else(online_ac_lo < 0, 0, online_ac_lo),
                                  lab_ac_hi = if_else(lab_ac_hi > 1, 1, lab_ac_hi),
                                  lab_ac_lo = if_else(lab_ac_lo < 0, 0, lab_ac_lo)
                                )

m_plot <- ggplot(memorability_ac_plot_df, aes(x = memorability_ac_lab, 
                                              y = memorability_ac_online)) +
  geom_abline() +
  geom_pointrange(aes(ymin = online_ac_lo, ymax = online_ac_hi)) +
  geom_pointrange(aes(xmin = lab_ac_lo, xmax = lab_ac_hi)) + 
  theme_classic() + 
  theme(text = element_text(size = 20)) +
  xlab("Prop. Correct (Lab)") +
  ylab('Prop. Correct (Online)') +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
m_plot

pg <- plot_grid(s_plot, m_plot, labels = "auto", label_size = 20)
pg
save_plot(file.path(fig_dir, 'behav_plots_online_lab.pdf'), pg, ncol = 2)

gm_recall <- glmer(recall_acc ~ test_location + (1 | sub_id) + (1 | lith_word), 
                   family = binomial, 
                   data = trial_word_df)
summary(gm_recall)

#singular
gm_recall <- glmer(recall_acc ~ test_location + (1 | sub_id) + (1 | lith_word/test_location), 
                   family = binomial, 
                   data = trial_word_df)
summary(gm_recall)

stan_gm_recall <- stan_glmer(recall_acc ~ test_location + (1 | sub_id) + (1 | lith_word/test_location), 
                   family = binomial, 
                   data = trial_word_df)
prior_summary(stan_gm_recall)
summary(stan_gm_recall, digits = 2, regex_pars = "igma|test_locationonline", probs = c(.025, .5, .975))
param_plot <- plot(stan_gm_recall, regex_pars = "igma|test_locationonline", prob_outer = 0.95) + 
  scale_y_discrete(
    labels = c("test_locationonline" = expression(beta),
               "Sigma[test_location:lith_word:(Intercept),(Intercept)]" = expression(sigma[gamma]),
               "Sigma[sub_id:(Intercept),(Intercept)]" = expression(sigma[theta]),
               "Sigma[lith_word:(Intercept),(Intercept)]" = expression(sigma[eta])
               )
)
ggsave(file.path(fig_dir, 'behav_params_rstanarm_online_lab.pdf'), param_plot)
