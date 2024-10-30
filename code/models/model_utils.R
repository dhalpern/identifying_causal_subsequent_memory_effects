add_model_preds <- function(mod, assessment_df, colname = "glmnet_preds", s = "lambda.1se", offset = NULL) {
  # Fit the model to the 90%
  if(!is.null(offset)) {
    newoffset <- assessment_df[[offset]]
  } else {
    newoffset <- NULL
  }
  pred_df <- assessment_df %>% mutate(
    {{colname}} := predict(mod, assessment_df, newoffset = newoffset, type = "response", s = s)[, 1]
  ) %>% select(lith_word, sub_id, recall_acc, glmer_preds, {{colname}})
}

create_jol_metric_df <- function(clf_tbl) {
  jol_df <- clf_tbl %>% filter(var == "jol") %>% 
    ungroup %>% 
    select(sub_id, model_type, glmnet_metrics)
  neural_df <- clf_tbl %>% filter(var != "jol")
  new_clf_tbl <- neural_df %>% inner_join(jol_df, by = c("sub_id", "model_type"), suffix = c(".neural", ".jol"))
  return(new_clf_tbl)
}

mod_pred_t.test <- function(fit_mod_df, metric, col1 = 'glmnet_metrics', col2 = 'glmer_metrics', 
                            alternative = "two.sided", paired = TRUE, 
                            mu = 0, ret_zero_if_const = FALSE) {
  col1_metric <- sym(str_glue("{col1}.metric"))
  fit_mod_df_auc <- fit_mod_df %>% filter((!!col1_metric) == metric)
  col1 <- fit_mod_df_auc %>% pull(str_glue("{col1}.estimate"))
  if(!is.null(col2)) {
    col2 <- fit_mod_df_auc %>% pull(str_glue("{col2}.estimate"))
  }
  res <- try(t.test(col1, col2, alternative = alternative,
                    paired = paired, mu = mu))
  if(!(typeof(res) == "list")) {
    if(ret_zero_if_const) {
      tibble(statistic = 0, estimate = mu, constant_error = TRUE)
    } else {
      return(NULL)
    }
  } else {
    return(tidy(res))
  }
}

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

fit_cvglmnet_ <- function(data, formula, offset, family, intercept, alpha, 
                          use_foldid = TRUE, standardize = TRUE) {
  if(use_foldid == TRUE) {
    foldid <- data$foldid
  } else {
    foldid <- NULL
  }
  if(length(all.vars(formula)) == 2) {
    formula <- update(formula, ~ . + fake_intercept)
  }
  
  if(is.null(offset)) {
    cv.glmnet(formula, data = data, family = family, 
              intercept = intercept, alpha = alpha, foldid = foldid, standardize = standardize)
  } else {
    cv.glmnet(formula, data = data, offset = data[[offset]], family = family, 
              intercept = intercept, alpha = alpha, foldid = foldid, standardize = standardize)
  }
}

add_glmer_offset <- function(df, causal_m, standard_m, offset_type, re.form = NULL) {
    ret_df <- df %>% mutate(
      causal_offset = predict(causal_m, newdata = df, re.form = re.form, allow.new.levels = TRUE),
      standard_offset = predict(standard_m, newdata = df, re.form = re.form, allow.new.levels = TRUE),
      glmer_preds = predict(causal_m, newdata = df, re.form = re.form, allow.new.levels = TRUE, type = "response")
    )
  return(ret_df)
}

add_glm_offset <- function(df, m) {
  # print(model.matrix(m, data = df))
  df %>% mutate(
      offset = predict(m, newdata = df),
      glm_preds = predict(m, newdata = df, type = "response")
    ) 
}

add_folds_col <- function(df, group_col_str, k) {
  #first permute 1:k, then create n_groups assignments then permute assignments
  groups <- df %>% pull(!!quo_name(group_col_str)) %>% unique()
  g_folds <- sample(rep_len(sample(k), length(groups)))
  fold_tbl <- tibble(!!quo_name(group_col_str) := groups,
                     foldid = g_folds)
  df %>% left_join(fold_tbl, by = group_col_str)
}

fit_cvglmnet <- function(formula, model_type, rs_obj, alpha = 1, 
                         use_offset = TRUE, 
                         s = "lambda.1se", standardize = TRUE) {
  if (use_offset) {
    intercept <- FALSE
    offset_col <- paste(model_type, 'offset', sep = '_')
  } else {
    intercept <- TRUE
    offset_col <- NULL
  }
  print(formula)
  print(offset_col)
  
  rs_obj %>% mutate(
    cvglmnet_fit = map(analysis_df, fit_cvglmnet_, formula, offset = offset_col, family = "binomial", 
                       intercept = intercept, alpha = alpha, standardize = standardize),
    pred_df = map2(cvglmnet_fit, assessment_df, add_model_preds, offset = offset_col, s = s)
  ) %>% select(-analysis_df, -assessment_df, -splits, -id, -causal_glmer_fit, -standard_glmer_fit)
}

check_conv <- function(mod) {
  if(!("code" %in% names(mod@optinfo$conv$lme4))) {
    return(TRUE)
  } else if((mod@optinfo$conv$opt == 0) & (mod@optinfo$conv$lme4$code == 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

fit_glmer <- function(formula, data, family = binomial, max_n_tries = 5) {
  print(paste('n_tries', 0))
  gm <- glmer(formula, family = family, data = data)
  conv_opt <- check_conv(gm)
  n_tries <- 1
  print(paste('converged?', conv_opt))
  print(n_tries)
  
  while(!(check_conv(gm)) & n_tries < max_n_tries) {
    ss <- getME(gm, c("theta","fixef"))
    gm <- update(gm, start=ss, control=glmerControl(optCtrl=list(maxfun=((1e4) * n_tries))))
    n_tries <- n_tries + 1
    print(n_tries)
    print(check_conv(gm))
  }
  print('end')
  print(check_conv(gm))
  return(gm)
}
