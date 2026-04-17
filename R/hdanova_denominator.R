.hda_collect_effect_structure <- function(formula, term_labels, model_frame, mixed_r,
                                          unrestricted, random_info = NULL, mom_anova = NULL){
  all_effects <- term_labels
  fixed_effects <- term_labels
  random_effects <- character(0)
  main_randoms_only_in_interactions <- character(0)

  if(mixed_r && is.list(random_info) && length(random_info) > 1){
    all_effects <- random_info$all
    fixed_effects <- random_info$fixed
    random_effects <- random_info$random
    main_randoms_only_in_interactions <- random_info$main.rands.only.inter
  }

  main_effects <- .hda_parse_main_effects(formula)
  n_levels <- numeric(length(main_effects))
  names(n_levels) <- main_effects
  for(i in seq_along(main_effects)){
    if(main_effects[i] %in% names(model_frame) && is.factor(model_frame[[main_effects[i]]]))
      n_levels[i] <- nlevels(model_frame[[main_effects[i]]])
  }

  list(
    formula = formula,
    term_labels = term_labels,
    all_effects = all_effects,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    main_randoms_only_in_interactions = main_randoms_only_in_interactions,
    main_effects = main_effects,
    n_levels = n_levels,
    mixed_r = mixed_r,
    unrestricted = unrestricted,
    mom_anova = mom_anova
  )
}


.hda_find_denominator_candidates <- function(structure){
  all_effects <- structure$all_effects
  n_effects <- length(all_effects)
  candidates <- vector("list", n_effects)
  names(candidates) <- all_effects

  for(i in seq_len(n_effects)){
    this_effect <- strsplit(all_effects[i], ":", fixed = TRUE)[[1]]
    which_contains <- integer(0)
    for(j in seq_len(n_effects)){
      effect_parts <- strsplit(all_effects[j], ":", fixed = TRUE)[[1]]
      if(i != j && all(this_effect %in% effect_parts) && length(effect_parts) > length(this_effect))
        which_contains <- union(which_contains, j)
    }
    candidates[[i]] <- sort(which_contains)
  }

  candidates
}


.hda_build_denominator_rules <- function(structure, candidates, mom_anova = NULL, df_res){
  term_labels <- structure$term_labels
  rules <- vector("list", length(term_labels))
  names(rules) <- term_labels

  if(!is.null(mom_anova) && !is.null(mom_anova$err.terms)){
    anova_tab <- mom_anova$anova
    if(is.list(anova_tab) && !is.data.frame(anova_tab))
      anova_tab <- anova_tab[[1]]
    err_terms <- mom_anova$err.terms
    if(is.null(names(err_terms)) && !is.null(anova_tab))
      names(err_terms) <- rownames(anova_tab)

    for(i in seq_along(term_labels)){
      term <- term_labels[i]
      den <- err_terms[[term]]
      if(is.null(den) && length(err_terms) >= i)
        den <- err_terms[[i]]
      if(is.null(den) || (length(den) == 1 && is.na(den))){
        rules[[i]] <- list(indices = length(term_labels) + 1,
                           coefficients = 1,
                           label = "Residuals",
                           df = df_res)
      } else {
        idx <- suppressWarnings(as.numeric(names(den)))
        labels <- ifelse(idx == length(term_labels) + 1, "Residuals", term_labels[idx])
        rules[[i]] <- list(indices = idx,
                           coefficients = as.numeric(den),
                           label = paste(unique(labels), collapse = "+"),
                           df = if(term %in% names(mom_anova$denom.df)) mom_anova$denom.df[term] else NA_real_)
      }
    }
    return(rules)
  }

  for(i in seq_along(term_labels)){
    rules[[i]] <- list(indices = length(term_labels) + 1,
                       coefficients = 1,
                       label = "Residuals",
                       df = df_res,
                       candidates = candidates[[i]])
  }
  rules
}


.hda_satterthwaite_df <- function(mean_squares, coefficients, dfs){
  numerator <- sum(mean_squares * coefficients)^2
  denominator <- sum(((mean_squares * coefficients)^2) / dfs)
  if(isTRUE(all.equal(denominator, 0)))
    return(NA_real_)
  numerator / denominator
}


.hda_apply_denominator_rules <- function(rules, MS_matrix, MS_res, df_terms, df_res, term_labels){
  n_terms <- length(term_labels)
  p <- ncol(MS_matrix)
  MS_error_matrix <- matrix(NA_real_, nrow = n_terms, ncol = p,
                            dimnames = list(term_labels, colnames(MS_matrix)))
  den_labels <- rep("Residuals", n_terms)
  den_dfs <- rep(df_res, n_terms)

  for(i in seq_along(term_labels)){
    rule <- rules[[i]]
    ms_term <- numeric(p)
    component_dfs <- numeric(length(rule$indices))
    component_ms <- numeric(length(rule$indices))

    for(j in seq_along(rule$indices)){
      idx <- rule$indices[j]
      coef_j <- rule$coefficients[j]
      if(idx == length(term_labels) + 1){
        ms_term <- ms_term + coef_j * MS_res
        component_dfs[j] <- df_res
        component_ms[j] <- mean(MS_res)
      } else if(idx >= 1 && idx <= length(term_labels)){
        ms_term <- ms_term + coef_j * MS_matrix[idx, ]
        component_dfs[j] <- df_terms[idx]
        component_ms[j] <- mean(MS_matrix[idx, ])
      }
    }

    if(all(ms_term == 0))
      ms_term <- MS_res
    MS_error_matrix[i, ] <- ms_term
    den_labels[i] <- rule$label
    den_dfs[i] <- if(!is.null(rule$df) && !is.na(rule$df)) {
      rule$df
    } else if(length(component_dfs) == 1) {
      component_dfs[1]
    } else {
      .hda_satterthwaite_df(component_ms, rule$coefficients, component_dfs)
    }
  }

  list(MS_error_matrix = MS_error_matrix,
       den_labels = den_labels,
       den_dfs = den_dfs)
}


.hda_parse_main_effects <- function(formula){
  rhs <- as.character(formula)[3]
  if(is.na(rhs) || rhs == "")
    return(character(0))
  terms_obj <- terms(stats::reformulate(rhs))
  labels <- attr(terms_obj, "term.labels")
  unique(unlist(lapply(labels, function(term) strsplit(term, ":", fixed = TRUE)[[1]])))
}