fast_asca <- function(formula, data, type = 3) {
  # --- 1. Setup & Formula Parsing ---
  clean_form <- mixlm::rparse(formula, REML = FALSE)
  orig_terms <- attr(terms(formula), "term.labels")
  is_random  <- grepl("r\\(", orig_terms)

  old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old_opts))

  resp_var <- all.vars(clean_form)[1]
  Y <- data[[resp_var]]; if (!is.matrix(Y)) Y <- as.matrix(Y)
  n <- nrow(Y); p <- ncol(Y)
  # Denominator for Expl.var(%)
  Y_centered <- sweep(Y, 2, colMeans(Y), "-")
  SS_total_global <- sum(Y_centered^2)

  rhs_form <- as.formula(paste("~", as.character(clean_form)[3]))
  X <- model.matrix(rhs_form, data = data)
  assign <- attr(X, "assign")
  term_labels <- attr(terms(rhs_form), "term.labels")
  n_terms <- length(term_labels)

  # --- 2. SS Calculations (Type I, II, or III) ---
  SS_matrix <- matrix(0, n_terms, p, dimnames = list(term_labels, NULL))
  Effect_Matrices <- list()
  QR_full <- qr(X)
  Y_hat_full <- qr.fitted(QR_full, Y)
  SSE_full <- colSums((Y - Y_hat_full)^2)

  for (i in 1:n_terms) {
    if (type == 1) {
      Y_hat_curr <- qr.fitted(qr(X[, which(assign <= i), drop=FALSE]), Y)
      Y_hat_prev <- qr.fitted(qr(X[, which(assign < i), drop=FALSE]), Y)
    } else if (type == 2) {
      is_higher <- grepl(paste0("^", term_labels[i], ":"), term_labels) |
        grepl(paste0(":", term_labels[i], "$"), term_labels)
      Y_hat_curr <- qr.fitted(qr(X[, which(!(assign %in% which(is_higher))), drop=FALSE]), Y)
      Y_hat_prev <- qr.fitted(qr(X[, which(!(assign %in% c(i, which(is_higher)))), drop=FALSE]), Y)
    } else { # Type III
      Y_hat_curr <- Y_hat_full
      Y_hat_prev <- qr.fitted(qr(X[, which(assign != i), drop=FALSE]), Y)
    }
    Effect_Matrices[[term_labels[i]]] <- Y_hat_curr - Y_hat_prev
    SS_matrix[i, ] <- colSums(Effect_Matrices[[term_labels[i]]]^2)
  }

  # --- 3. Degrees of Freedom & MS ---
  df_res <- n - QR_full$rank
  df_terms <- as.vector(table(assign)[-1])
  names(df_terms) <- term_labels
  MS_matrix <- sweep(SS_matrix, 1, df_terms, "/")
  MS_res <- SSE_full / df_res

  # --- 4. MoM Denominator Mapping & P-values ---
  F_matrix <- matrix(NA, n_terms, p)
  P_matrix <- matrix(NA, n_terms, p)
  MS_error_matrix <- matrix(NA, n_terms, p) # The missing "projected error"
  den_labels <- character(n_terms)
  den_dfs <- numeric(n_terms)

  for (i in 1:n_terms) {
    candidates <- which(is_random & grepl(term_labels[i], term_labels) & (1:n_terms != i))
    if (length(candidates) > 0) {
      target <- candidates[which.max(attr(terms(rhs_form), "order")[candidates])]
      MS_error_matrix[i, ] <- MS_matrix[target, ]
      den_dfs[i] <- df_terms[target]
      den_labels[i] <- term_labels[target]
    } else {
      MS_error_matrix[i, ] <- MS_res
      den_dfs[i] <- df_res
      den_labels[i] <- "Residuals"
    }
    F_matrix[i, ] <- MS_matrix[i, ] / MS_error_matrix[i, ]
    P_matrix[i, ] <- pf(F_matrix[i, ], df_terms[i], den_dfs[i], lower.tail = FALSE)
  }

  # --- 5. Global Summary Table ---
  total_ss_vec <- c(rowSums(SS_matrix), sum(SSE_full))
  summary_tab <- data.frame(
    Df = c(df_terms, df_res),
    Sum.Sq = total_ss_vec,
    `Expl.var.(%)` = (total_ss_vec / SS_total_global) * 100,
    Error.Term = c(den_labels, "None"),
    stringsAsFactors = FALSE
  )
  colnames(summary_tab) <- c("Df", "Sum.Sq.", "Expl.var.(%)", "Error.Term")
  rownames(summary_tab) <- c(term_labels, "Residuals")

  return(list(Summary = summary_tab, Effects = Effect_Matrices,
              P_values = P_matrix, MS_errors = MS_error_matrix,
              df = list(terms = df_terms, res = df_res, den = den_dfs),
              Y_residuals = Y - Y_hat_full))
}

candies2 <- candies[-c(1,7,23),]

# 9 variabler
SStype <- 3
reps <- 10
cat("------------\n")
T <- system.time(for(i in 1:reps)mod <- hdanova(assessment ~ candy * assessor, data=candies2, SStype = SStype))
summary(mod)
Tq <- system.time(for(i in 1:reps)modq <- fast_asca(assessment ~ candy * assessor, candies2, type = SStype))
cat(T[3]/Tq[3], "x raskere\n")
modq$Summary
Tm <- system.time(for(i in 1:reps)modm <- hdanova(assessment ~ candy * r(assessor), data=candies2, SStype = SStype))
summary(modm)
Tqm <- system.time(for(i in 1:reps)modqm <- fast_asca(assessment ~ candy * r(assessor), candies2, type = SStype))
cat(Tm[3]/Tqm[3], "x raskere\n")
modqm$Summary


# 1632 variabler
SStype <- 3
reps <- 5
cat("------------\n")
T <- system.time(for(i in 1:reps)mod <- hdanova(transcriptome ~ strain*growthrate, data=Lactobacillus, SStype = SStype))
summary(mod)
Tq <- system.time(for(i in 1:reps)modq <- fast_asca(transcriptome ~ strain*growthrate, data=Lactobacillus, type = SStype))
cat(T[3]/Tq[3], "x raskere", paste0("(",round(T[3],4),",",round(Tq[3],4),")\n"))
modq$Summary
Tm <- system.time(for(i in 1:reps)modm <- hdanova(transcriptome ~ strain*r(growthrate), data=Lactobacillus, SStype = SStype))
summary(modm)
Tqm <- system.time(for(i in 1:reps)modqm <- fast_asca(transcriptome ~ strain*r(growthrate), data=Lactobacillus, type = SStype))
cat(Tm[3]/Tqm[3], "x raskere", paste0("(",round(Tm[3],4),",",round(Tqm[3],4),")\n"))
modqm$Summary


# -----------------------------------------------------------
scoreplot_asca(modqm, "candy", candies2, "candy")
# Corrected scoreplot function with accurate axes and dropped residuals
scoreplot_asca <- function(res, term = "candy", data, factor_name = "candy", comps = c(1, 2)) {
  # 1. Base PCA on the Effect Matrix
  E <- res$Effects[[term]]
  svd_res <- svd(E)

  # Scores of the Level Means (Centers) and Loadings
  level_scores <- svd_res$u %*% diag(svd_res$d)
  loadings <- svd_res$v

  # 2. Project ONLY the Mixed Model Interaction (No Residuals)
  err_label <- res$Summary[term, "Error.Term"]
  if (err_label == "Residuals" || is.na(err_label)) {
    Error_Matrix <- matrix(0, nrow(E), ncol(E))
  } else {
    Error_Matrix <- res$Effects[[err_label]]
  }

  # Full Observation Scores: (Means + Interaction) %*% Loadings
  obs_scores <- (E + Error_Matrix) %*% loadings

  # 3. Standard ASCA Explained Variance
  # Based on eigenvalues (d^2) of the Effect Matrix only
  eig_values <- svd_res$d^2
  perc_expl <- (eig_values / sum(eig_values)) * 100

  # 4. Plotting
  plot(obs_scores[, comps],
       xlab = paste0("PC", comps[1], " (", round(perc_expl[comps[1]], 2), "%)"),
       ylab = paste0("PC", comps[2], " (", round(perc_expl[comps[2]], 2), "%)"),
       main = paste("ASCA Score Plot: ", term),
       pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3), cex = 0.8)

  # Overlay Level Centers using Factor Levels
  fac <- as.factor(data[[factor_name]])
  lvls <- levels(fac)
  for(i in seq_along(lvls)) {
    idx <- which(fac == lvls[i])
    # Plot the center of the observations for this level
    points(mean(obs_scores[idx, comps[1]]), mean(obs_scores[idx, comps[2]]),
           pch = 21, bg = i+1, col = "black", cex = 1.8)
    text(mean(obs_scores[idx, comps[1]]), mean(obs_scores[idx, comps[2]]),
         labels = lvls[i], pos = 3, font = 2)
  }

  abline(h = 0, v = 0, lty = 2, col = "gray")
}
