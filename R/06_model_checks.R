source("R/01_data.R", local = TRUE)

RESP <- c(`Fruit set` = "cbind(fruits, fruits_failed)", `Seed set` = "cbind(seeds, ovules_unfilled)")
RE <- "(1 | plant_id) + (1 | plant_year) + (1 | pair_id)"
u <- mutate(units, year = factor(year), index_c = index - mean(plant_years$index))
fit <- function(f, family = binomial, data = u) glmmTMB(as.formula(f), family = family, data = data)
treatment_or <- function(m) {
  e <- summary(contrast(emmeans(m, ~ treatment), list(`Supplemental / Control` = c(-1, 1))), infer = TRUE, type = "response")
  tibble(estimate = e$odds.ratio, lwr = e$asymp.LCL, upr = e$asymp.UCL, p = e$p.value)
}

# Simulated residuals for every saved model ----
residual_check <- function(m, data) {
  if (inherits(m, "glm"))
    return(list(stats = tibble(family = m$family$family, n = nobs(m), converged = m$converged,
                               dispersion_ratio = sum(residuals(m, "pearson")^2) / m$df.residual)))
  set.seed(SEED)
  yy <- model.response(m$frame); y <- if (is.matrix(yy)) yy[, 1] else as.numeric(yy)
  sims <- simulate(m, nsim = N_SIM)
  if (is.matrix(sims[[1]]) || is.data.frame(sims[[1]])) {
    S <- sapply(sims, function(s) as.matrix(s)[, 1]); size <- rowSums(as.matrix(sims[[1]]))
  } else { S <- as.matrix(as.data.frame(sims)); size <- 1 }
  mu <- predict(m, newdata = data, re.form = NA, type = "response") * size
  obs <- sum((y - mu)^2); sim <- colSums((S - mu)^2)
  pit <- vapply(seq_along(y), function(i) mean(S[i, ] < y[i]) + runif(1) * mean(S[i, ] == y[i]), 0)
  ks_py <- NA_real_
  if ("plant_year" %in% names(m$frame)) {        # units of a plant-year are not independent
    g <- m$frame$plant_year; Y <- tapply(y, g, sum); SS <- apply(S, 2, function(s) tapply(s, g, sum))
    pg <- vapply(seq_along(Y), function(i) mean(SS[i, ] < Y[i]) + runif(1) * mean(SS[i, ] == Y[i]), 0)
    ks_py <- suppressWarnings(ks.test(pg, "punif")$p.value)
  }
  list(stats = tibble(family = family(m)$family, n = length(y), converged = isTRUE(m$sdr$pdHess) && m$fit$convergence == 0,
                      dispersion_ratio = obs / mean(sim), p_dispersion = min(1, 2 * min(mean(sim >= obs), mean(sim <= obs))),
                      ks_p_units = suppressWarnings(ks.test(pit, "punif")$p.value), ks_p_plant_years = ks_py),
       pit = pit, mu = mu)
}
files <- list.files(DIR_MOD, "\\.rds$", full.names = TRUE)
checks <- lapply(files, function(f) { x <- readRDS(f); residual_check(x$model, x$data) })
names(checks) <- sub("\\.rds$", "", basename(files))
save_table(bind_rows(lapply(checks, `[[`, "stats"), .id = "model"), "residual_checks", DIR_CHK)

pdf(file.path(DIR_CHK, "residual_plots.pdf"), width = 7, height = 3.6)
for (k in names(checks)) {
  x <- checks[[k]]; if (is.null(x$pit)) next
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), cex = 0.8)
  plot(ppoints(length(x$pit)), sort(x$pit), pch = 16, cex = 0.5, xlab = "Expected", ylab = "Observed", main = paste(k, "(QQ)"))
  abline(0, 1, col = "grey50")
  plot(rank(x$mu) / length(x$mu), x$pit, pch = 16, cex = 0.5, ylim = c(0, 1),
       xlab = "Model prediction (rank)", ylab = "Scaled residual", main = "Residuals vs predicted")
  abline(h = c(0.25, 0.5, 0.75), lty = 2, col = "grey50")
}
invisible(dev.off())

# Binomial vs beta-binomial error distribution ----
errdist <- list()
for (r in names(RESP)) {
  mb  <- fit(paste(RESP[[r]], "~ treatment * year +", RE))
  mbb <- fit(paste(RESP[[r]], "~ treatment * year +", RE), betabinomial)
  mb0 <- fit(paste(RESP[[r]], "~ treatment + year +", RE), betabinomial)
  errdist[[r]] <- tibble(model = paste("Pollen limitation,", tolower(r)), aic_binomial = AIC(mb), aic_betabinomial = AIC(mbb),
                         effect = "treatment odds ratio (beta-binomial)", treatment_or(mbb),
                         p_treatment_x_year = anova(mb0, mbb)[["Pr(>Chisq)"]][2])
}
f_dist <- paste(RESP[["Seed set"]], "~ treatment * index_c + treatment * year +", RE)
mb <- fit(f_dist); mbb <- fit(f_dist, betabinomial)
errdist$dist <- tibble(model = "Disturbance, seed set", aic_binomial = AIC(mb), aic_betabinomial = AIC(mbb),
                       effect = "treatment x disturbance, ratio of odds ratios (beta-binomial)",
                       select(wald(mbb, "treatmentSupplemental:index_c", "", exp = TRUE), -test))
save_table(bind_rows(errdist) |> mutate(delta_aic = aic_betabinomial - aic_binomial, .after = aic_betabinomial),
           "error_distribution", DIR_CHK)

# Unit size ----
save_table(u |> group_by(year, treatment) |>
             summarise(n_units = n(), median_flowers = median(flowers), q1_flowers = quantile(flowers, 0.25),
                       q3_flowers = quantile(flowers, 0.75), median_ovules = median(ovules)),
           "unit_size_medians", DIR_CHK)
save_table(bind_rows(lapply(names(RESP), function(r) {
  m0 <- fit(paste(RESP[[r]], "~ treatment * year +", RE))
  m1 <- fit(paste(RESP[[r]], "~ treatment * year + log2(flowers) +", RE))
  bind_rows(treatment_or(m0) |> mutate(term = "treatment odds ratio, without unit size"),
            treatment_or(m1) |> mutate(term = "treatment odds ratio, with unit size"),
            wald(m1, "log2(flowers)", "", exp = TRUE) |> select(-test) |> mutate(term = "odds ratio per doubling of flowers")) |>
    mutate(response = r, aic = c(AIC(m0), AIC(m1), NA), .before = 1)
})) |> relocate(term, .after = response), "unit_size_models", DIR_CHK)

# Field scores of observer 1 instead of photograph scores of observer 2 ----
pyf <- plant_years |> filter(in_experiment) |> mutate(year = factor(year))
pv  <- mutate(plant_years, year = factor(year))
field <- bind_rows(lapply(c(observer2 = "index", observer1 = "index_field"), function(ix) {
  uu <- mutate(u, D = .data[[ix]]); pp <- mutate(pyf, D = .data[[ix]]); vv <- mutate(pv, D = .data[[ix]])
  bind_rows(
    wald(fit("bees ~ D + year + offset(log(flowers_observed)) + (1 | plant_id)", nbinom2, vv), "D",
         "Bee visits vs disturbance", exp = TRUE),
    bind_rows(lapply(names(RESP), function(r) wald(fit(paste(RESP[[r]], "~ treatment * D + treatment * year +", RE), binomial, uu),
                                                   "treatmentSupplemental:D", paste(r, "- treatment x disturbance"), exp = TRUE))),
    wald(fit("pli_fs ~ D + year + (1 | plant_id)", gaussian, pp), "D", "pli_fs vs disturbance"),
    wald(fit("pli_ss ~ D + year + (1 | plant_id)", gaussian, pp), "D", "pli_ss vs disturbance"))
}), .id = "scores")
save_table(arrange(field, test, scores), "field_scores", DIR_CHK)

# Robustness of the disturbance results ----
fs_int <- function(dd) wald(fit(paste(RESP[["Fruit set"]], "~ treatment * D + treatment * year +", RE), binomial,
                                mutate(dd, D = index - mean(index))), "treatmentSupplemental:D", "", exp = TRUE)
pli_sl <- function(dd) wald(fit("pli_fs ~ index + year + (1 | plant_id)", gaussian, dd), "index", "")
loo <- function(d, f) bind_rows(lapply(setNames(unique(d$plant_year), unique(d$plant_year)), function(p) f(filter(d, plant_year != p))), .id = "dropped")
hi <- max(pyf$index) - 1
li <- loo(u, fs_int); lp <- loo(pyf, pli_sl)
ri <- fs_int(filter(u, index < hi)); rp <- pli_sl(filter(pyf, index < hi))
rob <- bind_rows(
  tibble(check = "Fruit set treatment x disturbance: leave one plant-year out", estimate_min = min(li$estimate),
         estimate_max = max(li$estimate), p_max = max(li$p), detail = paste("largest P without", li$dropped[which.max(li$p)])),
  tibble(check = paste0("Fruit set treatment x disturbance: plant-years with index < ", hi, " only"), estimate_min = ri$estimate,
         estimate_max = ri$estimate, p_max = ri$p, detail = paste(n_distinct(filter(u, index < hi)$plant_year), "plant-years")),
  tibble(check = "PLI fruit set slope: leave one plant-year out", estimate_min = min(lp$estimate),
         estimate_max = max(lp$estimate), p_max = max(lp$p), detail = paste("largest P without", lp$dropped[which.max(lp$p)])),
  tibble(check = paste0("PLI fruit set slope: plant-years with index < ", hi, " only"), estimate_min = rp$estimate,
         estimate_max = rp$estimate, p_max = rp$p, detail = paste(sum(pyf$index < hi), "plant-years")),
  bind_rows(lapply(levels(pyf$year), function(y) {
    co <- summary(lm(pli_fs ~ index, data = filter(pyf, year == y)))$coefficients["index", ]
    tibble(check = paste("PLI fruit set slope within", y), estimate_min = co[[1]], estimate_max = co[[1]], p_max = co[[4]],
           detail = paste(sum(pyf$year == y), "plant-years"))
  })))
save_table(rob, "robustness", DIR_CHK)
