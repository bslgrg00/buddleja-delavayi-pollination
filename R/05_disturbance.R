source("R/01_data.R", local = TRUE)

RESP <- c(`Fruit set` = "cbind(fruits, fruits_failed)", `Seed set` = "cbind(seeds, ovules_unfilled)")
RE <- "(1 | plant_id) + (1 | plant_year) + (1 | pair_id)"
centre <- mean(plant_years$index)
grid <- round(seq(min(plant_years$index), max(plant_years$index), by = 0.1), 1)

# Agreement between observers ----
wkappa <- function(a, b, k = 3) {
  O <- table(factor(a, 0:(k - 1)), factor(b, 0:(k - 1))) / length(a); E <- outer(rowSums(O), colSums(O))
  W <- 1 - abs(outer(0:(k - 1), 0:(k - 1), "-")) / (k - 1); pe <- sum(W * E)
  if (pe >= 1) NA_real_ else (sum(W * O) - pe) / (1 - pe)
}
set.seed(SEED)
agree <- bind_rows(lapply(c("herbivory", "cutting"), function(comp) {
  a <- disturbance[[paste0(comp, "_observer1")]]; b <- disturbance[[paste0(comp, "_observer2")]]
  boot <- replicate(2000, { i <- sample(length(a), replace = TRUE); wkappa(a[i], b[i]) })
  tibble(component = comp, n_plant_years = length(a), exact_agreement = mean(a == b), within_one = mean(abs(a - b) <= 1),
         weighted_kappa = wkappa(a, b), kappa_lwr = quantile(boot, 0.025, na.rm = TRUE), kappa_upr = quantile(boot, 0.975, na.rm = TRUE),
         observer1_higher = sum(a > b), observer2_higher = sum(a < b))
}))
rho <- cor.test(plant_years$index_field, plant_years$index, method = "spearman", exact = FALSE)
agree <- bind_rows(agree, tibble(component = "index", n_plant_years = nrow(plant_years),
                                 exact_agreement = mean(plant_years$index_field == plant_years$index),
                                 within_one = mean(abs(plant_years$index_field - plant_years$index) <= 1),
                                 spearman_rho = unname(rho$estimate), spearman_p = rho$p.value))
save_table(agree, "interobserver_agreement")

pairs_c <- list(c("road", "herbivory"), c("herbivory", "cutting"), c("road", "cutting"))
save_table(bind_rows(lapply(pairs_c, function(pr) {
  ct <- cor.test(plant_years[[pr[1]]], plant_years[[pr[2]]], method = "spearman", exact = FALSE)
  tibble(pair = paste(pr, collapse = " - "), rho = unname(ct$estimate), p = ct$p.value, n = nrow(plant_years))
})), "disturbance_components")

tests <- list(); preds <- list()

# Visitation ----
pv <- plant_years |> mutate(index_c = index - centre, year = factor(year), total = bees + butterflies + others)
for (resp in c("bees", "total")) {
  m <- glmmTMB(as.formula(paste(resp, "~ index_c + year + offset(log(flowers_observed)) + (1 | plant_id)")), family = nbinom2, data = pv)
  save_model(m, paste0("disturbance_", resp, "_visits"), pv)
  tests[[resp]] <- wald(m, "index_c", paste(if (resp == "bees") "Bee" else "Total", "visits vs disturbance"), exp = TRUE) |>
    mutate(effect = "rate ratio per index point", bh_family = TRUE)
  if (resp == "bees")
    preds$bees <- as_tibble(summary(emmeans(m, ~ index_c | year, at = list(index_c = grid - centre), offset = log(100)), type = "response")) |>
      transmute(estimate = response, lwr = asymp.LCL, upr = asymp.UCL, response = "Bee visitation",
                group = as.character(year), index = index_c + centre)
}

# Reproduction: treatment x disturbance ----
u <- mutate(units, year = factor(year), index_c = index - centre)
slopes <- function(m, var, label, bh) {
  tr <- summary(emtrends(m, ~ treatment, var = var), infer = TRUE)
  bind_rows(
    tibble(test = paste(label, "- slope in", tr$treatment), estimate = exp(tr[[2]]), lwr = exp(tr$asymp.LCL),
           upr = exp(tr$asymp.UCL), p = tr$p.value, effect = "odds ratio per index point", bh_family = FALSE),
    wald(m, paste0("treatmentSupplemental:", var), paste(label, "- treatment x disturbance"), exp = TRUE) |>
      mutate(effect = "ratio of odds ratios", bh_family = bh))
}
for (r in names(RESP)) {
  m <- glmmTMB(as.formula(paste(RESP[[r]], "~ treatment * index_c + treatment * year +", RE)), family = binomial, data = u)
  save_model(m, paste0("disturbance_", sub(" ", "_", tolower(r))), u)
  tests[[r]] <- slopes(m, "index_c", r, bh = TRUE)
  preds[[r]] <- as_tibble(summary(emmeans(m, ~ treatment | index_c, at = list(index_c = grid - centre)), type = "response")) |>
    transmute(response = r, group = as.character(treatment), index = index_c + centre, estimate = prob, lwr = asymp.LCL, upr = asymp.UCL)
  for (v in c("road", "herbivory", "cutting", "damage")) {
    ms <- glmmTMB(as.formula(paste(RESP[[r]], "~ treatment *", v, "+ treatment * year +", RE)), family = binomial, data = u)
    tests[[paste(r, v)]] <- slopes(ms, v, paste0(r, " [", v, "]"), bh = FALSE)
  }
}

# Pollen limitation per plant-year ----
pyf <- plant_years |> filter(in_experiment) |> mutate(year = factor(year))
for (v in c("pli_fs", "pli_ss")) {
  m <- glmmTMB(as.formula(paste(v, "~ index + year + (1 | plant_id)")), data = pyf)
  save_model(m, paste0("disturbance_", v), pyf)
  tests[[v]] <- wald(m, "index", paste(v, "vs disturbance")) |> mutate(effect = "change per index point", bh_family = FALSE)
  e <- as_tibble(summary(emmeans(m, ~ index, at = list(index = grid))))
  preds[[v]] <- tibble(response = v, group = NA_character_, index = e$index, estimate = e$emmean,
                       lwr = e[[grep("LCL|lower", names(e))[1]]], upr = e[[grep("UCL|upper", names(e))[1]]])
}

# Loss of experiments ----
loss <- mutate(plant_years, lost = !in_experiment, cut2 = cutting == 2)
ft <- fisher.test(table(loss$cut2, loss$lost))
tests$loss <- tibble(test = "Experiment lost vs cutting score 2 (Fisher)", estimate = unname(ft$estimate),
                     lwr = ft$conf.int[1], upr = ft$conf.int[2], p = ft$p.value, effect = "odds ratio", bh_family = TRUE)
save_table(loss |> count(cutting, lost) |> mutate(lost = if_else(lost, "lost", "retained")) |>
             pivot_wider(names_from = lost, values_from = n, values_fill = 0), "loss_by_cutting")

save_table(bind_rows(tests) |> relocate(effect, .after = test), "disturbance_tests")
save_table(bind_rows(preds) |> transmute(response, group, index = round(index, 1), estimate, lwr, upr), "disturbance_predictions")
