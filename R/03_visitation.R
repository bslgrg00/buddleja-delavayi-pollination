source("R/01_data.R", local = TRUE)

# Summary ----
set.seed(SEED)
vis_summary <- vis |>
  group_by(year, visitor_group) |>
  summarise(n_plants = n(), plants_with_zero = sum(visits == 0),
            pooled_rate = 100 * sum(visits) / sum(flowers_observed), describe(rate), visits = sum(visits)) |>
  group_by(year) |>
  mutate(share_of_visits = visits / sum(visits)) |>
  ungroup()
save_table(vis_summary, "visitation_summary")

# Negative-binomial GLMM ----
d <- mutate(vis, year = factor(year))
m <- glmmTMB(visits ~ visitor_group * year + offset(log(flowers_observed)) + (1 | plant_id) + (1 | plant_year),
             family = nbinom2, data = d)
save_model(m, "visitation", d)

em <- emmeans(m, ~ visitor_group | year, offset = log(100))
ct <- as_tibble(summary(contrast(em, "pairwise"), infer = TRUE, type = "response"))
rates <- as_tibble(summary(em, type = "response")) |>
  transmute(year, visitor_group, rate = response, lwr = asymp.LCL, upr = asymp.UCL)
rates$letter <- unlist(lapply(levels(d$year), function(y) {
  pw <- filter(ct, year == y)
  pairs <- tibble(g1 = sub(" / .*$", "", pw$contrast), g2 = sub("^.* / ", "", pw$contrast), p = pw$p.value)
  r <- filter(rates, year == y)
  cld_letters(GROUPS, pairs, setNames(r$rate, as.character(r$visitor_group)))[GROUPS]
}))
save_table(rates, "visitation_rates")

between_years <- as_tibble(summary(contrast(emmeans(m, ~ year | visitor_group, offset = log(100)), "revpairwise"),
                                   infer = TRUE, type = "response"))
save_table(bind_rows(
  transmute(ct, comparison = as.character(year), contrast, rate_ratio = ratio, lwr = asymp.LCL, upr = asymp.UCL, p = p.value),
  transmute(between_years, comparison = paste(visitor_group, "2026 vs 2025"), contrast, rate_ratio = ratio,
            lwr = asymp.LCL, upr = asymp.UCL, p = p.value)), "visitation_contrasts")

# Exploratory: pollen limitation and bee visitation ----
pv <- plant_years |> filter(in_experiment) |> mutate(bee_rate = 100 * bees / flowers_observed, year = factor(year))
tests <- list(); preds <- list()
for (v in c("pli_fs", "pli_ss")) {
  m <- glmmTMB(as.formula(paste(v, "~ log(bee_rate) + year + (1 | plant_id)")), data = pv)
  save_model(m, paste0("visitation_", v), pv)
  tests[[v]] <- wald(m, "log(bee_rate)", paste(v, "vs log bee visitation rate"))
  preds[[v]] <- bind_rows(lapply(levels(pv$year), function(y) {
    r <- range(pv$bee_rate[pv$year == y])
    nd <- tibble(bee_rate = exp(seq(log(r[1]), log(r[2]), length.out = 60)), year = factor(y, levels(pv$year)), plant_id = NA)
    p <- predict(m, newdata = nd, re.form = NA, allow.new.levels = TRUE, se.fit = TRUE)
    mutate(nd, variable = v, estimate = p$fit, lwr = p$fit - qnorm(0.975) * p$se.fit, upr = p$fit + qnorm(0.975) * p$se.fit)
  }))
}
save_table(bind_rows(tests), "visitation_pli_tests")
save_table(bind_rows(preds) |> select(variable, year, bee_rate, estimate, lwr, upr), "visitation_pli_predictions")
