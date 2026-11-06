source("R/01_data.R", local = TRUE)

VARS <- c(fs_c = "Fruit set, control", fs_s = "Fruit set, supplemental",
          ss_c = "Seed set, control", ss_s = "Seed set, supplemental",
          pli_fs = "PLI, fruit set", pli_ss = "PLI, seed set",
          lrr_fs = "ln(S/C), fruit set", lrr_ss = "ln(S/C), seed set")
summarise_vars <- function(d, period)
  bind_rows(lapply(names(VARS), function(v) describe(d[[v]]) |> mutate(period = period, variable = VARS[[v]], .before = 1)))

# Plant means over the years sampled ----
py <- filter(plant_years, in_experiment)
plants <- py |>
  group_by(plant_id) |>
  summarise(n_years = n(), across(c(fs_c, fs_s, ss_c, ss_s), mean)) |>
  mutate(pli_fs = 1 - fs_c / fs_s, pli_ss = 1 - ss_c / ss_s, lrr_fs = log(fs_s / fs_c), lrr_ss = log(ss_s / ss_c))

set.seed(SEED); both <- summarise_vars(plants, "Both years")
set.seed(SEED); by_year <- bind_rows(lapply(c(2025, 2026), function(y) summarise_vars(filter(py, year == y), as.character(y))))
save_table(bind_rows(by_year, both), "pollen_limitation_summary")

paired <- bind_rows(lapply(list(c("fs_s", "fs_c", "Fruit set"), c("ss_s", "ss_c", "Seed set")), function(v) {
  s <- plants[[v[1]]]; c <- plants[[v[2]]]
  w <- suppressWarnings(wilcox.test(s, c, paired = TRUE, exact = TRUE)); tt <- t.test(s, c, paired = TRUE)
  tibble(response = v[3], n_plants = length(s), mean_control = mean(c), mean_supplemental = mean(s),
         mean_difference = mean(s - c), diff_lwr = tt$conf.int[1], diff_upr = tt$conf.int[2],
         n_supplemental_higher = sum(s > c), wilcoxon_v = unname(w$statistic), wilcoxon_p = w$p.value,
         t = unname(tt$statistic), t_p = tt$p.value)
}))
save_table(paired, "pollen_limitation_paired_tests")

# Binomial GLMMs ----
RESP <- c(`Fruit set` = "cbind(fruits, fruits_failed)", `Seed set` = "cbind(seeds, ovules_unfilled)")
RE <- "(1 | plant_id) + (1 | plant_year) + (1 | pair_id)"
u <- mutate(units, year = factor(year))
glmm <- list(); within <- list(); year_tests <- list()
for (r in names(RESP)) {
  m  <- glmmTMB(as.formula(paste(RESP[[r]], "~ treatment * year +", RE)), family = binomial, data = u)
  m0 <- glmmTMB(as.formula(paste(RESP[[r]], "~ treatment + year +", RE)), family = binomial, data = u)
  save_model(m, paste0("pollen_limitation_", sub(" ", "_", tolower(r))), u)

  em <- emmeans(m, ~ treatment)
  pr <- summary(em, type = "response")
  rt <- summary(contrast(regrid(em, transform = "log"), "pairwise"), infer = TRUE, type = "response")
  glmm[[r]] <- tibble(response = r, control = pr$prob[1], supplemental = pr$prob[2],
                      pli = 1 - rt$ratio, pli_lwr = 1 - rt$asymp.UCL, pli_upr = 1 - rt$asymp.LCL, p = rt$p.value)

  emy <- emmeans(m, ~ treatment | year)
  ct <- summary(contrast(emy, list(`Supplemental / Control` = c(-1, 1))), by = NULL, adjust = "sidak", infer = TRUE, type = "response")
  py_means <- as_tibble(summary(emy, type = "response")) |> select(year, treatment, prob) |>
    pivot_wider(names_from = treatment, values_from = prob, names_prefix = "model_") |> rename_with(tolower)
  within[[r]] <- as_tibble(ct) |>
    transmute(response = r, year, odds_ratio = odds.ratio, lwr = asymp.LCL, upr = asymp.UCL, z = z.ratio, p_sidak = p.value) |>
    left_join(py_means, by = "year")

  rty <- summary(contrast(regrid(emy, transform = "log"), "pairwise"), infer = TRUE, type = "response")
  year_tests[[r]] <- bind_rows(
    lrt(m0, m, paste0(r, ": treatment x year")) |> mutate(bh_family = TRUE),
    tibble(test = paste0(r, ": model PLI ", rty$year), estimate = 1 - rty$ratio,
           lwr = 1 - rty$asymp.UCL, upr = 1 - rty$asymp.LCL, bh_family = FALSE))
}
save_table(bind_rows(glmm), "pollen_limitation_glmm")
save_table(bind_rows(within), "pollen_limitation_within_year")

# Differences between years ----
pyf <- mutate(py, year = factor(year))
for (v in c("pli_fs", "pli_ss")) {
  m <- glmmTMB(as.formula(paste(v, "~ year + (1 | plant_id)")), data = pyf)
  save_model(m, paste0("years_", v), pyf)
  year_tests[[v]] <- wald(m, "year2026", paste0(VARS[[v]], ": 2026 - 2025")) |> mutate(bh_family = TRUE)
}
save_table(bind_rows(year_tests) |> relocate(bh_family, .after = last_col()), "year_tests")

# Plants sampled in both years ----
rep_long <- py |>
  filter(plant_id %in% plant_id[duplicated(plant_id)]) |>
  select(plant_id, year, all_of(names(VARS)[1:6])) |>
  pivot_longer(-c(plant_id, year), names_to = "variable") |>
  pivot_wider(names_from = year, values_from = value, names_prefix = "value_")
repeated <- rep_long |>
  group_by(variable) |>
  summarise(n_plants = n(), mean_2025 = mean(value_2025), mean_2026 = mean(value_2026),
            mean_change = mean(value_2026 - value_2025),
            n_increase = sum(value_2026 > value_2025), n_decrease = sum(value_2026 < value_2025),
            wilcoxon_p = suppressWarnings(wilcox.test(value_2026, value_2025, paired = TRUE, exact = TRUE)$p.value),
            spearman_rho = cor(value_2025, value_2026, method = "spearman")) |>
  mutate(label = VARS[variable], .after = variable) |>
  arrange(match(variable, names(VARS)))
save_table(repeated, "repeated_plants")
