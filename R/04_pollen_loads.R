source("R/01_data.R", local = TRUE)

# Summary ----
pol_summary <- pollen |>
  group_by(visitor_group) |>
  summarise(n_insects = n(), insects_2025 = sum(year == 2025), insects_2026 = sum(year == 2026),
            on_study_plants = sum(!is.na(plant_id)), n_with_pollen = sum(total > 0),
            n_carriers = sum(carries_bd), pct_carriers = 100 * mean(carries_bd),
            pct_loads_mostly_bd = 100 * mean(bd_share[total > 0] >= 0.5),
            total_median = median(total), total_q1 = quantile(total, 0.25), total_q3 = quantile(total, 0.75), total_mean = mean(total),
            bd_median = median(bd), bd_q1 = quantile(bd, 0.25), bd_q3 = quantile(bd, 0.75), bd_mean = mean(bd),
            bd_share_median = median(bd_share, na.rm = TRUE), bd_grains = sum(bd)) |>
  mutate(share_of_bd_grains = bd_grains / sum(bd_grains))
save_table(pol_summary, "pollen_load_summary")

# Models ----
estimates <- function(m, response, back = identity) {
  e <- as_tibble(summary(emmeans(m, ~ visitor_group), type = "response"))
  lc <- grep("LCL|lower", names(e), value = TRUE)[1]; uc <- grep("UCL|upper", names(e), value = TRUE)[1]
  tibble(response = response, visitor_group = e$visitor_group, estimate = back(e[[2]]), lwr = back(e[[lc]]), upr = back(e[[uc]]))
}
with_bd  <- filter(pollen, bd > 0) |> mutate(y = log10(bd))
with_pol <- filter(pollen, total > 0) |> mutate(y = log10(total))

m_carrier <- glmmTMB(carries_bd ~ visitor_group, family = binomial, data = pollen)
m_bd      <- glmmTMB(y ~ visitor_group, data = with_bd)
m_total   <- glmmTMB(y ~ visitor_group, data = with_pol)
m_share   <- glm(bd_share ~ visitor_group, family = quasibinomial, data = with_pol)
save_model(m_carrier, "pollen_carrier", pollen); save_model(m_bd, "pollen_bd_load", with_bd)
save_model(m_total, "pollen_total_load", with_pol); save_model(m_share, "pollen_bd_share", with_pol)

f_share <- anova(glm(bd_share ~ 1, family = quasibinomial, data = with_pol), m_share, test = "F")
pol_tests <- bind_rows(
  lrt(glmmTMB(carries_bd ~ 1, family = binomial, data = pollen), m_carrier, "Pollen loads: carrying B. delavayi pollen") |> mutate(n = nrow(pollen)),
  lrt(glmmTMB(y ~ 1, data = with_bd), m_bd, "Pollen loads: B. delavayi grains per carrier") |> mutate(n = nrow(with_bd)),
  lrt(glmmTMB(y ~ 1, data = with_pol), m_total, "Pollen loads: grains per insect with pollen") |> mutate(n = nrow(with_pol)),
  tibble(test = "Pollen loads: B. delavayi share of the load", F = f_share$F[2], df = f_share$Df[2],
         p = f_share$`Pr(>F)`[2], n = nrow(with_pol))) |>
  mutate(bh_family = TRUE)
save_table(pol_tests, "pollen_load_tests")
save_table(bind_rows(
  estimates(m_carrier, "Carrying B. delavayi pollen"),
  estimates(m_bd, "B. delavayi grains per carrier (geometric mean)", function(v) 10^v),
  estimates(m_total, "Grains per insect with pollen (geometric mean)", function(v) 10^v),
  estimates(m_share, "B. delavayi share of the load")), "pollen_load_estimates")

# Pollen-transport index ----
rates_wide <- vis |> select(plant_year, visitor_group, rate) |> pivot_wider(names_from = visitor_group, values_from = rate)
loads <- split(pollen$bd, pollen$visitor_group)
transport <- function(v, l) { x <- sapply(GROUPS, function(g) mean(v[[g]]) * mean(l[[g]])); x / sum(x) }
set.seed(SEED)
B <- replicate(N_BOOT %/% 2, {
  v <- rates_wide[sample(nrow(rates_wide), replace = TRUE), ]
  l <- lapply(loads, function(x) x[sample.int(length(x), replace = TRUE)])
  transport(v, l)
})
save_table(tibble(visitor_group = GROUPS, mean_rate = sapply(GROUPS, function(g) mean(rates_wide[[g]])),
                  mean_bd_load = sapply(loads, mean), n_insects = lengths(loads),
                  share = transport(rates_wide, loads),
                  share_lwr = apply(B, 1, quantile, 0.025), share_upr = apply(B, 1, quantile, 0.975)),
           "pollen_transport_index")

# Composition of the loads ----
long <- pollen |>
  select(insect_no, visitor_group, all_of(POLLEN_TYPES)) |>
  pivot_longer(all_of(POLLEN_TYPES), names_to = "pollen_type", values_to = "grains")
save_table(long |> group_by(pollen_type) |> summarise(grains = sum(grains)) |>
             mutate(share = grains / sum(grains), label = pollen_label(pollen_type)) |> arrange(desc(grains)),
           "pollen_type_totals")

W <- as.matrix(pollen[, POLLEN_TYPES]); keep <- rowSums(W) > 0
P <- W[keep, ] / rowSums(W[keep, ]); g <- droplevels(pollen$visitor_group[keep])
comp <- as_tibble(P) |>
  mutate(visitor_group = g) |>
  pivot_longer(-visitor_group, names_to = "pollen_type", values_to = "prop") |>
  group_by(visitor_group, pollen_type) |>
  summarise(mean_pct = 100 * mean(prop), pct_insects_carrying = 100 * mean(prop > 0)) |>
  left_join(count(tibble(visitor_group = g), visitor_group, name = "n_insects_with_pollen"), by = "visitor_group") |>
  mutate(label = pollen_label(pollen_type)) |>
  arrange(visitor_group, desc(mean_pct))
save_table(comp, "pollen_composition")

# PERMANOVA and PERMDISP on Bray-Curtis dissimilarities
D <- as.matrix(dist(P, method = "manhattan")) / 2
permanova <- function(D, g, nperm = N_PERM) {
  D2 <- D^2; n <- length(g); k <- nlevels(g); sst <- sum(D2[upper.tri(D2)]) / n
  ssw <- function(gg) sum(vapply(split(seq_len(n), gg), function(ix)
    if (length(ix) > 1) sum(D2[ix, ix][upper.tri(D2[ix, ix])]) / length(ix) else 0, 0))
  fstat <- function(gg) { s <- ssw(gg); ((sst - s) / (k - 1)) / (s / (n - k)) }
  f0 <- fstat(g); fp <- replicate(nperm, fstat(sample(g)))
  tibble(pseudo_F = f0, R2 = (sst - ssw(g)) / sst, p = (sum(fp >= f0) + 1) / (nperm + 1), n = n)
}
permdisp <- function(D, g, nperm = N_PERM) {
  pc <- suppressWarnings(cmdscale(D, k = nrow(D) - 1, eig = TRUE))   # negative eigenvalues are dropped
  X <- pc$points[, pc$eig[seq_len(ncol(pc$points))] > 1e-10, drop = FALSE]
  z <- vapply(seq_len(nrow(X)), function(i) sqrt(sum((X[i, ] - colMeans(X[g == g[i], , drop = FALSE]))^2)), 0)
  fstat <- function(gg) anova(lm(z ~ gg))$`F value`[1]
  f0 <- fstat(g); fp <- replicate(nperm, fstat(sample(g)))
  tibble(pseudo_F = f0, p = (sum(fp >= f0) + 1) / (nperm + 1), n = length(g),
         note = paste(levels(g), round(tapply(z, g, mean), 3), sep = " mean distance to centroid ", collapse = "; "))
}
set.seed(SEED)
pm <- permanova(D, g) |> mutate(test = "PERMANOVA: composition ~ visitor group", .before = 1)
dp <- permdisp(D, g) |> mutate(test = "PERMDISP: dispersion ~ visitor group", .before = 1)
pw <- bind_rows(lapply(combn(levels(g), 2, simplify = FALSE), function(pr) {
  ix <- g %in% pr
  permanova(D[ix, ix], droplevels(g[ix])) |> mutate(test = paste("Pairwise:", pr[1], "vs", pr[2]), .before = 1)
})) |> mutate(p_holm = p.adjust(p, "holm"))
save_table(bind_rows(pm, pw, dp), "pollen_composition_tests")
