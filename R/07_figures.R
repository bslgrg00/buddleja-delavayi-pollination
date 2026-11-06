source("R/01_data.R", local = TRUE)

PT <- 1 / .pt                                  # converts font size in pt to ggplot text size
jit <- function(w = 0.1) position_jitter(width = w, height = 0, seed = SEED)
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s.")))
py <- filter(plant_years, in_experiment) |> mutate(year = factor(year))
pts_year <- list(scale_shape_manual(values = c(`2025` = 21, `2026` = 21)),
                 scale_fill_manual(values = c(`2025` = "white", `2026` = "#5E5E5E")))
RATE_LAB <- "Visits per 100 flowers"

# Figure 1: visitor rings, 2025 (to be placed on the map) ----
R_MAX <- 9; R_MIN <- 1.9; CELL <- 30; NCOL <- 6     # ring radius (mm) for the highest rate and the smallest ring; grid cell (mm)
ring <- function(x, y, r, shares, id, inner = 0.55, n = 90) {
  shares <- shares / sum(shares); a0 <- pi / 2; out <- list()
  for (k in seq_along(GROUPS)) {
    if (shares[k] <= 0) next
    th <- seq(a0, a0 - 2 * pi * shares[k], length.out = max(4, ceiling(n * shares[k])))
    out[[k]] <- tibble(id = paste(id, k), group = factor(GROUPS[k], GROUPS),
                       x = c(x + r * cos(th), rev(x + inner * r * cos(th))), y = c(y + r * sin(th), rev(y + inner * r * sin(th))))
    a0 <- a0 - 2 * pi * shares[k]
  }
  bind_rows(out)
}
circle <- function(x, y, r, n = 100) tibble(x = x + r * cos(seq(0, 2 * pi, length.out = n)), y = y + r * sin(seq(0, 2 * pi, length.out = n)))

f1 <- plant_years |>
  filter(year == 2025) |>
  arrange(as.integer(gsub("\\D", "", field_label))) |>
  mutate(total = 100 * (bees + butterflies + others) / flowers_observed, k = row_number() - 1)
K <- R_MAX / sqrt(max(f1$total))
H1 <- CELL * ceiling(nrow(f1) / NCOL) + 34
f1 <- mutate(f1, r = pmax(K * sqrt(total), R_MIN), x = CELL / 2 + CELL * (k %% NCOL), y = H1 - CELL / 2 - 2 - CELL * (k %/% NCOL))
rings <- bind_rows(lapply(seq_len(nrow(f1)), function(i) with(f1[i, ], ring(x, y, r, c(bees, butterflies, others), plant_id))))
lost <- bind_rows(lapply(which(!f1$in_experiment), function(i) circle(f1$x[i], f1$y[i], f1$r[i] + 1.5) |> mutate(id = f1$plant_id[i])))
key_v <- c(10, 40, 100); kx <- 22; ky <- 6
key <- bind_rows(lapply(key_v, function(v) circle(kx, ky + K * sqrt(v), K * sqrt(v)) |> mutate(v = v)))
leg <- tibble(group = factor(GROUPS, GROUPS), x = 70 + 30 * (0:2), y = 14)
fig1 <- ggplot() +
  geom_polygon(data = rings, aes(x, y, group = id, fill = group), colour = "white", linewidth = 0.25) +
  geom_path(data = lost, aes(x, y, group = id), colour = INK2, linewidth = 0.35, linetype = "22") +
  geom_text(data = filter(f1, r >= 5), aes(x, y, label = round(total)), size = 6 * PT, colour = INK) +
  geom_text(data = f1, aes(x, y - R_MAX - 3, label = plant_id), size = 7 * PT, colour = INK, fontface = "bold") +
  geom_text(data = f1, aes(x, y - R_MAX - 6, label = field_label), size = 6 * PT, colour = MUTED) +
  geom_path(data = key, aes(x, y, group = v), colour = MUTED, linewidth = 0.3) +
  geom_segment(data = tibble(v = key_v), aes(x = kx, xend = kx + K * sqrt(100) + 3, y = ky + 2 * K * sqrt(v), yend = ky + 2 * K * sqrt(v)),
               colour = MUTED, linewidth = 0.2) +
  geom_text(data = tibble(v = key_v), aes(x = kx + K * sqrt(100) + 4, y = ky + 2 * K * sqrt(v), label = v), size = 6 * PT, hjust = 0, colour = INK2) +
  annotate("text", x = kx, y = ky - 3, label = RATE_LAB, size = 6.5 * PT, colour = INK2) +
  geom_tile(data = leg, aes(x, y, fill = group), width = 3.2, height = 3.2) +
  geom_text(data = leg, aes(x + 3, y, label = group), size = 7 * PT, hjust = 0, colour = INK) +
  { if (nrow(lost)) list(geom_path(data = circle(70, 5, 1.6), aes(x, y), colour = INK2, linewidth = 0.35, linetype = "22"),
                         annotate("text", x = 73, y = 5, label = "Experiment lost", size = 7 * PT, hjust = 0, colour = INK)) } +
  scale_fill_manual(values = COL_GROUP, guide = "none") +
  coord_equal(xlim = c(0, 180), ylim = c(0, H1), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA))
save_fig(fig1, "Fig1_visitor_rings_2025", 180, H1)

# Figure 2: pollen limitation ----
OFF <- 0.75 / 4
brackets <- function(b) list(
  geom_segment(data = b, aes(x = x - OFF, xend = x + OFF, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3, colour = INK2),
  geom_segment(data = b, aes(x = x - OFF, xend = x - OFF, y = y, yend = y - h), inherit.aes = FALSE, linewidth = 0.3, colour = INK2),
  geom_segment(data = b, aes(x = x + OFF, xend = x + OFF, y = y, yend = y - h), inherit.aes = FALSE, linewidth = 0.3, colour = INK2),
  geom_text(data = b, aes(x = x, y = y + 0.5 * h, label = lab, size = ifelse(lab == "n.s.", 2.5, 3.2)), inherit.aes = FALSE,
            vjust = 0, colour = INK), scale_size_identity())
dodged_box <- function(d, x, fill, b, values) ggplot(d, aes({{ x }}, value, fill = {{ fill }})) +
  geom_boxplot(width = 0.62, outlier.shape = NA, colour = LINE, linewidth = 0.35, position = position_dodge(0.75)) +
  geom_point(aes(colour = {{ fill }}), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75, seed = SEED),
             size = 1.1, shape = 21, fill = "white", stroke = 0.3, show.legend = FALSE) +
  brackets(b) + scale_fill_manual(values = values) +
  scale_colour_manual(values = setNames(rep(INK2, length(values)), names(values)), guide = "none") + theme_pub()
within <- load_table("pollen_limitation_within_year")
box_trt <- function(resp, cols, tag) {
  d <- py |> select(year, C = all_of(cols[1]), S = all_of(cols[2])) |> pivot_longer(c(C, S)) |>
    mutate(treatment = factor(if_else(name == "C", "Control", "Supplemental"), names(COL_TRT)))
  rng <- diff(range(d$value))
  b <- d |> group_by(year) |> summarise(top = max(value)) |>
    left_join(filter(within, response == resp) |> mutate(year = factor(year)), by = "year") |>
    transmute(x = as.integer(year), y = top + 0.07 * rng, h = 0.03 * rng, lab = stars(p_sidak))
  dodged_box(d, year, treatment, b, COL_TRT) +
    scale_y_continuous(labels = percent_format(1), expand = expansion(mult = c(0.03, 0.12))) + labs(x = NULL, y = resp, tag = tag)
}
yt <- load_table("year_tests")
pli <- py |> select(year, `Fruit set` = pli_fs, `Seed set` = pli_ss) |> pivot_longer(-year, names_to = "response") |>
  mutate(response = factor(response))
rngC <- diff(range(c(0, pli$value)))
bC <- pli |> group_by(response) |> summarise(top = max(value)) |>
  mutate(p = yt$p[match(paste0("PLI, ", tolower(response), ": 2026 - 2025"), yt$test)]) |>
  transmute(x = as.integer(response), y = top + 0.07 * rngC, h = 0.03 * rngC, lab = stars(p))
p2C <- dodged_box(pli, response, year, bC, c(`2025` = "#EDEDED", `2026` = "#A6A6A6")) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "22", colour = MUTED) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) + labs(x = NULL, y = "Pollen limitation index", tag = "C")
fig2 <- (box_trt("Fruit set", c("fs_c", "fs_s"), "A") + box_trt("Seed set", c("ss_c", "ss_s"), "B") +
           plot_layout(guides = "collect") & theme(legend.position = "top")) | p2C
save_fig(fig2 + plot_layout(widths = c(2, 1)), "Fig2_pollen_limitation", 180, 88)

# Figure 3: plants sampled in both years ----
rp <- py |> filter(plant_id %in% plant_id[duplicated(plant_id)]) |> select(plant_id, year, pli_fs, pli_ss) |>
  pivot_longer(c(pli_fs, pli_ss), names_to = "variable", values_to = "pli") |>
  mutate(variable = c(pli_fs = "Fruit set", pli_ss = "Seed set")[variable])
rp <- bind_rows(rp, rp |> group_by(year, variable) |> summarise(pli = mean(pli)) |> mutate(plant_id = "Mean")) |>
  mutate(plant_id = factor(plant_id, rev(c(sort(unique(rp$plant_id)), "Mean"))))
seg <- pivot_wider(rp, names_from = year, values_from = pli, names_prefix = "y")
fig3 <- ggplot(rp, aes(y = plant_id)) +
  geom_segment(data = seg, aes(x = y2025, xend = y2026, yend = plant_id), colour = "#BDBDBD", linewidth = 0.8) +
  geom_point(aes(x = pli, shape = year, fill = year), size = 2.3, colour = LINE, stroke = 0.4) +
  geom_hline(yintercept = 1.5, colour = "#D9D9D9", linewidth = 0.3) +
  facet_wrap(~ variable, nrow = 1) + pts_year +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Pollen limitation index (PLI)", y = NULL) + theme_pub() +
  theme(strip.text = element_text(face = "plain", hjust = 0), panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.3))
save_fig(fig3, "Fig3_repeated_plants", 180, 72)

# Figure 4: visitation and pollen transport ----
vr <- load_table("visitation_rates") |> mutate(visitor_group = factor(visitor_group, GROUPS))
vtop <- vis |> group_by(year, visitor_group) |> summarise(ymax = max(rate)) |> left_join(vr, by = c("year", "visitor_group"))
p4A <- ggplot(vis, aes(visitor_group, rate, fill = visitor_group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, colour = LINE, linewidth = 0.35, alpha = 0.85) +
  geom_point(position = jit(0.12), size = 1.1, shape = 21, fill = "white", colour = INK2, stroke = 0.3) +
  geom_text(data = vtop, aes(y = ymax * 1.6 + 0.3, label = letter), size = 2.9) +
  facet_wrap(~ year) + scale_fill_manual(values = COL_GROUP, guide = "none") +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.5, base = 10), breaks = c(0, 1, 10, 100)) +
  labs(x = NULL, y = RATE_LAB, tag = "A") + theme_pub()

totals <- load_table("pollen_type_totals")
shown <- totals$pollen_type[totals$share >= 0.003 | totals$pollen_type %in% c("b_delavayi", "unidentified")]
net <- pollen |> select(visitor_group, all_of(POLLEN_TYPES)) |>
  pivot_longer(-visitor_group, names_to = "pollen_type", values_to = "grains") |>
  mutate(pollen_type = if_else(pollen_type %in% shown, pollen_type, "other_plants")) |>
  group_by(visitor_group, pollen_type) |> summarise(grains = sum(grains))
types <- net |> group_by(pollen_type) |> summarise(G = sum(grains)) |>
  arrange(desc(pollen_type == "b_delavayi"), desc(if_else(pollen_type == "other_plants", -Inf, G))) |>
  mutate(y = rev(seq_len(n())), label = pollen_label(pollen_type), italic = !pollen_type %in% c("unidentified", "other_plants"),
         G_lab = ifelse(G >= 1e6, paste0(formatC(G / 1e6, format = "f", digits = 2), "M"), paste0(formatC(G / 1e3, format = "f", digits = 1), "k")),
         expr = paste0(ifelse(italic, paste0("italic('", label, "')"), paste0("'", label, "'")), "~~'(", G_lab, ")'"))
grp <- net |> group_by(visitor_group) |> summarise(G = sum(grains)) |> left_join(count(pollen, visitor_group, name = "n"), by = "visitor_group") |>
  mutate(y = rescale(rev(seq_along(GROUPS))[match(visitor_group, GROUPS)], to = range(types$y) + c(0.6, -0.6)))
links <- net |> filter(grains > 0) |> left_join(select(grp, visitor_group, y1 = y), by = "visitor_group") |>
  left_join(select(types, pollen_type, y2 = y), by = "pollen_type") |> mutate(bd = pollen_type == "b_delavayi")
p4B <- ggplot() +
  geom_segment(data = links, aes(x = 0, xend = 1, y = y1, yend = y2, linewidth = grains, colour = bd), alpha = 0.55, lineend = "round") +
  geom_point(data = grp, aes(0, y, size = G, fill = visitor_group), shape = 21, colour = LINE, stroke = 0.4) +
  geom_point(data = types, aes(1, y, size = G, fill = pollen_type == "b_delavayi"), shape = 21, colour = LINE, stroke = 0.3) +
  geom_text(data = grp, aes(-0.08, y, label = paste0(visitor_group, "\n(n = ", n, ")")), hjust = 1, size = 2.6, lineheight = 0.9) +
  geom_text(data = types, aes(1.07, y, label = expr), parse = TRUE, hjust = 0, size = 2.4) +
  scale_linewidth(range = c(0.15, 4), guide = "none") + scale_size_area(max_size = 9, guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "#4A6378", `FALSE` = "#C8C8C8"), breaks = c("TRUE", "FALSE"),
                      labels = c(expression(italic("B. delavayi") ~ "pollen"), "Other pollen"), name = NULL) +
  scale_fill_manual(values = c(COL_GROUP, `TRUE` = "#4A6378", `FALSE` = "#E3E3E3"), guide = "none") +
  scale_x_continuous(limits = c(-0.55, 2.0)) + labs(tag = "B") + theme_void(base_family = "sans") +
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 11))
save_fig(p4A / p4B + plot_layout(heights = c(1, 1.25)), "Fig4_visitors_pollen", 180, 181)

# Figure 5 and S2: disturbance ----
pr <- load_table("disturbance_predictions")
idx_x <- scale_x_continuous(breaks = 0:max(plant_years$index))
bp <- filter(pr, response == "Bee visitation") |> mutate(year = factor(group))
p5A <- ggplot(mutate(plant_years, rate = 100 * bees / flowers_observed, year = factor(year)), aes(index, rate)) +
  geom_ribbon(data = bp, aes(y = estimate, ymin = lwr, ymax = upr, group = year), fill = "#E6E6E6", alpha = 0.7) +
  geom_line(data = bp, aes(y = estimate, linetype = year), colour = LINE, linewidth = 0.5) +
  geom_point(aes(shape = year, fill = year), position = jit(), size = 1.6, colour = LINE, stroke = 0.35) +
  geom_text(data = bp |> group_by(year) |> slice_max(index, n = 1), aes(x = index + 0.12, y = estimate, label = year),
            hjust = 0, size = 2.5, colour = INK2) +
  scale_y_log10() + pts_year + scale_linetype_manual(values = c(`2025` = "dashed", `2026` = "solid"), guide = "none") +
  scale_x_continuous(breaks = 0:max(plant_years$index), expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = "Disturbance index", y = paste("Bee", tolower(RATE_LAB)), tag = "A") + theme_pub()
trt_dist <- function(resp, cols, tag) {
  pts <- py |> select(index, C = all_of(cols[1]), S = all_of(cols[2])) |> pivot_longer(c(C, S), names_to = "treatment") |>
    mutate(treatment = factor(if_else(treatment == "C", "Control", "Supplemental"), names(COL_TRT)))
  p <- filter(pr, response == resp) |> mutate(treatment = factor(group, names(COL_TRT)))
  ggplot(pts, aes(index, value)) +
    geom_ribbon(data = p, aes(y = estimate, ymin = lwr, ymax = upr, group = treatment), fill = "#E6E6E6", alpha = 0.7) +
    geom_point(aes(fill = treatment), shape = 21, position = jit(), size = 1.6, colour = LINE, stroke = 0.35) +
    geom_line(data = p, aes(y = estimate, colour = treatment, linetype = treatment), linewidth = 0.6) +
    geom_text(data = p |> group_by(treatment) |> slice_max(index, n = 1), aes(x = index + 0.12, y = estimate, label = treatment, colour = treatment),
              hjust = 0, size = 2.5) +
    scale_fill_manual(values = COL_TRT, guide = "none") + scale_colour_manual(values = LINE_TRT, guide = "none") +
    scale_linetype_manual(values = c(Control = "dashed", Supplemental = "solid"), guide = "none") +
    scale_y_continuous(labels = percent_format(1)) +
    scale_x_continuous(breaks = 0:max(plant_years$index), expand = expansion(mult = c(0.04, 0.30))) +
    labs(x = "Disturbance index", y = resp, tag = tag) + theme_pub()
}
pli_dist <- function(v, resp, tag) {
  p <- filter(pr, response == v)
  ggplot(py, aes(index, .data[[v]])) +
    geom_ribbon(data = p, aes(y = estimate, ymin = lwr, ymax = upr), fill = "#E6E6E6") +
    geom_line(data = p, aes(y = estimate), colour = LINE, linewidth = 0.5) +
    geom_point(aes(shape = year, fill = year), position = jit(), size = 1.6, colour = LINE, stroke = 0.35) +
    pts_year + idx_x + labs(x = "Disturbance index", y = paste("PLI,", tolower(resp)), tag = tag) + theme_pub()
}
loss <- load_table("loss_by_cutting") |> mutate(total = retained + lost, pct = lost / total)
p5D <- ggplot(loss, aes(factor(cutting), pct)) +
  geom_col(width = 0.55, fill = "#BDBDBD", colour = LINE, linewidth = 0.3) +
  geom_text(aes(label = paste0(lost, "/", total)), vjust = -0.5, size = 2.6) +
  scale_y_continuous(labels = percent_format(1), limits = c(0, max(loss$pct) * 1.35 + 0.05), expand = expansion(mult = c(0, 0))) +
  labs(x = "Cutting score", y = "Plant-years with experiment lost", tag = "D") + theme_pub()
fig5 <- (p5A + trt_dist("Fruit set", c("fs_c", "fs_s"), "B")) / (pli_dist("pli_fs", "Fruit set", "C") + p5D) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(fig5, "Fig5_disturbance", 180, 155)
save_fig(trt_dist("Seed set", c("ss_c", "ss_s"), "A") + pli_dist("pli_ss", "Seed set", "B") +
           plot_layout(guides = "collect") & theme(legend.position = "bottom"), "FigS2_disturbance_seed_set", 180, 83)

# Figure S1: sampling design ----
s1 <- plant_years |> mutate(plant_id = factor(plant_id, rev(sort(unique(plant_id)))), year = factor(year),
                            label = paste0(field_label, " (", coalesce(n_pairs, 0L), ")"))
pS1 <- ggplot(s1, aes(year, plant_id)) +
  geom_line(data = filter(s1, repeated), aes(group = plant_id), colour = "#CFCFCF") +
  geom_point(aes(fill = index, shape = in_experiment), size = 2.8, colour = LINE, stroke = 0.3) +
  geom_text(aes(label = label, hjust = ifelse(year == "2025", 1.25, -0.25)), size = 2.1, colour = INK2) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), labels = c(`TRUE` = "experiment", `FALSE` = "visitation only"), name = NULL) +
  scale_y_discrete(limits = levels(s1$plant_id)) +
  scale_fill_gradient(low = "#F2F2F2", high = "#4D4D4D", name = "Disturbance index") +
  scale_x_discrete(expand = expansion(add = 0.9)) + labs(x = NULL, y = NULL) + theme_pub(8) +
  theme(legend.title = element_text(size = 7, colour = INK2))
save_fig(pS1, "FigS1_sampling_design", 180, 150)

# Figure S3: pollen limitation and bee visitation ----
vp <- load_table("visitation_pli_predictions") |> mutate(year = factor(year))
pv <- mutate(py, bee_rate = 100 * bees / flowers_observed)
vis_panel <- function(v, resp, tag) ggplot(pv, aes(bee_rate, .data[[v]])) +
  geom_ribbon(data = filter(vp, variable == v), aes(y = estimate, ymin = lwr, ymax = upr, group = year), fill = "#E6E6E6", alpha = 0.7) +
  geom_line(data = filter(vp, variable == v), aes(y = estimate, linetype = year), colour = LINE, linewidth = 0.5) +
  geom_point(aes(shape = year, fill = year), size = 1.6, colour = LINE, stroke = 0.35) +
  scale_x_log10() + pts_year + scale_linetype_manual(values = c(`2025` = "dashed", `2026` = "solid"), guide = "none") +
  labs(x = paste("Bee", tolower(RATE_LAB)), y = paste("PLI,", tolower(resp)), tag = tag) + theme_pub()
save_fig(vis_panel("pli_fs", "Fruit set", "A") + vis_panel("pli_ss", "Seed set", "B") +
           plot_layout(guides = "collect") & theme(legend.position = "bottom"), "FigS3_visitation_vs_PLI", 180, 83)

# Figure S4: composition of pollen loads ----
cmp <- load_table("pollen_composition")
top <- cmp |> filter(!pollen_type %in% c("b_delavayi", "unidentified")) |> group_by(pollen_type) |>
  summarise(m = mean(mean_pct)) |> slice_max(m, n = 3) |> pull(pollen_type)
lev <- c("b_delavayi", top, "unidentified", "other_plants")
c5 <- cmp |> mutate(cat = if_else(pollen_type %in% lev, pollen_type, "other_plants")) |>
  group_by(visitor_group, n_insects_with_pollen, cat) |> summarise(pct = sum(mean_pct)) |> ungroup() |>
  mutate(cat = factor(cat, rev(lev)), ylab = paste0(visitor_group, "\n(n = ", n_insects_with_pollen, ")"),
         ylab = factor(ylab, rev(unique(ylab[order(match(visitor_group, GROUPS))]))))
labs5 <- setNames(ifelse(lev %in% c("unidentified", "other_plants"), paste0("'", pollen_label(lev), "'"),
                         paste0("italic('", pollen_label(lev), "')")), lev)
pS4 <- ggplot(c5, aes(pct, ylab, fill = cat)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  geom_text(data = filter(c5, pct >= 7), aes(label = paste0(round(pct), "%"), colour = cat == "b_delavayi"),
            position = position_stack(vjust = 0.5), size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_manual(values = setNames(c("#4A6378", "#8C8C8C", "#ADADAD", "#CBCBCB", "#E2E2E2", "#F3F3F3"), lev),
                    breaks = lev, labels = parse(text = labs5[lev]), name = NULL) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Mean share of an insect's pollen load", y = NULL) + theme_pub(8) +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7, hjust = 0)) +
  guides(fill = guide_legend(nrow = 2))
save_fig(pS4, "FigS4_pollen_load_composition", 180, 72)
