suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales); library(glmmTMB); library(emmeans)
})
options(dplyr.summarise.inform = FALSE, readr.show_col_types = FALSE)
emm_options(msg.interaction = FALSE)

# Settings ----
SEED   <- 20261102
N_BOOT <- 5000
N_SIM  <- 1000
N_PERM <- 9999
GROUPS <- c("Bees", "Butterflies", "Others")

DIR_DATA <- "data"
DIR_TAB  <- file.path("outputs", "tables")
DIR_MOD  <- file.path("outputs", "models")
DIR_CHK  <- file.path("outputs", "model_checks")
DIR_FIG  <- file.path("outputs", "figures")
for (d in c(DIR_TAB, DIR_MOD, DIR_CHK, DIR_FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Helpers ----
save_table <- function(x, name, dir = DIR_TAB) write_csv(x, file.path(dir, paste0(name, ".csv")), na = "")
save_model <- function(m, name, data) saveRDS(list(model = m, data = data), file.path(DIR_MOD, paste0(name, ".rds")))
load_table <- function(name, dir = DIR_TAB) read_csv(file.path(dir, paste0(name, ".csv")))

boot_ci <- function(x, fun = mean, n = N_BOOT) {
  b <- replicate(n, fun(sample(x, replace = TRUE)))
  quantile(b, c(0.025, 0.975), names = FALSE)
}
describe <- function(x) {
  x <- x[is.finite(x)]; ci <- boot_ci(x)
  tibble(n = length(x), mean = mean(x), sd = sd(x), median = median(x),
         q1 = quantile(x, 0.25, names = FALSE), q3 = quantile(x, 0.75, names = FALSE),
         min = min(x), max = max(x), boot_lwr = ci[1], boot_upr = ci[2])
}
lrt <- function(m0, m1, test) {
  chisq <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
  df <- attr(logLik(m1), "df") - attr(logLik(m0), "df")
  tibble(test = test, chisq = chisq, df = df, p = pchisq(chisq, df, lower.tail = FALSE))
}
wald <- function(m, term, test, exp = FALSE) {
  co <- summary(m)$coefficients$cond[term, ]
  ci <- co[[1]] + c(-1, 1) * qnorm(0.975) * co[[2]]
  tr <- if (exp) base::exp else identity
  tibble(test = test, estimate = tr(co[[1]]), lwr = tr(ci[1]), upr = tr(ci[2]), p = co[[4]])
}

pollen_label <- function(type) {
  u <- unique(type)
  lab <- case_when(u == "b_delavayi" ~ "Buddleja delavayi", u == "unidentified" ~ "Unidentified",
                   u == "other_plants" ~ "Other plants", TRUE ~ gsub("_", " ", sub("^c[0-9]+_", "", u)))
  lab <- paste0(toupper(substr(lab, 1, 1)), substring(lab, 2))
  dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
  lab[dup] <- paste0(lab[dup], " (", toupper(sub("_.*", "", u[dup])), ")")
  lab[match(type, u)]
}

# Compact letter display from pairwise P values (insert-and-absorb, Piepho 2004)
cld_letters <- function(levels, pairs, means, alpha = 0.05) {
  L <- matrix(TRUE, length(levels), 1, dimnames = list(levels, NULL))
  sig <- filter(pairs, p < alpha)
  for (r in seq_len(nrow(sig))) {
    i <- sig$g1[r]; j <- sig$g2[r]
    for (k in rev(which(L[i, ] & L[j, ]))) {
      a <- L[, k]; a[i] <- FALSE; b <- L[, k]; b[j] <- FALSE
      L <- cbind(L[, -k, drop = FALSE], a, b)
    }
    keep <- vapply(seq_len(ncol(L)), function(c) !any(vapply(seq_len(ncol(L)), function(d)
      d != c && all(L[, c] <= L[, d]) && (any(L[, c] < L[, d]) || d < c), logical(1))), logical(1))
    L <- L[, keep, drop = FALSE]
  }
  first <- apply(L, 2, function(col) min(match(levels[col], levels[order(means[levels])])))
  L <- L[, order(first), drop = FALSE]
  vapply(levels, function(g) paste(letters[which(L[g, ])], collapse = ""), character(1))
}

# Figure style ----
COL_GROUP <- c(Bees = "#D9A94A", Butterflies = "#9A8BB5", Others = "#506B5C")
COL_TRT   <- c(Control = "#DCDCDC", Supplemental = "#8FA9C2")
LINE_TRT  <- c(Control = "#7A7A7A", Supplemental = "#4A6378")
INK <- "#1A1A1A"; INK2 <- "#4D4D4D"; MUTED <- "#8C8C8C"; LINE <- "#3D3D3D"

theme_pub <- function(size = 9) {
  theme_classic(base_size = size, base_family = "sans") +
    theme(text = element_text(colour = INK), axis.text = element_text(colour = INK2, size = size - 1),
          axis.line = element_line(colour = INK2, linewidth = 0.35), axis.ticks = element_line(colour = INK2, linewidth = 0.35),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = size),
          legend.position = "top", legend.title = element_blank(), legend.key.size = unit(3.2, "mm"),
          legend.text = element_text(size = size - 1), plot.tag = element_text(face = "bold", size = size + 2),
          plot.background = element_rect(fill = "white", colour = NA))
}
save_fig <- function(p, name, w = 180, h = 110) {
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", bg = "white")
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
}
