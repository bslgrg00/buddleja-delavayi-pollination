scripts <- c("02_pollen_limitation.R", "03_visitation.R", "04_pollen_loads.R", "05_disturbance.R",
             "06_model_checks.R", "07_figures.R", "08_tables.R")
rscript <- file.path(R.home("bin"), "Rscript")
for (s in scripts) {
  message("Running ", s)
  if (system2(rscript, file.path("R", s)) != 0) stop("Stopped at ", s, call. = FALSE)
}
message("Done. Results are in outputs/")
