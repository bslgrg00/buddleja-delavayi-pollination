source("R/01_data.R", local = TRUE)

# Benjamini-Hochberg correction across the tests of years, pollen loads and disturbance ----
bh <- bind_rows(load_table("year_tests"), load_table("pollen_load_tests"), load_table("disturbance_tests")) |>
  filter(bh_family) |>
  select(test, any_of(c("effect", "estimate", "lwr", "upr", "chisq", "F", "df")), p) |>
  mutate(q = p.adjust(p, "BH"))
save_table(bh, "tests_BH")

# Identifiers and experimental pairs per plant ----
save_table(plant_years |> select(plant_id, year, field_label, n_pairs) |>
             pivot_wider(names_from = year, values_from = c(field_label, n_pairs)) |> arrange(plant_id),
           "plant_identifiers")

writeLines(capture.output(sessionInfo()), file.path("outputs", "session_info.txt"))
