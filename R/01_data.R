source("R/00_setup.R", local = TRUE)

experiment  <- read_csv(file.path(DIR_DATA, "pollination_experiment.csv"))
visitation  <- read_csv(file.path(DIR_DATA, "visitation.csv"))
disturbance <- read_csv(file.path(DIR_DATA, "disturbance.csv"))
pollen      <- read_csv(file.path(DIR_DATA, "pollen_loads.csv"))

# Input checks ----
check <- function(ok, msg) if (!all(ok)) stop(msg, call. = FALSE)
unique_key <- function(d, ...) !anyDuplicated(select(d, ...))
check(unique_key(disturbance, year, plant_id), "disturbance.csv: a plant-year appears twice")
check(unique_key(visitation, year, plant_id, visitor_group), "visitation.csv: a plant-year and visitor group appear twice")
check(unique_key(experiment, year, plant_id, pair, treatment), "pollination_experiment.csv: a unit appears twice")
check(unique_key(pollen, insect_no), "pollen_loads.csv: insect_no is not unique")
check(count(visitation, year, plant_id)$n == length(GROUPS), "visitation.csv: each plant-year needs one row per visitor group")
check(count(experiment, year, plant_id, pair)$n == 2, "pollination_experiment.csv: each pair needs a control and a supplemental unit")
check(nrow(anti_join(visitation, disturbance, by = c("year", "plant_id", "field_label"))) == 0 &&
      nrow(anti_join(disturbance, visitation, by = c("year", "plant_id", "field_label"))) == 0,
      "visitation.csv and disturbance.csv must hold the same plant-years and labels")
check(nrow(anti_join(experiment, disturbance, by = c("year", "plant_id", "field_label"))) == 0,
      "pollination_experiment.csv: plant-year or label not found in disturbance.csv")
check(c(visitation$visitor_group, pollen$visitor_group) %in% GROUPS, "visitor_group must be Bees, Butterflies or Others")
check(experiment$treatment %in% c("Control", "Supplemental"), "treatment must be Control or Supplemental")
check(as.matrix(select(disturbance, road:cutting_observer2)) %in% 0:2, "disturbance scores must be 0, 1 or 2")
check(with(experiment, fruits <= flowers & seeds <= ovules), "pollination_experiment.csv: fruits > flowers or seeds > ovules")
POLLEN_TYPES <- names(pollen)[match("b_delavayi", names(pollen)):ncol(pollen)]
check(c(as.matrix(experiment[, c("flowers", "fruits", "ovules", "seeds")]), visitation$visits,
        as.matrix(pollen[, POLLEN_TYPES])) >= 0, "counts must not be negative")

# Experimental units ----
units <- experiment |>
  mutate(plant_year = paste(plant_id, year, sep = "_"), pair_id = paste(plant_year, pair, sep = "_"),
         treatment = factor(treatment, c("Control", "Supplemental")),
         fruits_failed = flowers - fruits, ovules_unfilled = ovules - seeds)

# Plant-years ----
ratio <- function(num, den, trt, level) sum(num[trt == level]) / sum(den[trt == level])
reproduction <- units |>
  group_by(year, plant_id, plant_year) |>
  summarise(n_pairs = n_distinct(pair),
            fs_c = ratio(fruits, flowers, treatment, "Control"), fs_s = ratio(fruits, flowers, treatment, "Supplemental"),
            ss_c = ratio(seeds, ovules, treatment, "Control"),   ss_s = ratio(seeds, ovules, treatment, "Supplemental"),
            .groups = "drop") |>
  mutate(pli_fs = 1 - fs_c / fs_s, pli_ss = 1 - ss_c / ss_s,
         lrr_fs = log(fs_s / fs_c), lrr_ss = log(ss_s / ss_c))

plant_years <- disturbance |>
  transmute(year, plant_id, field_label, plant_year = paste(plant_id, year, sep = "_"),
            road, herbivory = herbivory_observer2, cutting = cutting_observer2,
            index = road + herbivory + cutting, damage = herbivory + cutting,
            index_field = road + herbivory_observer1 + cutting_observer1) |>
  left_join(visitation |> select(year, plant_id, flowers_observed, visitor_group, visits) |>
              pivot_wider(names_from = visitor_group, values_from = visits) |> rename_with(tolower, any_of(GROUPS)),
            by = c("year", "plant_id")) |>
  left_join(select(reproduction, -plant_year), by = c("year", "plant_id")) |>
  mutate(in_experiment = !is.na(n_pairs), repeated = plant_id %in% plant_id[duplicated(plant_id)])

units <- left_join(units, select(plant_years, plant_year, road, herbivory, cutting, index, damage, index_field), by = "plant_year")

# Visitation ----
vis <- visitation |>
  mutate(plant_year = paste(plant_id, year, sep = "_"), visitor_group = factor(visitor_group, GROUPS),
         rate = 100 * visits / flowers_observed) |>
  left_join(select(plant_years, plant_year, index), by = "plant_year")

# Pollen loads ----
pollen <- pollen |>
  mutate(visitor_group = factor(visitor_group, GROUPS), bd = b_delavayi,
         total = rowSums(across(all_of(POLLEN_TYPES))), carries_bd = bd > 0,
         bd_share = if_else(total > 0, bd / total, NA_real_))
