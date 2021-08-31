# Denote Plate Lanes

# A 'lane' is a vertical or horizontal tract of wells that in combination define
# 'sections'. These sections are used to denote regions in which primers will
# lie.

# Denote Vertical Lanes

# Vertical Lanes do not require any knowledge of experimental specifics
denote_vlane <- function(plate) {
  dims <- useable_plate_dims(plate)
  n_lanes <- dims$cols %/% reps
  wells_per_lane <- dims$rows * reps
  lanes <- rep(1:n_lanes, each = wells_per_lane)
  length(lanes) <- nrow(plate) # Includes 'unusable' dims (ie borders of borderless)
  arrange(plate, col, desc(row)) |>
    mutate(lane_v = lanes)
}

# Horizontal lanes depend on the number of samples
denote_hlane <- function(plate, n_samples) {
  dims <- useable_plate_dims(plate)
  total_samples <- n_samples + ntc
  num_lanes <-  dims$rows %/% total_samples
  wells_per_lane <- dims$cols * total_samples
  lanes <- rep(1:num_lanes, each = wells_per_lane)
  length(lanes) <- nrow(plate)
  arrange(plate, desc(row), col) |>
    mutate(lane_h = lanes)
}

flow_lanes <- function(plate, n_primers, n_samples) {
  total_samples <- n_samples + ntc
  primers <- rep(1:n_primers, each = total_samples * reps)
  length(primers) <- nrow(plate)
  arrange(plate, lane_v, desc(row)) |>
    mutate(primer = primers,
           primer = if_else(primer <= n_primers,
                            primer,
                            NA_integer_),
           available_well = TRUE)
}

section_lanes <- function(plate, n_primers) {
  plate |>
    mutate(primer = lane_h * 100 + lane_v) |> # Numerical hierarchy for ordering
    arrange(primer) |>
    group_by(primer) |>
    mutate(primer = if_else(cur_group_id() <= n_primers,
                            cur_group_id(),
                            NA_integer_),
           available_well = TRUE) |>
    ungroup()
  str(plate)
  plate
}



## Add names to lanes...
add_sample_names <- function(plate, sample_names) {
  sample_names <- c(sample_names, "NTC")
  plate |>
    group_by(primer) |>
    mutate(sample = (row_number() + 2) %/% reps) |> 
    group_by(sample) |>
    nest() |> 
    ungroup() |> 
    mutate(sample_name = sample_names) |>
    unnest(cols = data) |>
    group_by(sample, primer) |>
    arrange(desc(row), col) |>
    mutate(sample_name = if_else(!is.na(primer) & row_number() == 2, sample_name, NA_character_))
}

add_primer_names <- function(plate, primer_names) {
  plate <- plate |>
    group_by(primer) |>
    nest() |>
    arrange(primer) |>
    ungroup()
  
  length(primer_names) <- nrow(plate) # Account for any NAs (that may or may not be present)
  
  plate |> 
    mutate(primer_name = primer_names) |>
    unnest(cols = data) |>
    group_by(primer) |> 
    arrange(desc(row), col) |>
    mutate(primer_name = if_else(row_number() == 2, primer_name, NA_character_))
}
