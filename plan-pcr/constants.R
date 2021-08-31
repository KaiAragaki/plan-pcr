# Experimental Constants

# Non-targeting control - add one to sample number
ntc <- 1

# Final [RNA], determined by protocol
final_rna_conc <- 5#ng/uL

# Perform in triplicate
reps <- 3

# Extra, since nothing is ever perfect
safety_reps <- 6

# Vol RNA/well
rna_per_well <- 2#uL 


# Plate Templates --------------------------------------------------------------

plate_template <- function(n_row, n_col) {
  expand_grid(col = 1:n_col,
              row = letters[1:n_row]) |> 
    mutate(row = as.factor(row) |> fct_rev(),
           col = as.factor(col))
}

plate_96 <- plate_template(8, 12)
plate_384 <- plate_template(16, 24)


## Borderless plates -----------------------------------------------------------

plate_96_nb <- plate_96 |> 
  filter(!(col %in% c(1, 12) | row %in% c("a", "h")))
plate_384_nb <- plate_384 |> 
  filter(!(col %in% c(1, 24) | row %in% c("a", "p")))
