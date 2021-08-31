# Utility Functions


# Usable plate dims ----------------------------------------------------------
# If the border is excluded, it will not be counted

## plate_dims ----------------------------------------------------------------
useable_plate_dims <- function(plate) {
  list(rows = length(unique(plate$row)),
       cols = length(unique(plate$col)))
}
