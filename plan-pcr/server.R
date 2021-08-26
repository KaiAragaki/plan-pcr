library(shiny)
library(tidyverse)
library(readr)

server <- function(input, output) {

  # Plate Templating -----------------------------------------------------------
  
  
  ## Empty plate templates -----------------------------------------------------
  # Upfront single initialization makes sense, minimal overhead
  make_plate_template <- function(n_row, n_col) {
    expand_grid(col = 1:n_col,
                row = letters[1:n_row]) |> 
      mutate(row = as.factor(row) |> fct_rev(),
             col = as.factor(col))
  }
  
  empty_plate_96 <- make_plate_template(8, 12)
  empty_plate_384 <- make_plate_template(16, 24)
  
  ### No border plates ---------------------------------------------------------
  empty_plate_96_nb <- empty_plate_96 |> 
    filter(!(col %in% c(1, 12) | row %in% c("a", "h"))) |> 
    mutate(across(everything(), fct_drop))
  empty_plate_384_nb <- empty_plate_384 |> 
    filter(!(col %in% c(1, 24) | row %in% c("a", "p"))) |> 
    mutate(across(everything(), fct_drop))
  
  # Lanes are horizontal and vertical tracts of wells that (in combination) form
  # a checkerboard pattern that outlines where primers and samples will go
  
  # Vertical lanes are experiment independent and can be made outright
  make_lanes_v <- function(plate) {
    num_lanes <- length(levels(plate$col)) %/% replicates
    wells_per_lane <- length(levels(plate$row)) * replicates
    lane <- rep(1:num_lanes, each = wells_per_lane) |> 
      as.factor()
    length(lane) <- nrow(plate)
    plate <- arrange(plate, col, desc(row))
    plate$lane_v <- lane
    plate
  }
  
  # Horizontal lanes depend on the number of samples and need to be calculated after user input.
  make_lanes_h <- function(plate) {
    num_lanes <- length(levels(plate$row)) %/% (n_samples() + ntc)
    if (num_lanes < 1){
      # do stuff
    } else {
      wells_per_lane <- length(levels(plate$col)) * (n_samples() + ntc)
      lane <- rep(1:num_lanes, each = wells_per_lane) |> 
        as.factor()
      length(lane) <- nrow(plate)
      plate <- arrange(plate, desc(row), col)
      plate$lane_h <- lane
    }
    plate
  }
    
    # Some instances that may be a problem:
    # - One primer, samples > nrow plate
    ## - Fairly simple solution - catch if n_samples > nrow plate and allow for wrapping
    ## Since it won't be tidy regardless, we can just wrap it all the way (no need to shift over after one primer has finished)
    # - many (> ncol_plate %/% 3) primers, all over 0.5 * nrow_plate but  < nrow_plate.
    ## More complex to catch, need to make a design decision to figure out what to do.
    ## In real life I would do the first (<= ncol_plate %/% 3) primers normally, then I would wrap the last one.
    # Many, many primers with only one or two samples
    ## Honestly this can probably just follow the same strategy as above: Move laterally until you can't anymore, then wrap the bottom ones.

  
  
  
  
  # Experimental constants -----------------------------------------------------
  
  # Non-targeting control - add one to samples
  ntc <- 1
  
  # Final RNA concentration, determined by protocol
  final_rna_conc <- 5#ng/uL
  
  # Perform in triplicate
  replicates <- 3
  
  # Extra, since nothing is ever perfect
  safety_replicates <- 6
  
  # Volume of RNA per well
  rna_per_well <- 2#uL
  
  # Final RNA volume per sample
  final_vol <- reactive({
    x <- ((input$primers * replicates) + safety_replicates) * rna_per_well
    as.integer(x)
  })
  
  # Read in data
  rna <- reactive({
    req(input$rna_data)
    x <- read_tsv(input$rna_data$datapath, show_col_types = FALSE)
    if(ncol(x) == 1) {
      x <- x |> 
        mutate(names = paste0("Sample ", 1:nrow(x))) |> 
        relocate(names)
    }
    x |> 
      setNames(c("names", "conc"))
  })
  
  get_best_factor <- function(vol_to_add) {
    if (vol_to_add < 1) {
      exact_factor <- 1 / vol_to_add
      best_factor <- ceiling(exact_factor/5) * 5 # Give something divisible by 5
    } else {
      best_factor <- 1
    }
    as.integer(best_factor)
  }
  
  rna_dil_factor <- reactive({
    x <- rna() |> 
      rowwise() |> 
      mutate(vol_to_add = final_rna_conc * final_vol() / conc,
             dilution_factor = get_best_factor(vol_to_add))
    x$dilution_factor
  })

  rna_table <- reactive({
    rna() |> 
      mutate(dilution_factor = rna_dil_factor(),
             diluted_concentration = conc/dilution_factor,
             final_vol = final_vol(),
             diluted_rna_to_add = final_rna_conc * final_vol / diluted_concentration,
             water_to_add = final_vol - diluted_rna_to_add)
  })
  
  n_samples <- reactive({
    nrow(rna_table())
  })
  
  # Format wells, then bind samples on the side
  # Might not be the best choice in terms of reactivity...
  # Do this as a working model then make it sexy later
  
  
  format_wells <- function(plate_structure) {
    if(input$exclude_border) {
      
    }
  }
  
  format_samples <- function(sample_table) {
    
  }
  
  
  
  output$rna_table <- renderTable({
    rna_table()
  })
  
  output$sample_layout <- renderPlot({
    if(input$plate_format == "96_well") {
      plate <- empty_plate_96 |> 
        make_lanes_h() |> 
        make_lanes_v() |> 
        arrange(lane_h, lane_v) |> 
        rowwise() |> 
        mutate(section = if_else(is.na(lane_h) || is.na(lane_v), NA_character_, paste0(lane_h, ", ", lane_v))) |> 
        group_by(section) |> 
        mutate(id = cur_group_id(),
               section = if_else(id > input$primers, NA_character_, section))
        
      size = 16
      
    } else {
      plate <- empty_plate_384 |> 
        make_lanes_h() |> 
        make_lanes_v() |> 
        rowwise() |> 
        mutate(section = if_else(is.na(lane_h) || is.na(lane_v), NA_character_, paste0(lane_h, ", ", lane_v))) |> 
        group_by(section) |> 
        mutate(id = cur_group_id(),
               section = if_else(id > input$primers, NA_character_, section))
      size = 8
    }
    ggplot(plate, aes(x = col, y = row, color = section)) + geom_point(size = size) + theme(legend.position = "none")
  }, width = 600)
  # Should eventually find a way to preserve given names in df
  
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end
  
  # Primer name input too
  # Be easy with primer naming. If too many primer names supplied, trim off the end
  # If too few, fill out the rest with dummy names.
  
  # Num of max primers set too low. Should depend on num samples.
  
  # Num of max primers should dep on whether border is excluded or not.
  
  # Error if too many samples * primers for plate selected
  
}
