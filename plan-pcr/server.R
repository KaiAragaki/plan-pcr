library(shiny)
library(tidyverse)
library(readr)

server <- function(input, output) {

  # Empty plate templates ------------------------------------------------------
  
  # Upfront single initialization makes sense, minimal overhead

  make_plate_template <- function(n_row, n_col) {
    expand_grid(col = 1:n_col,
                row = letters[1:n_row]) |> 
      mutate(row = as.factor(row) |> fct_rev(),
             col = as.factor(col))
  }
  
  empty_plate_96 <- make_plate_template(8, 12)
  empty_plate_384 <- make_plate_template(16, 24)
  
  # No border plates:
  empty_plate_96_nb <- empty_plate_96 |> 
    filter(!(col %in% c(1, 12) | row %in% c("a", "h"))) |> 
    mutate(across(everything(), fct_drop))
  empty_plate_384_nb <- empty_plate_384 |> 
    filter(!(col %in% c(1, 24) | row %in% c("a", "p"))) |> 
    mutate(across(everything(), fct_drop))
  
  # A 'lane' is a tract of wells that will encompass our replicates (three wide,
  # extending down the plate)
  make_lanes <- function(plate) {
    num_lanes <- length(levels(plate$col)) %/% 3
    num_rows <- length(levels(plate$row))
    lane <- rep(paste0("lane_", 1:num_lanes), each = num_rows * 3)
    length(lane) <- nrow(plate)
    plate$lane <- lane
  }
  
  
  
  
  
  
  # Experimental constants -----------------------------------------------------
  
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
      empty_plate <- empty_plate_96
    } else {
      empty_plate <- empty_plate_384
    }
    ggplot(empty_plate, aes(x = col, y = row)) + geom_point(size = size)
  }, width = 600)
  # Should eventually find a way to preserve given names in df
  
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end
  
  # Primer name input too
  # Be easy with primer naming. If too many primer names supplied, trim off the end
  # If too few, fill out the rest with dummy names.
  
  # Num of max primers set too low. Should depend on num samples.
  
  # Num of max primers should dep on whether border is excluded or not.
  
  
}
