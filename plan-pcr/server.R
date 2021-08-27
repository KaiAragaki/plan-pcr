library(shiny)
library(tidyverse)
library(readr)

server <- function(input, output) {
  
  # Experimental Constants -----------------------------------------------------
  
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
  
  # Final RNA volume per sample
  final_vol <- reactive({
    x <- ((input$primers * reps) + safety_reps) * rna_per_well
    as.integer(x)
  })
  
  
  # Plate Templating -----------------------------------------------------------
  
  ## Plate templates -----------------------------------------------------------
  # Upfront single initialization makes sense, minimal overhead
  
  ### Full plates --------------------------------------------------------------
  plate_template <- function(n_row, n_col) {
    expand_grid(col = 1:n_col,
                row = letters[1:n_row]) |> 
      mutate(row = as.factor(row) |> fct_rev(),
             col = as.factor(col))
  }
  
  plate_96 <- plate_template(8, 12)
  plate_384 <- plate_template(16, 24)
  
  
  ### Borderless plates ---------------------------------------------------------
  
  plate_96_nb <- plate_96 |> 
    filter(!(col %in% c(1, 12) | row %in% c("a", "h")))
  plate_384_nb <- plate_384 |> 
    filter(!(col %in% c(1, 24) | row %in% c("a", "p")))
  
  ## Initialize plate ----------------------------------------------------------
  plate_init <- reactive({
    if (input$plate_format == "96_well") {
      if(input$exclude_border) plate_96_nb else plate_96
    } else {
      if(input$exclude_border) plate_384_nb else plate_384
    }
  })
  
  
  ## Make lanes ---------------------------------------------------------------- 
  # Lanes are horizontal and vertical tracts of wells that (in combination)
  # form a checkerboard pattern that outlines where primers and samples will
  # go
  
  should_flow <- reactive({
    input$primers > max_sections_before_flow()
  }) 
  

  
  ### Vertical lanes -----------------------------------------------------------
  # Vertical lanes only depends on plate type (+ border status)
  plate_vlane <- reactive({
    num_lanes <- n_plate_cols() %/% reps
    wells_per_lane <- n_plate_rows() * reps
    lane <- rep(1:num_lanes, each = wells_per_lane) |> as.factor()
    length(lane) <- nrow(plate_init())
    arrange(plate_init(), col, desc(row)) |> 
      mutate(lane_v = lane)
  })
  
  ### Flowing lanes ------------------------------------------------------------
  plate_flow <- reactive({
    if(input$plate_format == "96_well") {
      full_plate <- plate_96
    } else {
      full_plate <- plate_384
    }
    primer <- rep(1:input$primers, each = (n_samples() + ntc) * reps)
    length(primer) <- nrow(plate_vlane())
    plate <- arrange(plate_vlane(), lane_v, desc(row)) |> 
      mutate(primer = primer) |> 
      mutate(primer = if_else(primer > input$primers, NA_character_, as.character(primer)),
             available_well = TRUE) |>
      right_join(full_plate)
    plate
  })
  
  ### Horizontal lanes ---------------------------------------------------------
  # Horizontal lanes depend on the number of samples and need to be calculated
  # after user input.
  plate_hlane <- reactive({
    num_lanes <- n_plate_rows() %/% (n_samples() + ntc)
    wells_per_lane <- n_plate_cols() * (n_samples() + ntc)
    lane <- rep(1:num_lanes, each = wells_per_lane) |> as.factor()
    length(lane) <- nrow(plate_vlane())
    arrange(plate_vlane(), desc(row), col) |> 
      mutate(lane_h = lane)
  })
  
  ## Make sections -------------------------------------------------------------
  # Sections are the checkerboard-like patterns made by the grid of lanes
  plate_section <- reactive({
    if(input$plate_format == "96_well") {
      full_plate <- plate_96
    } else {
      full_plate <- plate_384
    }
    plate_hlane() |> 
      arrange(lane_h, lane_v) |> 
      rowwise() |> 
      mutate(primer = if_else(is.na(lane_h) || is.na(lane_v), 
                              NA_character_, 
                              paste0(lane_h, ", ", lane_v))) |> 
      group_by(primer) |> 
      mutate(primer = if_else(cur_group_id() > input$primers, NA_character_, as.character(cur_group_id())),
             available_well = TRUE) |> 
      right_join(full_plate)
  })
  
  # Read in Data ---------------------------------------------------------------
  rna <- reactive({
    req(input$rna_data)
    x <- read_tsv(input$rna_data$datapath, show_col_types = FALSE)
    if(ncol(x) == 1) {
      x <- x |> 
        mutate(names = paste("Sample", 1:nrow(x))) |> 
        relocate(names)
    }
    x |> 
      setNames(c("names", "conc"))
  })
  
  
  # Data Accessors -------------------------------------------------------------
  ## Get number of samples -----------------------------------------------------
  n_samples <- reactive({
    nrow(rna())
  })
  
  ## Get (useable) dimension of plate ------------------------------------------
  # If the border is excluded, it will not be counted
  
  ### cols ---------------------------------------------------------------------
  n_plate_cols <- reactive({
    plate <- plate_init()
    length(unique(plate$col))
  })
  
  
  ### rows ---------------------------------------------------------------------
  n_plate_rows <- reactive({
    plate <- plate_init()
    length(unique(plate$row))
  })
  
  ## Get sections size ---------------------------------------------------------
  section_size <- reactive({
    reps * (n_samples() + ntc)
  })
  
  ## Get max # sections --------------------------------------------
  ### ...before flow -----------------------------------------------------------
  max_sections_before_flow <- reactive({
    max_wide <- n_plate_cols() %/% reps
    max_tall <- n_plate_rows() %/% (n_samples() + ntc)
    max_wide * max_tall
  })
  
  ### ...after flow ------------------------------------------------------------
  max_sections_after_flow <- reactive({
    max_wide <- n_plate_cols() %/% reps
    (max_wide * n_plate_rows()) %/% (n_samples() + ntc)
  })
  
  ### ...theoretically possible ------------------------------------------------
  max_sections_theoretical <- reactive({
    (n_plate_rows() * n_plate_cols()) %/% section_size()
  })
  
  # Checks ---------------------------------------------------------------------
  
  ## Current layout exceeds max # sections ----------------------------
  is_over <- reactive({
    if(max_sections_theoretical() < input$primers) {
      stop("This experiment requires too many wells.")
    }
  })
  
  
  # Calculate Best Dilution Factor ---------------------------------------------
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
  
  # Make RNA Table -------------------------------------------------------------
  rna_table <- reactive({
    req(input$primers)
    rna() |> 
      mutate(dilution_factor = rna_dil_factor(),
             diluted_concentration = conc/dilution_factor,
             final_vol = final_vol(),
             diluted_rna_to_add = final_rna_conc * final_vol / diluted_concentration,
             water_to_add = final_vol - diluted_rna_to_add)
  })
  
  output$rna_table <- renderTable({
    rna_table()
  })
  
  # Make Sample Layout ---------------------------------------------------------
  output$sample_layout <- renderPlot({
    req(input$rna_data)
    req(input$primers)
    is_over()
    if(input$plate_format == "96_well") {
      full_plate <- plate_96
      size = 16
    } else {
      full_plate <- plate_384
      size = 8
    }
    
    plate <- if(should_flow()) plate_flow() else plate_section()
    
    ggplot(plate, aes(x = col, y = row, color = primer)) + geom_point(size = size) + theme(legend.position = "none")
  }, width = 600)
  
  # Need to find a way to preserve plate structure but still plot without border
  ## Might be as simple as dropping the first and last factors of row and col
  ## You'd have to catch to make sure you weren't putting samples in the 'forbidden NA zone' though
  
  # Should eventually find a way to preserve given names in df
  
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end
  
  # Primer name input too
  # Be easy with primer naming. If too many primer names supplied, trim off the end
  # If too few, fill out the rest with dummy names.

  # Weird one that works but shouldn't:
  
  # 1 sample (2 with ntc), 10 primers, 96 well, no border
  
  # Clean up code after going through all that drama
  
  # Need to be better about not evaling when there isn't a datapoint
  
}
