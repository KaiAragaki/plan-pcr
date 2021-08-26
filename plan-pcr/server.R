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
  replicates <- 3
  
  # Extra, since nothing is ever perfect
  safety_replicates <- 6
  
  # Vol RNA/well
  rna_per_well <- 2#uL 
  
  # Final RNA volume per sample
  final_vol <- reactive({
    x <- ((input$primers * replicates) + safety_replicates) * rna_per_well
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
    filter(!(col %in% c(1, 12) | row %in% c("a", "h"))) |> 
    mutate(across(everything(), fct_drop))
  plate_384_nb <- plate_384 |> 
    filter(!(col %in% c(1, 24) | row %in% c("a", "p"))) |> 
    mutate(across(everything(), fct_drop))
  
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

  n_plate_cols <- reactive({
    plate <- plate_init()
    length(levels(plate$col))
  })
  
  n_plate_rows <- reactive({
    plate <- plate_init()
    length(levels(plate$row))
  })
  
  ### Vertical lanes -----------------------------------------------------------
  # Vertical lanes only depends on plate type (+ border status)
  plate_vlane <- reactive({
    num_lanes <- n_plate_cols() %/% replicates
    wells_per_lane <- n_plate_rows() * replicates
    lane <- rep(1:num_lanes, each = wells_per_lane) |> 
      as.factor()
    length(lane) <- nrow(plate_init())
    plate <- arrange(plate_init(), col, desc(row))
    plate$lane_v <- lane
    plate
  })

  
  ### Horizontal lanes ---------------------------------------------------------
  # Horizontal lanes depend on the number of samples and need to be calculated
  # after user input.
  plate_hlane <- reactive({
    num_lanes <- n_plate_rows() %/% (n_samples() + ntc)
    if (num_lanes < 1){
      stop("oops")
    } else {
      wells_per_lane <- n_plate_cols() * (n_samples() + ntc)
      lane <- rep(1:num_lanes, each = wells_per_lane) |> 
        as.factor()
      length(lane) <- nrow(plate_vlane())
      plate <- arrange(plate_vlane(), desc(row), col)
      plate$lane_h <- lane
    }
    plate
  })

    # Some instances that may be a problem:
    # - One primer, samples > nrow plate
    ## - Fairly simple solution - catch if n_samples > nrow plate and allow for wrapping
    ## Since it won't be tidy regardless, we can just wrap it all the way (no need to shift over after one primer has finished)
    # - many (> ncol_plate %/% 3) primers, all over 0.5 * nrow_plate but  < nrow_plate.
    ## More complex to catch, need to make a design decision to figure out what to do.
    ## In real life I would do the first (<= ncol_plate %/% 3) primers normally, then I would wrap the last one.
    # Many, many primers with only one or two samples
    ## Honestly this can probably just follow the same strategy as above: Move laterally until you can't anymore, then wrap the bottom ones.


  ## Make sections -------------------------------------------------------------
  # Sections are the checkerboard-like patterns made by the grid of lanes
  plate_section <- reactive({
    if(input$plate_format == "96_well") {
      full_plate <- plate_96
    } else {
      full_plate <- plate_384
    }
    plate <- plate_hlane() |> 
      arrange(lane_h, lane_v) |> 
      rowwise() |> 
      mutate(section = if_else(is.na(lane_h) || is.na(lane_v), 
                               NA_character_, 
                               paste0(lane_h, ", ", lane_v))) |> 
      group_by(section) |> 
      mutate(id = cur_group_id(),
             section = if_else(id > input$primers, NA_character_, section),
             available_well = TRUE)
    if(input$exclude_border) {
      plate <- left_join(full_plate, plate)
    }
    plate
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
  
  ## Get number of samples -----------------------------------------------------
  n_samples <- reactive({
    nrow(rna())
  })
  
  ## Get sections size ---------------------------------------------------------
  section_size <- reactive({
    replicates * (n_samples() + ntc)
  })
  
  max_sections_before_flow <- reactive({
    plate <- plate_section()
    max_wide <- group_by(plate, col) |> 
      summarize(sum = sum(available_well, na.rm = TRUE),
                available_col = sum > 0)
    max_wide <- sum(max_wide$available_col) %/% replicates
    max_tall <- group_by(plate, row) |> 
      summarize(sum = sum(available_well, na.rm = TRUE),
                available_row = sum > 0)
    max_tall <- sum(max_tall$available_row) %/% (n_samples() + ntc)
    max <- max_wide * max_tall
    print(paste0("msbf = ", max))
  })
  
  max_sections_after_flow <- reactive({
    plate <- plate_section()
    max_wide <- group_by(plate, col) |> 
      summarize(sum = sum(available_well, na.rm = TRUE),
                available_col = sum > 0)
    max_wide <- sum(max_wide$available_col) %/% replicates
    max_tall <- group_by(plate, row) |> 
      summarize(sum = sum(available_well, na.rm = TRUE),
                available_row = sum > 0)
    max <- (max_wide * sum(max_tall$available_row)) %/% (n_samples() + ntc)
    print(paste0("msaf = ", max))
  })
  
  max_sections_theoretical <- reactive({
    x <- sum(plate_section()$available_well, na.rm = TRUE) %/% section_size()
    print(x)
  })
  
  
  ## Check if current layout exceeds max # sections ----------------------------
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
    is_over()
    max_sections_before_flow()
    max_sections_after_flow()
    max_sections_theoretical()
    if(input$plate_format == "96_well") {
      full_plate <- plate_96
      size = 16
    } else {
      full_plate <- plate_384
      size = 8
    }
    ggplot(plate_section(), aes(x = col, y = row, color = section)) + geom_point(size = size) + theme(legend.position = "none")
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
  
  # Need catch for too many primers (should be in lanes_v or plotting, not sure.
  # Probably plotting though, lanes_v fills out as many as possible and then
  # stops)
  
  # Error if too many samples * primers for plate selected
 
  # Current implementation uses sections which will not work when I eventually
  # break down the idea of sections using wraparounds when sample numbers get
  # high with few primers.
   
  # Should just see if it can fit using a pure area calculation
 
  #If it can't stack side by side anymore (see in case of 384, no plate border 6
  # samples with 18 primers) then the next one should trigger a 'wrapping' event.
  ## The problem with this is that it almost requires a means of knowing about
  ## contiguity. 
  # One way of getting around doing this is to measure how big a
  # section will be, then calculate the largest rectangular section that can be
  # made from that (ie it can be two sections tall and 8 sections wide, in the
  # case above). Once that total number of sections is exceeded, switch to a flow
  # strategy.
  
  # A case (which likely makes up a family of cases) that I think it's ok to not
  # rearrange everything in order to make it all fit, even though it could:
  
  ## 10 -> 11 primers, 96 well, excluding border, 1 sample + ntc
   
  ## This can probably be caught by preventing changing the max number of wells
  ## wide, but allowing for flexibility of section hight (or flow or whatever)
  
  
  # MSBF and MSAF only need to be calculated once per input
}
