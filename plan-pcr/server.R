library(shiny)
library(tidyverse)
library(readr)
library(readxl)

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
  
  full_plate <- reactive({
    if(input$plate_format == "96_well") plate_96 else plate_384
  })
  
  
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
  
  
  # Primer Names ---------------------------------------------------------------
  
  primer_names <- reactive({
    pn <- paste("Primer", 1:input$primers)
    if (input$primer_names != "") {
      user_names <- unlist(strsplit(input$primer_names, split = "; ", ))
      user_names <- user_names[1:input$primers]
      pn[1:length(user_names)] <- user_names
    }
    pn 
  })
  
  ## Make lanes ---------------------------------------------------------------- 
  # Lanes are horizontal and vertical tracts of wells that (in combination)
  # form a checkerboard pattern that outlines where primers and samples will
  # go
  
  # If needed, the 'grid structure' will be abandoned in favor of a
  # top-to-bottom, left-to-right filling of wells. This will be employed if
  # using this method can fit more primers, and if more primers are needed to be
  # fit.
  
  should_flow <- reactive({
    input$primers > max_sections_before_flow()
  }) 
  
  
  ### Vertical lanes -----------------------------------------------------------
  # Vertical lanes only depends on plate type (+ border status)
  plate_vlane <- reactive({
    num_lanes <- n_plate_cols() %/% reps
    wells_per_lane <- n_plate_rows() * reps
    lane <- rep(1:num_lanes, each = wells_per_lane)
    length(lane) <- nrow(plate_init())
    arrange(plate_init(), col, desc(row)) |> 
      mutate(lane_v = lane)
  })
  
  
  ### Flowing lanes ------------------------------------------------------------
  plate_flow <- reactive({
    primer <- rep(1:input$primers, each = (n_samples() + ntc) * reps)
    primer_names <- rbind(NA, primer_names(), matrix(nrow = (n_samples() + ntc) * reps - 2, ncol = input$primers)) |> c()
    sample_names <- c(rna()$names, "NTC")
    blanks <- rep(NA, times = n_samples() + ntc)
    sample <- rbind(blanks, sample_names, blanks) |> c() |> rep(times = input$primers)
    sample_color <- rep(sample_names, each = 3) |> c() |>  rep(times = input$primers)
    length(primer) <- nrow(plate_vlane())
    length(primer_names) <- nrow(plate_vlane())
    length(sample) <- nrow(plate_vlane())
    length(sample_color) <- nrow(plate_vlane())
    plate <- arrange(plate_vlane(), lane_v, desc(row)) |> 
      mutate(primer = primer,
             primer_name = primer_names) |> 
      mutate(primer = if_else(primer > input$primers, NA_character_, as.character(primer)),
             primer_name = if_else(is.na(primer), NA_character_, primer_name),
             available_well = TRUE) |>
      arrange(primer) |> 
      mutate(sample = sample,
             sample_color = sample_color) |> 
      right_join(full_plate(), by = c("col", "row"))
    plate
  })
  
  
  ### Horizontal lanes ---------------------------------------------------------
  # Horizontal lanes depend on the number of samples and need to be calculated
  # after user input.
  plate_hlane <- reactive({
    num_lanes <- n_plate_rows() %/% (n_samples() + ntc)
    wells_per_lane <- n_plate_cols() * (n_samples() + ntc)
    lane <- rep(1:num_lanes, each = wells_per_lane)
    length(lane) <- nrow(plate_vlane())
    arrange(plate_vlane(), desc(row), col) |> 
      mutate(lane_h = lane)
  })
  
  
  ### Make sections -------------------------------------------------------------
  # Sections are the checkerboard-like patterns made by the grid of lanes
  plate_section <- reactive({
    primer_names <- rbind(NA, primer_names(), matrix(nrow = (n_samples() + ntc) * reps - 2, ncol = input$primers)) |> c()
    sample_names <- c(rna()$names, "NTC")
    

    blanks <- rep(NA, times = n_samples() + ntc)
    sample <- rbind(blanks, sample_names, blanks) |> c() |> rep(times = input$primers)
    sample_color <- rep(sample_names, each = 3) |> c() |> rep(times = input$primers)
    length(sample) <- nrow(plate_vlane())
    length(sample_color) <- nrow(plate_vlane())
    length(primer_names) <- nrow(plate_vlane())
    plate <- plate_hlane() |> 
      arrange(lane_h, lane_v) |> 
      rowwise() |> 
      mutate(primer = if_else(is.na(lane_h) || is.na(lane_v), 
                              NA_real_, 
                              lane_h * 100 + lane_v)) |> 
      arrange(primer) |> 
      group_by(primer) |> 
      mutate(primer = if_else(cur_group_id() > input$primers, NA_character_, as.character(cur_group_id())),
             available_well = TRUE) |>
      ungroup()

    plate <- plate |> 
      mutate(primer_name = primer_names) |> 
      arrange(primer) |> 
      mutate(sample = sample, 
             sample_color = sample_color) |> 
      right_join(full_plate(), by = c("col", "row"))
    plate
  })
  
  
  ## Finalize Plate ------------------------------------------------------------
  
  plate <- reactive({
    if(should_flow()) plate_flow() else plate_section()
  })
  
  
  # Read in Data ---------------------------------------------------------------
  rna <- reactive({
    ext <- tools::file_ext(input$rna_data$datapath)
    req(input$rna_data)
    validate(need(ext %in% c("tsv", "csv", "xls", "xlsx"), "Only .tsv, .csv, .xls(x) files supported"))
    if (ext == "tsv"){
      x <- read_tsv(input$rna_data$datapath, show_col_types = FALSE)
    } else if (ext == "csv"){
      x <- read_csv(input$rna_data$datapath, show_col_types = FALSE)
    } else {
      x <- read_excel(input$rna_data$datapath)
    }
    if(ncol(x) == 1) {
      x <- x |> 
        mutate(names = paste("Sample", 1:nrow(x))) |> 
        relocate(names)
    }
    x |> 
      setNames(c("names", "conc"))
  })
  
  
  # Mastermix Calculation ------------------------------------------------------
  mm <- reactive({
    tibble(reagent = c("2X RT-PCR Buffer", "Primer", "25X Enzyme Mix", "Nuclease Free H2O"),
           vol = c(6.25, .625, .5, 3.125) * (n_samples() + ntc + 2) * reps)
  })
  
  # DATA ACCESSORS ~~~~~~~~~~~~~~~~~~~~-----------------------------------------
  # n_samples ------------------------------------------------------------------
  n_samples <- reactive({
    nrow(rna())
  })
  
  # Usable plate dims ----------------------------------------------------------
  # If the border is excluded, it will not be counted
  
  ## n_plate_cols --------------------------------------------------------------
  n_plate_cols <- reactive({
    plate <- plate_init()
    length(unique(plate$col))
  })
  
  
  ## n_plate_rows --------------------------------------------------------------
  n_plate_rows <- reactive({
    plate <- plate_init()
    length(unique(plate$row))
  })
  
  # section_area ---------------------------------------------------------------
  section_area <- reactive({
    reps * (n_samples() + ntc)
  })
  
  # Max # sections -------------------------------------------------------------
  ## ...before flow -----------------------------------------------------------
  max_sections_before_flow <- reactive({
    max_wide <- n_plate_cols() %/% reps
    max_tall <- n_plate_rows() %/% (n_samples() + ntc)
    max_wide * max_tall
  })
  
  ## ...after flow ------------------------------------------------------------
  max_sections_after_flow <- reactive({
    max_wide <- n_plate_cols() %/% reps
    (max_wide * n_plate_rows()) %/% (n_samples() + ntc)
  })
  
  ## ...theoretically possible ------------------------------------------------
  max_sections_theoretical <- reactive({
    (n_plate_rows() * n_plate_cols()) %/% section_area()
  })
  
  # Checks ---------------------------------------------------------------------
  
  ## User should include control
  
  has_control <- observeEvent(input$primers, {
    shinyFeedback::feedbackWarning("primers", input$primers == 1, "Ensure you include an endogenous control, like GAPDH.")
  })
  
  ## Current layout exceeds max # sections -------------------------------------
  is_over <- reactive({
    if(max_sections_theoretical() < input$primers) {
      validate("This experiment requires too many wells.")
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
  
  # OUTPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~-----------------------------------------
  
  # Sample Preparation ---------------------------------------------------------
  sample_prep <- reactive({
    req(input$primers)
    rna() |> 
      mutate(dilution_factor = rna_dil_factor(),
             diluted_concentration = conc/dilution_factor,
             final_vol = final_vol(),
             diluted_rna_to_add = final_rna_conc * final_vol / diluted_concentration,
             water_to_add = final_vol - diluted_rna_to_add)
  })
  
  output$sample_prep <- renderTable({
    sample_prep()
  })
  
  # Mastermix Preparation ------------------------------------------------------
  output$mm_prep <- renderTable(
    mm()
  )
  
  # Mastermix Layout -----------------------------------------------------------
  output$mm_layout <- renderPlot({
    req(input$rna_data)
    req(input$primers)
    is_over()
    if(input$plate_format == "96_well") {
      size = 24
      size_text = 8
    } else {
      size = 12
      size_text = 6
    }

    ggplot(plate(), aes(x = col, y = row, color = primer, label = primer_name)) + 
      geom_point(size = size) + 
      geom_text(color = "black", size = size_text, na.rm = TRUE) + 
      theme(legend.position = "none", plot.background = element_blank())
  }, width = 800, height = 500)

  # Sample Layout --------------------------------------------------------------
  output$sample_layout <- renderPlot({
    req(input$rna_data)
    req(input$primers)
    is_over()
    if(input$plate_format == "96_well") {
      size = 24
      size_text = 8
    } else {
      size = 12
      size_text = 6
    }

    ggplot(plate(), aes(x = col, y = row, color = sample_color, label = sample)) + 
      geom_point(size = size) + 
      geom_text(color = "black", size = size_text, na.rm = TRUE) +
      theme(legend.position = "none", plot.background = element_blank())
    
  }, width = 800, height = 500) 
  

  
  # Should eventually find a way to preserve given names in df
  
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end

  
  # Clean up code after going through all that drama
  
  
  # Truncate sample names that are too long for plot (stringr)
  
  
  # Cycle through color limited distinct colors
  
  
  # Integration with Plate? Allow for a late with specific layers to be loaded in, autocalculation?
  
  # Report generation and download
  
  }
