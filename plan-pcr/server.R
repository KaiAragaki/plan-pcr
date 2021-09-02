library(shiny)
library(tidyverse)
library(readr)
library(readxl)
library(sass)
source("constants.R")
source("utils.R")
source("denote_lane.R")


server <- function(input, output) {
  
  # Final RNA volume per sample
  final_vol <- reactive({
    (n_primers()*reps + safety_reps) * rna_per_well |> 
      as.integer()
  })
  
  # 'Full plate' pronoun -------------------------------------------------------
  full_plate <- reactive({
    if(input$plate_format == "96_well") plate_96 else plate_384
  })
  
  # Initialize plate -----------------------------------------------------------
  this_plate <- reactive({
    if (input$plate_format == "96_well") {
      if(input$exclude_border) plate_96_nb else plate_96
    } else {
      if(input$exclude_border) plate_384_nb else plate_384
    }
  })
  
  # Primer Names ---------------------------------------------------------------
  primer_names <- reactive({
    pn <- paste("Primer", 1:n_primers())
    if (input$primer_names != "") {
      user_names <- unlist(strsplit(input$primer_names, split = "; ", ))
      user_names <- user_names[1:n_primers()]
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
    n_primers() > max_sections_before_flow()
  }) 
  
  ### Vertical lanes -----------------------------------------------------------
  # Vertical lanes only depends on plate type and border status
  plate_vlane <- reactive({
    denote_vlane(this_plate())
  })
  
  ### Horizontal lanes ---------------------------------------------------------
  # Horizontal lanes depend on the number of samples and need to be calculated
  # after user input.
  plate_hlane <- reactive({
    denote_hlane(plate_vlane(), n_samples())
  })
  
  ### Flow lanes ---------------------------------------------------------------
  plate_flow <- reactive({
    flow_lanes(plate_vlane(), n_primers(), n_samples())
  })
  
  ### Section lanes ------------------------------------------------------------
  plate_sect <- reactive({
    section_lanes(plate_hlane(), n_primers())
  })
  
  
  ## Determine sectioning method -----------------------------------------------
  plate <- reactive({
    if(should_flow()) plate_flow() else plate_sect()
  })
  
  ### Add Sample Name Labels ---------------------------------------------------
  
  plate_sn <- reactive({
    add_sample_names(plate(), sample_names())
  })
  
  ### Add Primer Name Labels ---------------------------------------------------
  
  plate_pn <- reactive({
    add_primer_names(plate_sn(), primer_names())
  })
  
  plate_final <- reactive({
    right_join(plate_pn(), full_plate(), by = c("col", "row"))
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
  
  # n_primers ------------------------------------------------------------------
  n_primers <- reactive({
    input$primers
  })
  
  # Usable plate dims ----------------------------------------------------------
  # If the border is excluded, it will not be counted
  plate_dims <- reactive({useable_plate_dims(this_plate())})
  
  
  # section_area ---------------------------------------------------------------
  section_area <- reactive({
    reps * (n_samples() + ntc)
  })
  
  # Max # sections -------------------------------------------------------------
  ## ...before flow -----------------------------------------------------------
  max_sections_before_flow <- reactive({
    max_wide <- plate_dims()$cols %/% reps
    max_tall <- plate_dims()$rows %/% (n_samples() + ntc)
    max_wide * max_tall
  })
  
  ## ...after flow ------------------------------------------------------------
  max_sections_after_flow <- reactive({
    max_wide <- plate_dims()$cols %/% reps
    (max_wide * plate_dims()$rows) %/% (n_samples() + ntc)
  })
  
  ## ...theoretically possible -------------------------------------------------
  max_sections_theoretical <- reactive({
    (plate_dims()$rows * plate_dims()$cols) %/% section_area()
  })
  
  # sample_names ---------------------------------------------------------------
  
  sample_names <- reactive({
    rna()$names
  })
  
  # Checks ---------------------------------------------------------------------
  
  ## User should include control
  
  has_control <- observeEvent(input$primers, {
    shinyFeedback::feedbackWarning("primers", input$primers == 1, 
                                   "Ensure you include an endogenous control, like GAPDH.")
  })
  
  ## Current layout exceeds max # sections -------------------------------------
  is_over <- reactive({
    if(max_sections_theoretical() < n_primers()) {
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
  
  size <- reactiveVal(24)
  size_text <- reactiveVal(8)
  
  observeEvent(input$plate_format, {
    if(input$plate_format == "96_well") {
      size(24)
      size_text(8)
    } else {
      size(12)
      size_text(6)
    }
  })
  
  output$mm_layout <- renderPlot({
    req(input$rna_data)
    req(input$primers)
    is_over()
    ggplot(plate_final(), aes(x = col, y = row, color = primer, label = primer_name)) + 
      geom_point(size = size()) + 
      geom_text(color = "black", size = size_text(), na.rm = TRUE) + 
      theme(legend.position = "none", plot.background = element_blank())
  }, width = 800, height = 500)
  
  # Sample Layout --------------------------------------------------------------
  output$sample_layout <- renderPlot({
    req(input$rna_data)
    req(input$primers)
    is_over()
    ggplot(plate_final(), aes(x = col, y = row, color = sample, label = sample_name)) + 
      geom_point(size = size()) + 
      geom_text(color = "black", size = size_text(), na.rm = TRUE) +
      theme(legend.position = "none", plot.background = element_blank())
    
  }, width = 800, height = 500) 
  
  
  # Should eventually find a way to preserve given names in df
  
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end
  
  
  # Truncate sample names that are too long for plot (stringr)
  
  
  # Cycle through color limited distinct colors
  
  
  # Integration with Plate? Allow for a late with specific layers to be loaded
  # in, autocalculation?
  
  # Report generation and download
  
}
