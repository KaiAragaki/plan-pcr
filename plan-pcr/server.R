library(shiny)
library(tidyverse)
library(readr)
library(readxl)
library(sass)
library(gt)
library(mop)
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
  
  # Sample Name ----------------------------------------------------------------
  sample_names <- reactive({
    sn <- rna()$names
    if (input$sample_names != "") {
      user_names <- unlist(strsplit(input$sample_names, split = "; ", ))
      if (length(user_names) > length(sn)) {
        user_names <- user_names[1:n_samples()]
      }
      sn[1:length(user_names)] <- user_names
    }
    sn 
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
    if(!isTruthy(input$rna_data)) {
      x <- read_tsv("./data/example_samples.tsv")
    } else {
      validate(need(ext %in% c("tsv", "csv", "xls", "xlsx"), "Only .tsv, .csv, .xls(x) files supported"))
      if (ext == "tsv"){
        x <- read_tsv(input$rna_data$datapath, show_col_types = FALSE)
      } else if (ext == "csv"){
        if (guess_encoding(input$rna_data$datapath)$encoding[1] == "UTF-16LE") {
          x <- read_nanodrop(input$rna_data$datapath) |> tidy_lab() |> scrub()
          if (input$use_adj_conc) {
            x <- mutate(
              x,
              conc = if_else(is.na(corrected_ngul), conc, corrected_ngul)
            )
          }
          x <- select(x, sample_name, conc)
        } else {
          x <- read_csv(input$rna_data$datapath, show_col_types = FALSE)
        }
      } else {
        x <- read_excel(input$rna_data$datapath)
      }

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
      mutate(names = sample_names(),
        dilution_factor = rna_dil_factor(),
             diluted_concentration = conc/dilution_factor,
             final_vol = final_vol(),
             diluted_rna_to_add = final_rna_conc * final_vol / diluted_concentration,
             water_to_add = final_vol - diluted_rna_to_add) |> 
      relocate(final_vol, .after = last_col()) |> 
      gt() |> 
      fmt_number(
        c("diluted_rna_to_add", "water_to_add")
      ) |> 
      cols_label(
        names = "Sample",
        conc = "[RNA]",
        dilution_factor = "Dil. Factor",
        diluted_concentration = "Dil. [RNA] ",
        diluted_rna_to_add = "Vol. dil. [RNA]",
        water_to_add = "Vol. H2O",
        final_vol = "Vol. final (uL)"
      ) |> 
      tab_style(
        style = list(
          cell_fill(color = "#2C5295"),
          cell_text(color = "#FFFFFF"),
          cell_borders(color = "#1C345F")
        ),
        locations = cells_body()
      ) |>  
      tab_style(
        style = list(
          cell_fill(color = "#1C345F"),
          cell_text(weight = "bold", color = "#FFFFFF"),
          cell_borders(color = "#1C345F")
        ),
        locations = cells_column_labels()
      ) |> 
      tab_options(table.align = "left") |> 
      tab_options(
        column_labels.border.bottom.color = "#1C345F"
      ) |> 
      opt_table_outline(color = "#1C345F")
  })
  
  output$sample_prep <- render_gt({
    sample_prep()
  })
  
  # Mastermix Preparation ------------------------------------------------------
  mm_prep <- reactive({
    mm() |> 
      gt() |> 
      fmt_number(
        c("vol")
      ) |> 
      cols_label(
        reagent = "Reagent",
        vol = "Volume (uL)"
      ) |> 
      tab_style(
        style = list(
          cell_fill(color = "#2C5295"),
          cell_text(color = "#FFFFFF"),
          cell_borders(color = "#1C345F")
        ),
        locations = cells_body()
      ) |>  
      tab_style(
        style = list(
          cell_fill(color = "#1C345F"),
          cell_text(weight = "bold", color = "#FFFFFF"),
          cell_borders(color = "#1C345F")
        ),
        locations = cells_column_labels()
      ) |> 
      tab_options(table.align = "left") |> 
      tab_options(
        column_labels.border.bottom.color = "#1C345F"
      ) |> 
      opt_table_outline(color = "#1C345F")
  })

  output$mm_prep <- render_gt({
    mm_prep()
  })
  
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
  
  mm_layout <- reactive({
    req(input$primers)
    is_over()
    ggplot(plate_final(), aes(x = col, y = row, color = primer, label = primer_name)) + 
      geom_point(size = size()) + 
      geom_text(color = "black", size = size_text(), na.rm = TRUE) + 
      theme(legend.position = "none", plot.background = element_blank())
  })
  
  output$mm_layout <- renderPlot({
    mm_layout()  
  }, width = 800, height = 500)
  
  # Sample Layout --------------------------------------------------------------
  
  sample_layout <- reactive({
    req(input$primers)
    is_over()
    ggplot(plate_final(), aes(x = col, y = row, color = sample, label = sample_name)) + 
      geom_point(size = size()) + 
      geom_text(color = "black", size = size_text(), na.rm = TRUE) +
      theme(legend.position = "none", plot.background = element_blank())
  })
  
  output$sample_layout <- renderPlot({
    sample_layout()
  }, width = 800, height = 500) 
  
  
  # Should eventually find a way to preserve given names in df
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end
  
  # Truncate sample names that are too long for plot (stringr)
  
  # Cycle through color limited distinct colors
  
  # Integration with Plate? Allow for a late with specific layers to be loaded
  # in, autocalculation?
  
  make_filename <- reactive({
    if(isTruthy(input$rna_date$name)){
      file_name <- input$rna_data$name |> 
        str_replace( "\\..*$", "_report.html")
    } else(
      file_name <- "example_report.html"
    )

  })
  
  output$get_report <- downloadHandler(
    filename = function(){
      make_filename()
    },
    content = function(file) {
      temp_report = file.path("report_template.rmd")
      params = list(file = input$rna_data,
                    sample_prep = sample_prep(),
                    mm_prep = mm_prep(),
                    mm_layout = mm_layout(),
                    sample_layout = sample_layout(),
                    primer_number = input$primers,
                    primer_names = input$primer_names,
                    plate_format = input$plate_format,
                    exclude_border = input$exclude_border)
      rmarkdown::render(temp_report, output_file = file, 
                        params = params, envir = new.env(parent = globalenv()))
    })
}
