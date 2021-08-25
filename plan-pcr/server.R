library(shiny)
library(tidyverse)
library(readr)

server <- function(input, output) {

  # Final RNA concentration, determined by protocol
  final_rna_conc <- 5#ng/uL
  
  # Perform in triplicate
  replicates <- 3
  
  # Extra, since nothing is ever perfect
  safety_replicates <- 6
  
  # Volume of RNA per well
  rna_per_well <- 2#uL
  
  # Final RNA volume per sample
  final_rna_vol <- function(num_primers) {
    ((num_primers * replicates) + safety_replicates) * rna_per_well
  }
  
  rna_dil_factor <- function(initial_rna_conc, num_primers) {
    vol_to_add <- final_rna_conc*final_rna_vol(num_primers)/initial_rna_conc
    best_factor <- 1
    if (vol_to_add < 1) {
      exact_factor <- 1/vol_to_add
      best_factor <- ceiling(exact_factor/5)*5 # Give something divisible by 5
    }
    best_factor
  }
  
  output$rna_table <- renderTable({
    rna <- read_tsv(input$rna_data$datapath, show_col_types = FALSE)
    if(ncol(rna) == 1) {
      rna <- rna |> 
        mutate(names = paste0("Sample ", 1:nrow(rna))) |> 
        relocate(names)
    }
    rna |> 
      setNames(c("names", "conc")) |> 
      rowwise() |> 
      mutate(dilution_factor = rna_dil_factor(conc, input$primers),
             diluted_concentration = conc/dilution_factor,
             final_vol = final_rna_vol(input$primers),
             diluted_rna_to_add = final_rna_conc * final_vol / diluted_concentration,
             water_to_add = final_vol - diluted_rna_to_add)
  })
  
  output$sample_layout <- renderPlot({
    if(input$plate_format == "96_well") {
      plate_rows <- 8
      plate_cols <- 12
      size <- 16
    } else if(input$plate_format == "384_well") {
      plate_rows <- 16
      plate_cols <- 24
      size <- 8
    }
    empty_plate <- expand_grid(row = letters[1:plate_rows], col = 1:plate_cols) |> 
      mutate(row = as.factor(row) |> fct_rev(),
             col = as.factor(col))
    ggplot(empty_plate, aes(x = col, y = row)) + geom_point(size = size)
  }, width = 600)
  # Should eventually find a way to preserve given names in df
  
  # Might be as simple/rudimentary as just caching the names and resetting them
  # at the end
  
  # Should delay render table until input
  
  # Primer name input too
  # Be easy with primer naming. If too many primer names supplied, trim off the end
  # If too few, fill out the rest with dummy names.
  
  # Should probably calculate final_rna_vol once (only when input$primers
  # changes) and then supply that elsewhere
  
  # Num of max primers set too low. Should depend on num samples.
  
  # Num of max primers should dep on whether border is excluded or not.
  
  
}
