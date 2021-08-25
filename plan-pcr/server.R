server <- function(input, output) {

  # Final RNA concentration, determined by protocol
  final_rna_conc <- 5#ng/uL
  
  # Perform in triplicate
  replicates <- 3
  
  # Extra, since nothing is ever perfect
  safety_replicates <- 6
  
  # Volume of RNA per well
  rna_per_well <- 2#uL
  
  ##### DEV STATIC VARS
  num_primers <- 3
  initial_rna_conc <- 200
  initial_rna_conc <- 40
  
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
    read_tsv(input$rna_data$datapath, show_col_types = FALSE) |> 
      setNames(c("names", "conc")) |> 
      rowwise() %>%
      mutate(df = rna_dil_factor(conc, input$primers))
  })
  # Should make it such that if input data has only one col, make names
  # Primer name input too
}
