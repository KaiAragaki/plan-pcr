server <- function(input, output) {
  
  output$rna_table <- renderTable({
    read_tsv(input$rna_data$datapath, show_col_types = FALSE)
  })
}
