ui <- fluidPage(
  titlePanel("PCR Experiment Planning"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("rna_data",
                "Upload RNA Concentrations",
                width = 250),
      numericInput("primers",
                  "Number of Primers:",
                  min = 1,
                  max = 64,
                  value = 2,
                  step = 1, 
                  width = 250),
      radioButtons("plate_format",
                   "Plate Format:",
                   c("96 Well" = "96_well", "384 Well" = "384_well")),
      radioButtons("exclude_border",
                   "Exclude Plate Border?",
                   c("Yes" = TRUE, "No" = FALSE), 
                   selected = TRUE),
      textInput("primer_names", "Primer Names (optional)", placeholder = "Primer 1; Primer 2..."),
      downloadButton("report",
                     "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sample Preparation", tableOutput("sample_prep")),
        tabPanel("Mastermix Preparation", tableOutput("mm_prep")),
        tabPanel("Mastermix Layout", plotOutput("mm_layout")),
        tabPanel("Sample Layout", plotOutput("sample_layout"))
      )
    )
  )
)
# Need to validate file ext (allow for tsv, csv, excel files)
# See https://shiny.rstudio.com/reference/shiny/latest/fileInput.html

# Allow user to supply code to select which samples they want? Security risk?
# Maybe just allow for a simple click to select?
