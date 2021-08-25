ui <- fluidPage(
  titlePanel("PCR Experiment Planning"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("primers",
                  "Number of Primers:",
                  min = 1,
                  max = 7,
                  value = 2),
      radioButtons("plate_format",
                   "Plate Format:",
                   c("96 Well" = "96_well", "384 Well" = "384_well")),
      radioButtons("exclude_border",
                   "Exclude Plate Border?",
                   c("Yes", "No"),
                   "Yes"),
      fileInput("rna_data",
                "Upload RNA Concentrations"),
      downloadButton("report",
                     "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sample Preparation", tableOutput("rna_table")),
        tabPanel("Mastermix Preparation"),
        tabPanel("Sample Layout", plotOutput("sample_layout"))
      )
    )
  )
)
# Need to validate file ext (allow for tsv, csv, excel files)
# See https://shiny.rstudio.com/reference/shiny/latest/fileInput.html

# Allow user to supply code to select which samples they want? Security risk?
# Maybe just allow for a simple click to select?
