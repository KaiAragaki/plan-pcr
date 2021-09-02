library(sass)
library(gt)
css <- sass(sass_file("./www/style.sass"))
ui <- fluidPage(

  tags$head(tags$style(css)),
  shinyFeedback::useShinyFeedback(),
  titlePanel("PCR Experiment Planning"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("rna_data",
                "Upload RNA Concentrations",
                accept = c(".csv", ".tsv", ".xls", ".xlsx"),
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
      textInput("primer_names", 
                "Primer Names (optional)", 
                placeholder = "Primer 1; Primer 2...",
                width = 250),
      downloadButton("report",
                     "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sample Preparation", p(), gt_output("sample_prep")),
        tabPanel("Mastermix Preparation", p(), gt_output("mm_prep")),
        tabPanel("Mastermix Layout", plotOutput("mm_layout")),
        tabPanel("Sample Layout", plotOutput("sample_layout"))
      )
    )
  )
)
