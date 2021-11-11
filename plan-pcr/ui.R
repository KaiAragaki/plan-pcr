library(sass)
library(gt)
css <- sass(sass_file("./www/style.sass"))
ui <- fluidPage(

  tags$head(tags$style(css)),
  shinyFeedback::useShinyFeedback(),
  titlePanel("PCR Experiment Planning"),
  
  sidebarLayout(
    sidebarPanel(
      p("File should have at least 1 column - RNA concentrations."),
      fileInput("rna_data",
                "Upload RNA Concentrations",
                accept = c(".csv", ".tsv", ".xls", ".xlsx"),
                width = 250),
      p("You should likely have at least two primers: One for your gene of interest, and one of your endogenous control, eg GAPDH."),
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
      p("Include borders if you need the extra space."),
      radioButtons("exclude_border",
                   "Exclude Plate Border?",
                   c("Yes" = TRUE, "No" = FALSE), 
                   selected = TRUE),
      textInput("primer_names", 
                "Primer Names (optional)", 
                placeholder = "Primer 1; Primer 2...",
                width = 250),
      p("The report will include all plots and tables listed here, plus the options set above."),
      downloadButton("get_report",
                     "Download Report")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sample Preparation", gt_output("sample_prep"), p("This table shows how you should dilute samples to get them at the proper concentration for doing PCR. First, dilute by the dilution factor (a dilution factor of 1 means no dilution). Then, mix the indicated volume of diluted RNA (column 5) with the indicated amount of H2O (column 6).")),
        tabPanel("Mastermix Preparation", gt_output("mm_prep"), p("This table shows how much master mix you should create *for each primer*. Combine the following reagents as well as *one* primer for each mastermix.")),
        tabPanel("Mastermix Layout", p("This plot shows you where to put each mastermix on the plate."), plotOutput("mm_layout")),
        tabPanel("Sample Layout", p("This plot shows you were to put each sample on the plate."), plotOutput("sample_layout"))
      )
    )
  )
)
