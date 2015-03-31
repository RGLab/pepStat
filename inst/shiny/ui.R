library(shiny)

## Protect against masking from other packages defining tags
tags <- shiny::tags

opts <- list(
  placeholder="Please upload a set of .gpr files"
)

source("common_functions.R")

shinyUI( fluidPage(

  ## OpenTip
  includeScript("www/opentip/opentip-jquery.min.js"),
  includeCSS("www/opentip/opentip.css"),

  ## Own
  includeCSS("styles.css"),
  includeScript("scripts.js"),
  includeScript("tooltips.js"),

  HTML( "<style>#controls + div { height: 520px; }</style>" ),

  titlePanel(""),

  sidebarLayout(

    sidebarPanel(

      ## Names are of the form '<function>_<arg>'
      tabsetPanel(id = "controls",
        tabPanel("Read Data",
          fileInputWithHelp("gpr_files", "Upload GenePix .gpr Files", multiple=TRUE),
          fileInputWithHelp("mapping_file", "Upload a Mapping File", multiple=FALSE),

          h3("File Status"),
          uiOutput("gpr_files_status"),
          uiOutput("mapping_file_status"),
          uiOutput("makePeptideSet"),
          HTML("<br/>"),
          uiOutput("pSet_status")
        ),
        tabPanel("Normalization",
                 selectizeInputWithHelp("summarizePeptides_summary", "Summary Method", choices=list("median", "mean")),
                 selectizeInputWithHelp("summarizePeptides_position", "Position Database",
                                        choices=list("pep_hxb2", "pep_mac239", "pep_hxb2JPT", "pep_m239smE543")
                 ),
                 numericInputWithHelp("slidingMean_width", "Sliding Mean Width", 9, 1, 100),
                 checkboxInputWithHelp("slidingMean_split_by_clade", "Split by Clade?", TRUE),
                 checkboxInputWithHelp("makePeptideSet_check_row_order", "Reduce slides to a common peptides set?",
                                       TRUE),
                 actionButton("do_summarizePeptides", "Normalize"),
                 HTML("<br />"),
                 uiOutput("summarize_status")
        ),
        tabPanel("Positivity Calls",
                 selectizeInputWithHelp("makeCalls_method", "Method", list(
                   `Absolute`="absolute",
                   `False Discovery Rate`="FDR"
                 )),
                 numericInputWithHelp("makeCalls_cutoff", "Cutoff", 1.2, 0, 10, 0.01),
                 HTML("<br />"),
                 selectizeInputWithHelp("makeCalls_group", "Group", ""),
                 HTML("<br />"),
                 actionButton("do_makeCalls", "Make Calls"),
                 downloadButton("download"),
                 uiOutput("makeCalls_status")
        )
      )
    ),

    mainPanel(
      tabsetPanel(id = "main_panel",
                  tabPanel("Array Images",
                           HTML("<br>"),
                           selectInput("arrayImageSelect", "Select an Array", ""),
                           plotOutput("arrayImage", height=250),
                           plotOutput("arrayResiduals", height=250)
                  ),
                  tabPanel("Calls",
                           HTML("<br>"),
                           uiOutput("call_status"),
                           dataTableOutput("calls")
                  ),
                  tabPanel("Visualization",
                           HTML("<br>"),
                           uiOutput("clades"),
                           uiOutput("rangeSlider"),
                           tags$div(style="overflow: auto;",
                                    #tags$div(style="float: left;", numericInput("Pviz_from", "From", 0)),
                                    #tags$div(style="float: left; margin-left: 20px;", numericInput("Pviz_to", "To", 0)),
                                    tags$div(style="float: right; margin-left: 20px;", actionButton("reset", "Reset"))
                           ),
                           plotOutput("Pviz_plot")
                  ),
                  tabPanel("Debug",
                           tags$div(
                             textInput("R_input", "Enter an R Command", ''),
                             tags$button(id = "R_send", type = "button", class = "btn action-button", "Submit",
                                         style = "margin-bottom: 10px;")
                           ),
                           tags$div(
                             verbatimTextOutput("debug")
                           )
                  )
      )
    )

  ) ) )
