## Give me a file
## Select a mapping file
## Run analysis
## Input for parameters

opts <- list(
  placeholder="Please upload a set of .gpr files"
)

source("common_functions.R")

shinyUI( fluidPage(
  
  includeCSS("styles.css"),
  includeScript("scripts.js"),
  
  ## OpenTip
  includeScript("www/opentip/opentip-jquery.min.js"),
  includeCSS("www/opentip/opentip.css"),
  
  HTML( "<style>#controls + div { height: 520px; }</style>" ),
  
  titlePanel(""),
  
  sidebarLayout(
  
    sidebarPanel(
      
      ## Names are of the form '<function>_<arg>'
      tabsetPanel(id = "controls",
        tabPanel("Read Data",
          fileInputWithHelp("gpr_files", "Upload GenePix .gpr Files", multiple=TRUE),
          fileInputWithHelp("mapping_file", "Upload a Mapping File", multiple=FALSE),
          uiOutput("file_status")
#           addHelp(selectizeInput("makePeptideSet_bgCorrect.method", "Method to be used for background correction", 
#             choices=c("normexp", "auto", "none", "subtract", "half", "min", "movingmin", "edwards")
#           )),
#           selectizeInput("makePeptideSet_rm.control.list", "Names of controls to be excluded", 
#             choices="", multiple=TRUE, options=opts),
#           selectizeInput("makePeptideSet_empty.control.list", "Names of empty controls",
#             choices="", multiple=TRUE, options=opts),
          
          # checkboxInput("makePeptideSet_log", "Perform a log2 transformation after BG correction?", TRUE),
#           checkboxInput("makePeptideSet_check_row_order", "Reduce slides to a common set of peptides?", TRUE),
#           actionButton("do_makePeptideSet", "Construct Peptide Set"),
#           uiOutput("do_makePeptideSet_status")
        ),
        tabPanel("Normalization",
          selectizeInputWithHelp("summarizePeptides_summary", "Summary Method", choices=list("median", "mean")),
          selectizeInputWithHelp("summarizePeptides_position", "Position Database",
            choices=list("pep_hxb2", "pep_mac239", "pep_hxb2JPT", "pep_m239smE543")
          ),
          numericInputWithHelp("slidingMean_width", "Sliding Mean Width", 5, 1, 100),
          checkboxInputWithHelp("slidingMean_split_by_space", "Split by Space?", TRUE),
          checkboxInputWithHelp("makePeptideSet_check_row_order", "Reduce slides to a common set of peptides?", TRUE),
          actionButton("do_summarizePeptides", "Normalize")
        ),
        tabPanel("Positivity Calls",
          selectInput("makeCalls_method", "Method", list(
            `Absolute`="absolute",
            `False Discovery Rate`="FDR"
          )),
          numericInput("makeCalls_cutoff", "Cutoff", 0.1, 0, 10, 0.01),
          selectInput("makeCalls_group", "Group", ""),
          HTML("<br />"),
          actionButton("do_makeCalls", "Make Calls"),
          downloadButton("export_calls"),
          uiOutput("makeCalls_status")
        )
      )
    ),
    
    mainPanel(
     tabsetPanel(
       tabPanel("Array Images",
         selectInput("arrayImageSelect", "Select an Array", ""),
         plotOutput("arrayImage", height=220),
         plotOutput("arrayResiduals", height=220)
       ),
       tabPanel("Calls",
         uiOutput("call_status"),
         dataTableOutput("calls")
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