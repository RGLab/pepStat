## Give me a file
## Select a mapping file
## Run analysis
## Input for parameters

shinyUI( pageWithSidebar(
  
  headerPanel("pepStat"),
  sidebarPanel(
    
#     selectInput("gpr_folder", "Select a Folder Containing GenePix Results",
#       choices=c(
#         system.file("extdata/gpr_samples", package = "PEP.db")
#       )
#     ),
#     
#     selectInput("mapping_file", "Select a Mapping File", 
#       choices=c(
#         system.file("extdata/mapping.csv", package = "PEP.db")
#       )
#     ),
    
    ## Names are of the form 'function_arg'
    h3("Peptide Set Construction"),
    textInput("makePeptideSet_rm.control.list", "Names of controls to be excluded", value=""),
    textInput("makePeptideSet_empty.control.list", "Names of empty controls", value=""),
    selectInput("makePeptideSet_bgCorrect.method", "Method to be used for background correction", 
      choices=c("normexp", "auto", "none", "subtract", "half", "min", "movingmin", "edwards")
    ),
    checkboxInput("makePeptideSet_log", "Perform a log2 transformation after BG correction?", TRUE),
    checkboxInput("makePeptideSet_check.row.order", "Reduce slides to a common set of peptides?", TRUE),
    actionButton("do_makePeptideSet", "Construct Peptide Set"),
    
    h3("Make Calls"),
    selectInput("makeCalls_method", "Method", list(
      `False Discovery Rate`="FDR",
      `Absolute`="absolute"
    )),
    numericInput("makeCalls_cutoff", "Cutoff", 0.1, 0, 10, 0.01),
    HTML("<br />"),
    actionButton("do_makeCalls", "Make Calls")
    
  ),
  
  mainPanel(
   tabsetPanel(
     tabPanel("Array Images",
       selectInput("arrayImageSelect", "Select an Array", ""),
       plotOutput("arrayImage")
     ),
     tabPanel("Array Residuals",
       selectInput("arrayResidualsSelect", "Select an Array", ""),
       plotOutput("arrayResiduals")
     ),
     tabPanel("Debug",
       tags$div(
         textInput("R_input", "Enter an R Command", ''),
         actionButton("R_send", label="Send Command")
       ),
       tags$div(
         verbatimTextOutput("debug")
       )
     )
   )
  )
  
) )