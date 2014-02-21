shinyServer( function(input, output, session) {
  
  library(pepStat)
  library(PEP.db)
  data(pep_hxb2)
  source("common_functions.R")
  
  onButtonPress <- function(buttonId, x, env = parent.frame(), quoted = FALSE) {
    fun <- exprToFunction(x, env, quoted)
    observe({
      if (length(input[[buttonId]]) && input[[buttonId]] > 0) {
        isolate( fun() )
      } else {
        return( invisible(NULL) )
      }
    })
  }
  
  ## Globals
  pSet   <- NULL ## peptide set
  psSet  <- NULL ## summarized peptide set
  pnSet  <- NULL ## normalized peptide set
  psmSet <- NULL ## smoothed data set
  calls  <- NULL ## output from makeCalls(ps)
  
  gpr_folder   <- system.file("extdata/gpr_samples", package = "PEP.db")
  mapping.file <- system.file("extdata/mapping.csv", package = "PEP.db")
  
  ## Globals, which are updates when the appropriate buttons are pressed
  onButtonPress("do_makePeptideSet", {
    
    rm.control.list <- parseTextField(input$makePeptideSet_rm.control.list)
    empty.control.list <- parseTextField(input$makePeptideSet_empty.control.list)
    bgCorrect.method <- input$makePeptideSet_bgCorrect.method
    log <- input$makePeptideSet_log
    check.row.order <- input$makePeptideSet_check.row.order
    
    pSet <<- makePeptideSet(
      path=gpr_folder, 
      mapping.file=mapping.file,
      rm.control.list=rm.control.list,
      empty.control.list=empty.control.list,
      bgCorrect.method=bgCorrect.method,
      log=log,
      check.row.order=check.row.order
    )
    
    ## Normalize as well
    psSet <<- summarizePeptides(pSet, summary="mean", position=pep_hxb2)
    pnSet <<- normalizeArray(psSet)
    psmSet <<- slidingMean(pnSet, width=9)
    
    ## Update inputs
    updateSelectInput(session, "arrayImageSelect", "Select an Array",
      choices=sampleNames( phenoData(pSet) )
    )
    updateSelectInput(session, "arrayResidualsSelect", "Select an Array",
      choices=sampleNames( phenoData(pSet) )
    )
    
  })
  
  onButtonPress("do_makeCalls", {
    
    method <- input$makeCalls_method
    cutoff <- input$makeCalls_cutoff
    
    if (!is.null(psmSet)) {
      calls <<- makeCalls(psmSet, cutoff=cutoff, method=method)
    }
    
  })
  
  output$arrayImage <- renderPlot({
    input$do_makePeptideSet
    if (!is.null(pSet)) {
      if (!input$arrayImageSelect == "") {
        nms <- sampleNames( phenoData(pSet) )
        ind <- match(input$arrayImageSelect, nms)
        plotArrayImage(pSet, array.index=ind)
      }
    } else {
      
      textPlot("Please construct a pSet to view plots.")
    }
  })
  
  output$arrayResiduals <- renderPlot({
    input$do_makePeptideSet
    if (!is.null(pSet)) {
      if (!input$arrayResidualsSelect == "") {
        nms <- sampleNames( phenoData(pSet) )
        ind <- match(input$arrayResidualsSelect, nms)
        plotArrayResiduals(pSet, array.index=ind, smooth=TRUE)
      }
    } else {
      textPlot("Please construct a pSet to view plots.")
    }
  })
  
  output$debug <- renderPrint({
    
    R_send <- input$R_send
    if (input$R_send == 0) {
      return( invisible(NULL) )
    }
    
    isolate({
      code <- input$R_input
      result <- eval( parse( text=code ) )
      return(result)
    })
    
  })
  
})