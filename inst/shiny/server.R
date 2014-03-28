## TODO
## Write calls matrix to file

library(data.table)

shinyServer( function(input, output, session) {
  
  library(pepStat)
  library(PEP.db)
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
  pSetSuccess <- FALSE
  
  ## Observer: updating the 'makePeptideSet_rm.control.list',
  ## 'makePeptideSet_empty.control.list' fields
  observe({
    gpr_files <- input$gpr_files
    
    ## If this is non-NULL, we just got some new files!
    if (is.null(gpr_files)) return(NULL)
    
    data <- lapply(gpr_files$datapath, function(x) {
      fread(x, skip=34)
    })
    
    ## The annotations are those not starting with a numeric value
    annotations <- Reduce(union, lapply(data, function(df) {
      unique(grep("^[^0-9]", df$Annotation, value=TRUE))
    }))
    
    updateSelectInput(session, "makePeptideSet_rm.control.list",
      "Names of controls to be excluded",
      choices=sort(annotations)
    )
    
    updateSelectInput(session, "makePeptideSet_empty.control.list",
      "Names of empty controls",
      choices=sort(annotations)
    )
    
  })
  
  ## Observer: start over
  observe({
    if (is.null(input$start_over)) return(NULL)
    pSet   <<- NULL ## peptide set
    psSet  <<- NULL ## summarized peptide set
    pnSet  <<- NULL ## normalized peptide set
    psmSet <<- NULL ## smoothed data set
    calls  <<- NULL ## output from makeCalls(ps)
    pSetSuccess <<- FALSE
  })
  
  onButtonPress("do_summarizePeptides", {
    
    gpr_files <- input$gpr_files
    mapping_file <- input$mapping_file
    summary <- input$summarizePeptides_summary
    position <- input$summarizePeptides_position
    
    if (is.null(gpr_files) || is.null(mapping_file)) {
      return(0)
    }
    
    ## Need to move around the files so that pepStat is happy. Shiny's defaults
    ## are odd, to say the least
    from <- gpr_files$datapath
    gpr_folder <- dirname(from)[1]
    to   <- file.path( dirname(from), gpr_files$name )
    file.copy(from, to)
    unlink(from)
    
    mapping.file <- mapping_file$datapath
    
    rm.control.list <- NULL #input$makePeptideSet_rm.control.list
    empty.control.list <- NULL #input$makePeptideSet_empty.control.list
    bgCorrect.method <- "normexp" #input$makePeptideSet_bgCorrect.method
    # log <- input$makePeptideSet_log ## apparently having this deselected breaks things
    check.row.order <- input$makePeptideSet_check_row_order
    
    pSet <<- makePeptideSet(
      path=gpr_folder, 
      mapping.file=mapping.file,
      rm.control.list=rm.control.list,
      empty.control.list=empty.control.list,
      bgCorrect.method=bgCorrect.method,
      # log=log,
      check.row.order=check.row.order
    )
    
    ## Update inputs
    updateSelectInput(session, "arrayImageSelect", "Select an Array",
      choices=sampleNames( phenoData(pSet) ),
      selected=sampleNames(phenoData(pSet))[1]
    )
    updateSelectInput(session, "arrayResidualsSelect", "Select an Array",
      choices=sampleNames( phenoData(pSet) ),
      selected=sampleNames(phenoData(pSet))[1]
    )
    
    groups <- names(pData(pSet))
    groups <- setdiff(groups, c("filename", "ptid", "visit"))
    updateSelectInput(session, "makeCalls_group", "Group",
      choices=c("None", groups)
    )
    
    ## Since 'position' is read in as a character string, we should get the
    ## actual data set it corresponds to
    call <- call("data", position)
    eval(call)
    
    psSet <<- summarizePeptides(pSet, summary=summary, position=get(position))
    pnSet <<- normalizeArray(psSet)
    psmSet <<- slidingMean(pnSet, width=9)
    
  })
  
  onButtonPress("do_makeCalls", {
    
    method <- input$makeCalls_method
    cutoff <- input$makeCalls_cutoff
    group <- input$makeCalls_group
    if (group == "None") group <- NULL
    
    if (!is.null(psmSet)) {
      calls <<- as.matrix(
        makeCalls(psmSet, cutoff=cutoff, method=method, group=group)
      )
    }
    
  })
  
  output$do_makePeptideSet_status <- renderUI({
    if (input$do_makePeptideSet == 0) return(invisible(NULL))
    print("do_makePeptideSet_status")
    if (is.null(psmSet)) {
      tagList(
        HTML("<br />"),
        p("Failed to construct and normalize peptide set.")
      )
    } else {
      tagList(
        HTML("<br />"),
        p("Peptide set constructed and normalized.")
      )
    }
  })
  
  output$arrayImage <- renderPlot({
    arrayImageSelect <- input$arrayImageSelect
    if (!is.null(pSet)) {
      if (arrayImageSelect != "") {
        nms <- sampleNames(phenoData(pSet))
        ind <- match(arrayImageSelect, nms)
        return(plotArrayImage(pSet, array.index=ind))
      }
    } else {
      textPlot("Please construct a peptide set to view plots.")
    }
  })
  
  output$arrayResiduals <- renderPlot({
    arrayImageSelect <- input$arrayImageSelect
    if (!is.null(pSet)) {
      if (input$arrayImageSelect != "") {
        nms <- sampleNames( phenoData(pSet) )
        ind <- match(arrayImageSelect, nms)
        return(plotArrayResiduals(pSet, array.index=ind, smooth=TRUE))
      }
    } else {
      invisible(NULL)
    }
  })
  
  output$call_status <- renderUI({
    input$do_makeCalls
    if (is.null(calls)) 
      tagList(
        p("No calls have been made yet."),
        HTML("<br />"),
        p("Please construct a peptide set from GenePix ",
        ".gpr files with the 'Read Data' tab, and then make calls using the interface ",
        "under the 'Make Calls' tab.")
      )
    else invisible(NULL)
  })
  
  output$calls <- renderDataTable(
    options=list(
      iDisplayLength = 10,
      iDisplayStart = 10
    ), {
    input$do_makeCalls
    if (!is.null(calls)) {
      output <- data.frame(Antigen=rownames(calls))
      output <- cbind(output, as.data.frame(calls))
      rownames(output) <- NULL
      return(output)
    } else {
      invisible(NULL)
    }
  })
  
  output$makeCalls_status <- renderUI({
    if (input$do_makeCalls == 0) return(invisible(NULL))
    if (is.null(calls)) {
      tagList(
        HTML("<br />"),
        p("Could not make calls! Have you constructed a peptide set yet?")
      )
    } else {
      tagList(
        HTML("<br />"),
        p("Calls have been constructed.")
      )
    }
  })
  
  output$export_calls <- downloadHandler(
    filename = function() {
      paste('calls-', Sys.Date(), ".csv", sep='')
    }, content = function(file) {
      if (is.null(calls)) {
        stop("No calls are available yet! Please navigate back to the Shiny application.")
      }
      write.table(calls, file,
        sep="\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE)
    }, contentType="text/csv"
  )
  
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
  
  output$file_status <- renderUI({
    gpr_files <- input$gpr_files
    mapping_file <- input$mapping_file
    
    gpr_text <- if (is.null(gpr_files)) {
      "No .gpr files have been uploaded."
    } else {
      paste(nrow(gpr_files), "files have been uploaded.")
    }
    
    mapping_text <- if (is.null(mapping_file)) {
      "No mapping file has been uploaded."
    } else {
      "A mapping file has been uploaded."
    }
    
    tagList(
      p(gpr_text),
      p(mapping_text)
    )
    
  })
  
})