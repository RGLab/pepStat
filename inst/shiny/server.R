library(data.table)
library(pepStat)
library(Pviz)
library(pepDat)
source("common_functions.R")

shinyServer( function(input, output, session) {

  onClick <- function(buttonId, x, env = parent.frame(), quoted = FALSE) {
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
  gpr_files_ready <- FALSE
  mapping_file_ready <- FALSE

  pos_db <- NULL ## the position database that gets loaded
  pSet   <- NULL ## peptide set
  psSet  <- NULL ## summarized peptide set
  pnSet  <- NULL ## normalized peptide set
  psmSet <- NULL ## smoothed data set
  calls  <- NULL ## output from makeCalls(ps)
  pSetSuccess <- FALSE
  restab <- NULL ## results table
  sbc    <- NULL ## did we split by clade?
  clades <- NULL ## what are the clades?

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

  ## Observer: make sure the GPR files are okay (TODO)
  observe({

    gpr_files <- input$gpr_files
    if (is.null(gpr_files)) {
      output$gpr_files_status <- renderUI({
        p("Please upload one or more GenePix .gpr files.")
      })
      return(NULL)
    }

    output$gpr_files_status <- renderUI({
      p( nrow(gpr_files), "GenePix .gpr files", if (nrow(gpr_files) == 1) "has" else "have",
        "been uploaded.")
    })

    gpr_files_ready <<- TRUE

  })

  ## Observer: make sure the mapping file is of correct format
  observe({

    mapping.file <- input$mapping_file
    if (is.null(mapping.file)) {
      output$mapping_file_status <- renderUI({
        p("Please upload a mapping file.")
      })
      return(NULL)
    }


    mp <- tryCatch(read.csv(mapping.file$datapath, header = TRUE), error=function(e) {
      output$mapping_file_status <- renderUI({
        p(style="color: red;", "Error: could not read mapping file!")
      })
      return(NULL)
    })

    if (!all( c("filename", "ptid", "visit") %in% names(mp))) {
      output$mapping_file_status <- renderUI({
        p(style="color: red;", "Error: the mapping file must include columns",
        "'filename', 'ptid', and 'visit'.")
      })
      return(NULL)
    }

    output$mapping_file_status <- renderUI({
      p("Mapping file successfully uploaded.")
    })

    mapping_file_ready <<- TRUE

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

  onClick("do_makePeptideSet", {

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

    tryCatch({
      pSet <<- makePeptideSet(
        path=gpr_folder,
        mapping.file=mapping.file,
        rm.control.list=rm.control.list,
        empty.control.list=empty.control.list,
        bgCorrect.method=bgCorrect.method,
        # log=log,
        check.row.order=check.row.order
      )
    }, error=function(e) {
      output$gpr_files_status <- renderUI({
        p(style="color: red;", "INTERNAL ERROR: Could not construct PeptideSet!")
      })
      output$mapping_file_status <- renderUI({
        tagList(
          pre( capture.output(e) )
        )
      })
      return(NULL)
    })

    if (is.null(pSet)) {
      return(NULL)
    }

    ## Update inputs
    updateSelectInput(session, "arrayImageSelect", "Select an Array",
      choices=sampleNames( phenoData(pSet) ),
      selected=sampleNames(phenoData(pSet))[1]
    )
    updateSelectInput(session, "arrayResidualsSelect", "Select an Array",
      choices=sampleNames( phenoData(pSet) ),
      selected=sampleNames(phenoData(pSet))[1]
    )

    output$pSet_status <- renderUI({
      p("PeptideSet successfully constructed.")
    })

  })

  onClick("do_summarizePeptides", {

    summary <- input$summarizePeptides_summary
    position <- input$summarizePeptides_position

    width <- input$slidingMean_width
    split.by.clade <- input$slidingMean_split_by_clade

    groups <- names(pData(pSet))
    groups <- setdiff(groups, c("filename", "ptid", "visit"))
    updateSelectInput(session, "makeCalls_group", "Group",
      choices=c("None", groups)
    )

    ## Since 'position' is read in as a character string, we should get the
    ## actual data set it corresponds to
    call <- call("data", position)
    eval(call)

    pos_db <<- get(position)
    psSet <<- summarizePeptides(pSet, summary=summary, position=pos_db)
    pnSet <<- normalizeArray(psSet)
    psmSet <<- slidingMean(pnSet, width=width, split.by.clade=split.by.clade)

    output$summarize_status <- renderUI({
      p("Peptide set successfully normalized.")
    })

  })

  onClick("do_makeCalls", {

    method <- input$makeCalls_method
    cutoff <- input$makeCalls_cutoff
    group <- input$makeCalls_group
    if (group == "None") group <- NULL

    if (!is.null(psmSet)) {
      calls <<- as.matrix(
        makeCalls(psmSet, cutoff=cutoff, method=method, group=group)
      )
      restab <<- restab(psmSet, calls)
      clades <<- sort(unique(unlist(strsplit(restab$clade, ",", fixed=TRUE))))

      ## Depending on whether we split by clade or not, we want to display
      ## different plots
      sbc <<- preproc(psmSet)$split.by.clade

      ## The 'split' case --> plot_clade
      if (!is.null(sbc) && isTRUE(sbc)) {

        ## Clade selection UI
        output$clades <- renderUI({
          selectizeInput("clades", "Clades", choices=clades, selected=NULL, multiple=TRUE)
        })

      } else {

        ## cleanup
        output$clades <- renderUI("")

      }

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
        p("No calls have been made yet.")
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
      output <- restab(psmSet, calls)
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
        p("Could not make calls! Have you constructed and normalized your peptide set yet?")
      )
    } else {
      tagList(
        HTML("<br />"),
        p("Calls have been made.")
      )
    }
  })

  output$download <- downloadHandler(

    ## restab
    ## exprs(psmSet)

    filename = function() {
      paste('pepStat-analysis-', Sys.Date(), ".zip", sep='')
    },

    content = function(file) {

      if (is.null(calls)) {
        stop("No calls are available yet! Please navigate back to the Shiny application.")
      }

      date <- Sys.Date()
      dir <- file.path( tempdir(), "pepStat" )
      if (!file.exists(dir) && !dir.create(dir)) {
        stop("Couldn't create an output directory in your temporary directory.\nPlease check your permissions and confirm that you have access to the directory pointed at:\n", dir, ".")
      }

      owd <- getwd()
      on.exit(setwd(owd))
      setwd(dir)

      ## write out restab
      restab_file <- paste0("results-", date, ".txt")

      write.table(restab, file=restab_file,
        sep="\t",
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE
      )

      ## write out expression matrix
      exprs_file <- paste0("exprs-", date, ".txt")
      exprs <- as.data.frame(exprs(psmSet))
      exprs <- cbind(
        peptide=gsub("_.*", "", rownames(exprs)),
        clade=gsub(".*_", "", rownames(exprs)),
        exprs,
        stringsAsFactors=FALSE
      )
      write.table(exprs, file=exprs_file,
        sep="\t",
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE
      )

      zipfile <- "results.zip"
      zip(zipfile, c(restab_file, exprs_file))
      file.rename("results.zip", file)

    },

    contentType = "application/zip"

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

  output$makePeptideSet <- renderUI({
    gpr_files <- input$gpr_files
    mapping_file <- input$mapping_file

    if (gpr_files_ready && mapping_file_ready) {
      actionButton("do_makePeptideSet", "Construct PeptideSet")
    }

  })

  output$Pviz_plot <- renderPlot({

    dmc <- input$do_makeCalls
    clades <- input$clades
    from <- input$Pviz_from
    to <- input$Pviz_to

    if (is.null(restab)) {
      grid.text("Please make calls before visualizing tracks")
      return(invisible(NULL))
    }

    if (from == 0 && to == 0) to <- max(restab$position)

    if (!isTRUE(sbc)) {
      Pviz::plot_inter(restab, sequence=metadata(pos_db)$sequence, from=from, to=to)
    } else {
      Pviz::plot_clade(restab, clades, sequence=metadata(pos_db)$sequence, from=from, to=to)
    }

  })

  onClick("reset", {
    updateNumericInput(session, "Pviz_from", "From", 0)
    updateNumericInput(session, "Pviz_to", "To", 0)
    if (!is.null(clades)) {
      updateSelectInput(session, "clades", "Clades", choices=clades, selected=NULL)
    }
  })

})
