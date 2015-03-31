library(grid)

textPlot <- function(text) {
  grid.newpage()
  grid.text(text)
}

parseTextField <- function(x) {
  if (!length(x)) return(NULL)
  if (identical(x, "")) return(NULL)
  return( scan( textConnection(x), what=character(), sep="," ) )
}

isPeptide <- function(x) {
  tmp <- gsub("[^[:alpha:]]", "", x, perl=TRUE)
  return( x == toupper(tmp) )
}

gpr_files_help <- "Select a set of GenePix results (.gpr) files you wish to analyze."

fileInputWithHelp <- function(inputId, label, multiple = FALSE, accept = NULL) {
  fi <- fileInput(inputId, label, multiple, accept)
  help_icon <- HTML(paste0('<i id=', inputId, '-help ',
                           'style="float: right; margin-top:18px;" ',
                           'class="fa fa-i fa fa-question-circle help-icon" ',
                           '></i>'))
  return(list(help_icon, fi))
}

selectizeInputWithHelp <- function(inputId, ..., options = NULL) {
  si <- selectizeInput(inputId, ..., options = options)
  help_icon <- HTML(paste0('<i id=', inputId, '-help ',
                           'style="float: right; margin-top:18px;" ',
                           'class="fa fa-i fa fa-question-circle help-icon" ',
                           '></i>'))
  return(list(help_icon, si))
}
#  icon <- list(
#    name="i",
#    attribs=list(
#      class="fa fa-question-circle help-icon",
#      style="float: right; margin-top: 3px;",
#      id=paste0(inputId, "-help")
#    ),
#    children=list()
#  )
#  class(icon) <- "shiny.tag"
#  output <- c(
#    list(icon),
#    si
#  )
#  class(output) <- c("shiny.tag.list", "list")
#  output
#}

numericInputWithHelp <- function(inputId, label, value, min = NA, max = NA, step = NA) {
  ni <- numericInput(inputId, label, value, min, max, step)
  help_icon <- HTML(paste0('<i id=', inputId, '-help ',
                           'style="float: right; margin-top:18px;" ',
                           'class="fa fa-i fa fa-question-circle help-icon" ',
                           '></i>'))
  return(list(help_icon, ni))
}
#  icon <- list(
#    name="i",
#    attribs=list(
#      class="fa fa-question-circle help-icon",
#      style="float: right; margin-top: 3px;",
#      id=paste0(inputId, "-help")
#    ),
#    children=list()
#  )
#  class(icon) <- "shiny.tag"
#  output <- c(
#    list(icon),
#    ni
#  )
#  class(output) <- c("shiny.tag.list", "list")
#  output
#}

checkboxInputWithHelp <- function(inputId, label, value = FALSE) {
  ci <- checkboxInput(inputId, label, value)
  help_icon <- HTML(paste0('<i id=', inputId, '-help ',
                           'style="float: right; margin-top:18px;" ',
                           'class="fa fa-i fa fa-question-circle help-icon" ',
                           '></i>'))
  return(list(help_icon, ci))
}
#  icon <- list(
#    name="i",
#    attribs=list(
#      class="fa fa-question-circle help-icon",
#      style="float: right; margin-top: 3px;",
#      id=paste0(inputId, "-help")
#    ),
#    children=list()
#  )
#  class(icon) <- "shiny.tag"
#  ni[[3]][[3]] <- icon
#  ni
#}
