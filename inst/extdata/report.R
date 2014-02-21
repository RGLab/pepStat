mk_project <- function(path="./", project_name="project", report_name="report"){
  #Create folder architecture
  message("Creating a new project in ", paste0(path, "/", project_name))
  invisible(sapply(paste0(path, "/", project_name,"/", c("reports","settings", paste0("data/", c("input", "output")))), dir.create, recursive=TRUE))
  #Add config.yaml, report.Rmd, README.txt
  file.copy(system.file("extdata/config.yaml", package="pepStat"), paste0(path, "/", project_name,"/settings/", report_name, ".yaml"))
  file.copy(system.file("extdata/report.Rmd", package="pepStat"), paste0(path, "/", project_name,"/reports/", report_name, ".Rmd"))  
  file.copy(system.file("extdata/README.txt", package="pepStat"), paste0(path, "/", project_name,"/"))  
  #TEMP dev
  file.copy(list.files(system.file("extdata/RVV/", package="pepStat"), full.names=TRUE), paste0(path, "/", project_name,"/data/input/"))
  file.copy(system.file("extdata/mapping.csv", package="pepStat"), paste0(path, "/", project_name,"/data/input/"))
  
  invisible(1)
}

run_report <- function(path="./", project_name="project", report_name="report"){
  if(!all(c("data", "reports", "settings", "README.txt") %in% list.files(paste0(path, "/", project_name))))
    warning("The given project does not seem to be complete.")
  #Need absolute path in projdir, otherwise, the report gets executed where it is located (./project/reports/).
  #However getwd + path wont work if the given path is already absolute
  projdir <- paste(file_path_as_absolute(path), project_name, sep="/")
  input <- paste0(projdir, "/reports/", report_name, ".Rmd")
  output <- paste0(projdir, "/data/output/", report_name, ".html")
  knit2html(input = input, output=output, quiet=FALSE)
}

rm_project <- function(path="./", project_name="project"){
  if(all(c("data", "reports", "settings", "README.txt") %in% list.files(paste0(path, "/", project_name)))){
    message("Removing the project ", paste0(path, "/", project_name))
    unlink(x=paste0(path, "/", project_name), recursive=TRUE)
  } else{
    stop("The given path and project_name do not look like a valid project.")
  }
}

clear_all_caches <- function(path="./", project_name="project"){
  to_remove <- list.files(paste0(path, "/", project_name, "/reports/"), pattern="^cache_", full.names=TRUE)
  unlink(x=to_remove, recursive=TRUE)
}
