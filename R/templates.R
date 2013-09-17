mk_project <- function(path="./", project_name="pepStat_template"){
  #Create folder architecture
  message("Creating a new project in ", paste0(path, "/", project_name))
  invisible(sapply(paste0(path, "/", project_name,"/", c("reports","settings", paste0("data/", c("input", "output")))), dir.create, recursive=TRUE))
  #Add config.yaml, report.Rmd, README.txt
  file.copy(system.file("extdata/config.yaml", package="pepStat"), paste0(path, "/", project_name,"/settings"))
  file.copy(system.file("extdata/report.Rmd", package="pepStat"), paste0(path, "/", project_name,"/reports"))  
  file.copy(system.file("extdata/README.txt", package="pepStat"), paste0(path, "/", project_name,"/"))  
  invisible(1)
}

rm_project <- function(path="./", project_name="pepStat_template"){
  message("Removing the project ", paste0(path, "/", project_name))
  unlink(x=paste0(path, "/", project_name), recursive=TRUE)
}
