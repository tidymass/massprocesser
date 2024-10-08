.onAttach <- function(...) {
  # needed <- core[!is_attached(core)]
  # if (length(needed) == 0)
  #   return()
  #
  crayon::num_colors(TRUE)
  massprocesser_attach()
  
  # if (!"package:conflicted" %in% search()) {
  #   x <- massprocesser_conflicts()
  #   msg(massprocesser_conflict_message(x), startup = TRUE)
  # }
  packageStartupMessage(paste0(
    "massprocesser ",
    massprocesser_version,
    " (",
    update_date,
    ')'
  ))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}


change_file_path <- function(object,
                             path = ".") {
  files <-
    object@processingData@files
  
  new_files <-
    list.files(
      path = path,
      pattern = '\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)',
      recursive = TRUE,
      full.names = TRUE,
      include.dirs = TRUE
    ) %>%
    normalizePath()
  
  if (length(new_files) != length(files)) {
    stop("Old save xdata is different with your mzXML data now,
             delete the xdata and run again.")
  }
  
  files_name <- basename(files)
  new_files_name <- basename(new_files)
  
  if (sum(files_name == new_files_name) != length(files_name)) {
    stop("Old save xdata is different with your mzXML data now,
             delete the xdata and run again.")
  }
  
  idx <-
    match(new_files_name, files_name)
  new_files <- new_files[idx]
  
  if (sum(files == new_files) != length(files)) {
    object@processingData@files <-
      new_files
  }
  return(object)
}

# filter_XCMSnExp <-
#   function(object, index){
#     object@processingData@files <-
#       object@processingData@files[index]
#     return(object)
#   }


MSnExpReqFvarLabels <- function() {
  # labels <- MSnbase:::.MSnExpReqFvarLabels
  c(
    "fileIdx",
    "spIdx",
    "acquisitionNum",
    "retentionTime",
    "msLevel",
    "precursorScanNum"
  )
}