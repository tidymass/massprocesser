.onAttach <- function(...) {
  # needed <- core[!is_attached(core)]
  # if (length(needed) == 0)
  #   return()
  #
  # crayon::num_colors(TRUE)
  # massprocesser_attach()
  
  # if (!"package:conflicted" %in% search()) {
  #   x <- massprocesser_conflicts()
  #   msg(massprocesser_conflict_message(x), startup = TRUE)
  # }
  msg(paste0("Version ", massprocesser_version, " (", update_date, ')'))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}
