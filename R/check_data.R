#' @title check_targeted_table
#' @description Check targeted_table for extract_eic function.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param targeted_table targeted_table, 
#' 3 columns, 1. variable_id, 2. mz, 3. rt.
#' @return Notice of data checking.
#' @export
#' @examples 
#' targeted_table = 
#'   data.frame(variable_id = letters[1:5],
#'              mz = c(100, 100.01, 100.02, 100.03, 100.05),
#'              rt = c(1,2,3,4,5))
#' check_targeted_table(targeted_table)

check_targeted_table <- function(targeted_table) {
  if (missing(targeted_table)) {
    stop("No targeted_table prodvided.")
  }
  
  if (all(!class(targeted_table) %in% c("matrix", "data.frame", "tbl"))) {
    stop("targeted_table must be a data.frame/matrix.")
  }
  
  if (ncol(targeted_table) < 3) {
    stop("The first column of targeted_table must be variable_id, 
         mz, and rt (s).")
  }
  
  if (sum(colnames(targeted_table)[seq_len(3)] == 
          c("variable_id", "mz", "rt")) != 3) {
    stop("The first 3 columns must be variable_id, mz, and rt (s).")
  }
  
  if (sum(is.na(targeted_table[, seq_len(3)])) > 0) {
    stop("No NA is allowed in the first 3 columns.")
  }
  
  if (sum(targeted_table[, seq_len(3)] == "") > 0) {
    stop("No space is allowed in the first 3 columns.")
  }
  
  if (!is.character(targeted_table$variable_id)) {
    stop("variable_id must all be character.")
  }
  
  if (!is.numeric(targeted_table$mz)) {
    stop("mz must all be numeric.")
  }
  
  if (!is.numeric(targeted_table$rt)) {
    stop("rt must all be numeric.")
  }
}
