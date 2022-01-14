
msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("massprocesser.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
  
}

#' List all packages in the massprocesser
#'
#' @param include_self Include massprocesser in the list?
#' @return massprocesser_packages
#' @export
#' @examples
#' massprocesser_packages()

massprocesser_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("massprocesser")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <- vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "massprocesser")
  }
  names
}

invert <- function(x) {
  if (length(x) == 0) return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(
    paste0(...),
    crayon::make_style(grDevices::grey(level), grey = TRUE)
  )
}


#' @title plot_chromatogram
#' @description Draw TIC or BPC.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object Object for tic.plot or bpc.plot.
#' @param title.size Font size of title..
#' @param lab.size Font size of lab title.
#' @param axis.text.size Font size of axis text.
#' @param alpha alpha.
#' @param title Title of the plot..
#' @param interactive interactive.
#' @param group_for_figure What groups to show EIC.
#' @return A ggplot2 object.
#' @export

plot_chromatogram <- function(
  object,
  title.size = 15,
  lab.size = 15,
  axis.text.size = 15,
  alpha = 0.5,
  title = "",
  interactive = FALSE,
  group_for_figure = "QC"
){
  options(warn = -1)
  info <- object@phenoData@data
  data <- object@.Data
  data <-
    data[1, which(info$sample_group %in% group_for_figure), drop = FALSE]
  info <-
    info[info$sample_group %in% group_for_figure, , drop = FALSE]
  
  if (nrow(info) == 0) {
    return(NULL)
  }
  
  if (nrow(info) > 10) {
    idx <- sort(sample(seq_len(nrow(info)), 10))
    info <- info[idx, , drop = FALSE]
    data <- data[, idx, drop = FALSE]
  }
  
  rm(list = c("object"))
  data <- apply(data, 2, function(x) {
    x <- x[[1]]
    x <-
      data.frame(
        "mz" = x@rtime,
        "intensity" = x@intensity,
        stringsAsFactors = FALSE
      )
    list(x)
  })
  
  data <- lapply(data, function(x) {
    x[[1]]
  })
  
  data <- mapply(
    FUN = function(x, y, z) {
      x <- data.frame(
        x,
        "group" = y,
        "sample" = z,
        stringsAsFactors = FALSE
      )
      list(x)
    },
    x = data,
    y = info[, 2],
    z = info[, 1]
  )
  
  # data <- lapply(data, function(x){
  #   x <- plyr::dlply(.data = x, .variables = plyr::.(sample))
  # })
  
  data <- do.call(rbind, args = data)

  plot <-
    ggplot2::ggplot(data = data,
                    ggplot2::aes(x = mz, y = intensity)) +
    ggplot2::geom_line(data = data,
                       mapping = ggplot2::aes(colour = sample, 
                                              group = sample)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention time (RT, second)", 
                  y = "Intensity", title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black",
        size = title.size,
        face = "plain",
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        color = "black",
        size = lab.size,
        face = "plain"
      ),
      axis.text = ggplot2::element_text(
        color = "black",
        size = axis.text.size,
        face = "plain"
      )
    )
  
  if (interactive) {
    plot <- plotly::ggplotly(plot)
  }
  
  return(plot)
  
}


#' @title plot_adjusted_rt
#' @description plot_adjusted_rt
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object object
#' @param title.size title.size
#' @param lab.size lab.size
#' @param axis.text.size axis.text.size
#' @return A ggplot2 object.

plot_adjusted_rt <- function(
  object,
  title.size = 15,
  lab.size = 15,
  axis.text.size = 15
){
  diffRt <- xcms::rtime(object, adjusted = TRUE) - xcms::rtime(object,
                                                               adjusted = FALSE)
  diffRt <- split(diffRt, MSnbase::fromFile(object))
  xRt <- xcms::rtime(object,
                     adjusted = TRUE,
                     bySample = TRUE)
  
  sample_name <- object@phenoData@data$sample_name
  sample_group <- object@phenoData@data$sample_group
  
  diffRt <- mapply(
    FUN = function(x, y) {
      list(data.frame(x, y, stringsAsFactors = FALSE))
    },
    x = diffRt,
    y = sample_name
  )
  
  xRt <- mapply(
    FUN = function(x, y) {
      list(data.frame(x, y, stringsAsFactors = FALSE))
    },
    x = xRt,
    y = sample_name
  )
  
  diffRt <- do.call(what = rbind, args = diffRt)
  xRt <- do.call(rbind, xRt)
  
  temp.data <-
    data.frame(xRt, diffRt, stringsAsFactors = FALSE)
  
  colnames(temp.data) <-
    c("rt", "sample.name", "diffRT", "sample.name2")
  rm(list = c("object", "xRt", "diffRt"))
  
  plot <-
    ggplot2::ggplot(data = temp.data, ggplot2::aes(x = rt, y = diffRT)) +
    ggplot2::geom_line(data = temp.data, ggplot2::aes(color = sample.name)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention time (second)", y = "RT deviation (second)") +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(
        color = "black",
        size = lab.size,
        face = "plain"
      ),
      axis.text = ggplot2::element_text(
        color = "black",
        size = axis.text.size,
        face = "plain"
      )
    )
  plot
}

