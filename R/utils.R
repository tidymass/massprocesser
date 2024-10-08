




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
  
  if (isTRUE(theme$dark))
    crayon::white(x)
  else
    crayon::black(x)
  
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
  names <-
    vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "massprocesser")
  }
  names
}

invert <- function(x) {
  if (length(x) == 0)
    return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(paste0(...),
                crayon::make_style(grDevices::grey(level), grey = TRUE))
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
#' @param title Title of the plot.
#' @param interactive interactive or not.
#' @param group_for_figure What groups to show EIC.
#' @param sample_for_figure sample_for_figure
#' @return A ggplot2 object.
#' @export

plot_chromatogram <- function(object,
                              title.size = 13,
                              lab.size = 13,
                              axis.text.size = 12,
                              alpha = 0.5,
                              title = "",
                              interactive = FALSE,
                              group_for_figure = "QC",
                              sample_for_figure = NULL) {
  if (is.null(object)) {
    return(NULL)
  }
  options(warn = -1)
  info <- object@phenoData@data
  data <- object@.Data
  if (group_for_figure == "all") {
    group_for_figure <-
      unique(info$sample_group)
  }
  
  if (!is.null(sample_for_figure)) {
    sample_for_figure <-
      sample_for_figure[sample_for_figure %in% info$sample_name]
    if (length(sample_for_figure) == 0) {
      sample_for_figure <- NULL
    }
  }
  
  if (!is.null(sample_for_figure)) {
    data <-
      data[1, which(info$sample_name %in% sample_for_figure), drop = FALSE]
    info <-
      info[info$sample_name %in% sample_for_figure, , drop = FALSE]
  } else{
    data <-
      data[1, which(info$sample_group %in% group_for_figure), drop = FALSE]
    info <-
      info[info$sample_group %in% group_for_figure, , drop = FALSE]
  }
  
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
  
  data$sample <-
    data$sample %>%
    stringr::str_replace("\\.mzXML", "") %>%
    stringr::str_replace("\\.mzML", "") %>%
    stringr::str_replace("\\.mzxml", "") %>%
    stringr::str_replace("\\.mzml", "")
  
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
  
  if (length(unique(data$sample)) > 10) {
    plot <-
      plot +
      theme(legend.position = "none")
  }
  
  if (interactive) {
    if (requireNamespace("plotly", quietly = TRUE)) {
      plot <- plotly::ggplotly(plot)
    } else{
      message("Please install plotly package first.")
    }
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
#' @param interactive interactive or not.
#' @param group_for_figure What groups to show EIC.
#' @param sample_for_figure sample_for_figure
#' @return A ggplot2 object.
#' @export

plot_adjusted_rt <- function(object,
                             title.size = 13,
                             lab.size = 13,
                             axis.text.size = 12,
                             interactive = FALSE,
                             group_for_figure = "QC",
                             sample_for_figure = NULL) {
  diffRt <- xcms::rtime(object, adjusted = TRUE) -
    xcms::rtime(object,
                adjusted = FALSE)
  diffRt <- split(diffRt, MSnbase::fromFile(object))
  xRt <- xcms::rtime(object,
                     adjusted = TRUE,
                     bySample = TRUE)
  
  sample_name <- object@phenoData@data$sample_name
  sample_group <- object@phenoData@data$sample_group
  
  diffRt <-
    purrr::map(.x = seq_along(sample_name), function(idx) {
      data.frame(
        diff_rt = diffRt[[idx]],
        sample_name = sample_name[idx],
        sample_group = sample_group[idx]
      )
    }) %>%
    dplyr::bind_rows()
  
  # diffRt <- mapply(
  #   FUN = function(x, y) {
  #     list(data.frame(x, y, stringsAsFactors = FALSE))
  #   },
  #   x = diffRt,
  #   y = sample_name
  # )
  
  xRt <-
    purrr::map(.x = seq_along(sample_name), function(idx) {
      data.frame(rt = xRt[[idx]],
                 sample_name = sample_name[idx],
                 sample_group = sample_group[idx])
    }) %>%
    dplyr::bind_rows()
  
  # xRt <- mapply(
  #   FUN = function(x, y) {
  #     list(data.frame(x, y, stringsAsFactors = FALSE))
  #   },
  #   x = xRt,
  #   y = sample_name
  # )
  #
  # diffRt <- do.call(what = rbind, args = diffRt)
  # xRt <- do.call(rbind, xRt)
  
  temp.data =
    cbind(xRt, diffRt[, 1, drop = FALSE])
  
  # temp.data <-
  #   data.frame(xRt, diffRt, stringsAsFactors = FALSE)
  
  colnames(temp.data) <-
    c("rt", "sample.name", "sample_group", "diffRT")
  
  rm(list = c("object", "xRt", "diffRt"))
  
  if (group_for_figure == "all") {
    group_for_figure <-
      unique(sample_group)
  }
  
  if (!is.null(sample_for_figure)) {
    sample_for_figure <-
      sample_for_figure[sample_for_figure %in% sample_name]
    if (length(sample_for_figure) == 0) {
      sample_for_figure <- NULL
    }
  }
  
  if (!is.null(sample_for_figure)) {
    temp.data <-
      temp.data %>%
      dplyr::filter(sample.name %in% sample_for_figure)
  } else{
    temp.data <-
      temp.data %>%
      dplyr::filter(sample_group %in% group_for_figure)
  }
  
  
  temp.data$sample.name <-
    temp.data$sample.name %>%
    stringr::str_replace("\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)", "")
  
  plot <-
    ggplot2::ggplot(data = temp.data, ggplot2::aes(x = rt, y = diffRT)) +
    ggplot2::geom_line(data = temp.data, ggplot2::aes(color = sample.name)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention time (second)", y = "RT deviation (second)") +
    ggplot2::theme(
      # legend.position = "none",
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
  
  if (length(unique(temp.data$sample.name)) > 10) {
    plot <-
      plot +
      theme(legend.position = "none")
  }
  
  if (interactive) {
    if (requireNamespace("plotly", quietly = TRUE)) {
      plot <- plotly::ggplotly(plot)
    } else{
      message("Please install plotly package first.")
    }
    
  }
  
  plot
}
