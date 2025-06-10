#' @title process_data
#' @description Raw MS data processing using xcms.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param path Work directory.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ppm see xcms.
#' @param peakwidth See xcms.
#' @param snthresh See xcms.
#' @param prefilter See xcms.
#' @param fitgauss see xcms.
#' @param integrate see xcms.
#' @param mzdiff see xcms.
#' @param noise See xcms.
#' @param threads Number of threads.
#' @param binSize see xcms.
#' @param bw see xcms.
#' @param output_tic Output TIC plot or not.
#' @param output_bpc Output BPC plot or not.
#' @param output_rt_correction_plot Output rt correction plot or not.
#' @param min_fraction See xcms.
#' @param fill_peaks Fill peaks NA or not.
#' @param group_for_figure Which group you want to use to output
#' TIC and BPC and EIC. Default is QC.
#' @param detect_peak_algorithm Use the algorithm from xcms or local from massprocesser
#' @return Peak table.
#' @export

# #debug
# setwd(masstools::get_project_wd())
# setwd("vignettes/massprocesser_demo_data/POS2/")
# library(xcms)
# library(MSnbase)
# library(mzR)
# 
# path = "."
# polarity = "positive"
# ppm = 15
# peakwidth = c(5, 30)
# snthresh = 10
# prefilter = c(3, 500)
# fitgauss = FALSE
# integrate = 2
# mzdiff = 0.01
# noise = 500
# threads = 20
# binSize = 0.025
# bw = 5
# output_tic = FALSE
# output_bpc = FALSE
# output_rt_correction_plot = FALSE
# min_fraction = 0.5
# fill_peaks = FALSE
# output_peak_eic = TRUE
# is.table = "is.xlsx"
# group_for_figure = "QC"
# 
# massprocesser::process_data(path = ".",
#             polarity = "positive",
#             ppm = 15,
#             peakwidth = c(5, 30),
#             snthresh = 10,
#             prefilter = c(3, 500),
#             fitgauss = FALSE,
#             integrate = 2,
#             mzdiff = 0.01,
#             noise = 500,
#             threads = 20,
#             binSize = 0.025,
#             bw = 5,
#             output_tic = TRUE,
#             output_bpc = TRUE,
#             output_rt_correction_plot = TRUE,
#             min_fraction = 0.5,
#             fill_peaks = FALSE,
#             group_for_figure = "QC"
# )

process_data <-
  function(path = ".",
           polarity = c("positive", "negative"),
           ppm = 15,
           peakwidth = c(5, 30),
           snthresh = 10,
           prefilter = c(3, 500),
           fitgauss = FALSE,
           integrate = 2,
           mzdiff = 0.01,
           noise = 500,
           threads = 6,
           binSize = 0.025,
           bw = 5,
           output_tic = TRUE,
           output_bpc = TRUE,
           output_rt_correction_plot = TRUE,
           min_fraction = 0.5,
           fill_peaks = FALSE,
           group_for_figure = "QC",
           detect_peak_algorithm = c("xcms", "massprocesser")) {
    polarity <- match.arg(polarity)
    detect_peak_algorithm <-
      match.arg(detect_peak_algorithm)
    output_path <- file.path(path, "Result")
    dir.create(output_path, showWarnings = FALSE)
    intermediate_data_path <-
      file.path(output_path, "intermediate_data")
    dir.create(intermediate_data_path, showWarnings = FALSE)
    if (masstools::get_os() == "windows") {
      bpparam <-
        BiocParallel::SnowParam(workers = threads,
                                progressbar = FALSE)
    } else{
      bpparam <- BiocParallel::MulticoreParam(workers = threads,
                                              progressbar = FALSE)
    }    
    ##parameters
    massprocesser_parameters <- new(
      Class = "tidymass_parameter",
      pacakge_name = "massprocesser",
      function_name = "process_data",
      parameter = list(
        path = path,
        polarity = polarity,
        ppm = ppm,
        peakwidth = peakwidth,
        snthresh = snthresh,
        prefilter = prefilter,
        fitgauss = fitgauss,
        integrate = integrate,
        mzdiff = mzdiff,
        noise = noise,
        threads = threads,
        binSize = binSize,
        bw = bw,
        output_tic = output_tic,
        output_bpc = output_bpc,
        output_rt_correction_plot = output_rt_correction_plot,
        min_fraction = min_fraction,
        fill_peaks = fill_peaks
      ),
      time = Sys.time()
    )
    
    save(
      massprocesser_parameters,
      file = file.path(intermediate_data_path, "massprocesser_parameters")
    )
    ##----------------------------------------------------
    #peak detection
    
    f.in <- list.files(
      path = path,
      pattern = '\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)',
      recursive = TRUE,
      full.names = TRUE
    )
    
    sample_group <-
      unlist(lapply(stringr::str_split(string = f.in, pattern = "/"),
                    function(x) {
                      x[length(x) - 1]
                    }))
    
    sample_group[sample_group == ""] <- "Subject"
    
    sample_group[grep("\\.(mz[X]{0,1}ML|cdf)", sample_group)] <-
      "Group0"
    
    if (missing(group_for_figure)) {
      group_for_figure <- "QC"
    }
    
    if (!group_for_figure %in% sample_group) {
      group_for_figure2 <-
        plyr::count(sample_group) %>%
        dplyr::filter(freq == min(freq) &
                        stringr::str_to_lower(x) != "blank") %>%
        pull(x) %>%
        `[`(1) %>%
        as.character()
      message(
        crayon::yellow(
          group_for_figure,
          "is not in you directory, so set is as",
          group_for_figure2
        )
      )
      group_for_figure <- group_for_figure2
    }
    
    pd <-
      data.frame(
        sample_name = sub(
          basename(f.in),
          pattern = '\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)',
          replacement = "",
          fixed = TRUE
        ),
        sample_group = sample_group,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    
    ## Define colors for different groups
    temp_color <-
      c(
        "#E41A1C",
        "#377EB8",
        "#4DAF4A",
        "#984EA3",
        "#FF7F00",
        "#FFFF33",
        "#A65628",
        "#F781BF",
        "#999999"
      )
    
    group_colors <-
      paste0(temp_color[seq_along(unique(sample_group))],
             "60")
    
    names(group_colors) <- unique(sample_group)
    
    message(crayon::green("Reading raw data, it will take a while..."))
    
    if (any(dir(intermediate_data_path) == "raw_data")) {
      message(crayon::yellow("Use old saved data in Result."))
      load(file.path(intermediate_data_path, "raw_data"))
      ###change the files if changed the path
      raw_data <-
        change_file_path(object = raw_data, path = path)
    } else{
      raw_data <- MSnbase::readMSData(
        files = f.in,
        pdata = new("AnnotatedDataFrame", pd),
        mode = "onDisk",
        verbose = TRUE
      )
      save(raw_data,
           file = file.path(intermediate_data_path, "raw_data"))
    }
    
    message(crayon::red("OK"))
    
    #----------------------------------------------------------------------------
    message(crayon::green("Detecting peaks..."))
    
    ###peak detection
    cwp <-
      xcms::CentWaveParam(
        ppm = ppm,
        prefilter = prefilter,
        integrate = integrate,
        peakwidth = peakwidth,
        snthresh = snthresh,
        mzdiff = mzdiff,
        noise = noise,
        fitgauss = fitgauss,
      )
    
    if (any(dir(intermediate_data_path) == "xdata")) {
      message(crayon::yellow("Use old saved data in Result."))
      load(file.path(intermediate_data_path, "xdata"))
      xdata <-
        change_file_path(object = xdata, path = path)
    } else{
      if (detect_peak_algorithm == "xcms") {
        xdata <-
          tryCatch(
            xcms::findChromPeaks(
              object = raw_data,
              param = cwp,
              BPPARAM = bpparam
            ),
            error = function(e) {
              NULL
            }
          )
      } else{
        xdata <-
          tryCatch(
            findChromPeaks_sxt(
              object = raw_data,
              param = cwp,
              BPPARAM = bpparam
            ),
            error = function(e) {
              NULL
            }
          )
      }
      
      if (is.null(xdata)) {        
        if (detect_peak_algorithm == "xcms") {
          xdata <-
            try(xcms::findChromPeaks(raw_data,
                                     param = cwp,
                                     BPPARAM = bpparam),
                silent = FALSE)
        } else{
          xdata <-
            try(findChromPeaks_sxt(raw_data,
                                   param = cwp,
                                   BPPARAM = bpparam),
                silent = FALSE)
        }
        
        
      }
      
      if (is(xdata, class2 = "try-error")) {
        stop("Error in xcms::findChromPeaks.")
      }
      
      save(xdata,
           file = file.path(intermediate_data_path, "xdata"))
    }
    
    rm(list = "raw_data")
    gc()
    message(crayon::red("OK"))
    
    #-------------------------------------------------------
    #retention time correction
    #Alignment
    message(crayon::green("Correcting rentention time..."))
    
    if (any(dir(intermediate_data_path) == "xdata2")) {
      message(crayon::yellow("Use old saved data in Result."))
      load(file.path(intermediate_data_path, "xdata2"))
      xdata2 <-
        change_file_path(object = xdata2, path = path)
    } else{
      xdata2 <- try(xcms::adjustRtime(xdata,
                                      param = xcms::ObiwarpParam(binSize = 0.5)),
                    silent = FALSE)
      save(xdata2,
           file = file.path(intermediate_data_path, "xdata2"))
    }
    
    message(crayon::red("OK"))
    
    if (is(xdata2, class2 = "try-error")) {
      xdata2 <- xdata
      save(xdata2,
           file = file.path(intermediate_data_path, "xdata2"))
    } else{
      ## Plot also the difference of adjusted to raw retention time.
      if (output_rt_correction_plot) {
        message(crayon::green("Drawing RT correction plot..."))
        rt_correction_plot <-
          plot_adjusted_rt(object = xdata2,
                           group_for_figure = group_for_figure)
        
        ggplot2::ggsave(
          filename = file.path(output_path, "RT correction plot.pdf"),
          plot = rt_correction_plot,
          width = 20,
          height = 7
        )
        rm(list = c("rt_correction_plot"))
        message(crayon::red("OK"))
      }
    }
    
    rm(list = "xdata")
    
    ###TIC
    if (masstools::get_os() == "windows") {
      bpparam <-
        BiocParallel::SnowParam(workers = threads,
                                progressbar = TRUE)
    } else{
      bpparam <- BiocParallel::MulticoreParam(workers = threads,
                                              progressbar = TRUE)
    }
    
    if (output_tic) {
      message(crayon::green("Drawing TIC plot..."))
      
      figure_sample_name <-
        pd %>%
        dplyr::filter(sample_group %in% group_for_figure) %>%
        dplyr::pull(sample_name)
      
      idx <-
        which(basename(xdata2@processingData@files) %in% figure_sample_name)
      
      if(length(idx) > 15){
        idx <- sort(sample(idx, 15, replace = FALSE))
      }
      
      tic.plot <-
        tryCatch(
          xcms::chromatogram(
            object = filterFile(object = xdata2, file = idx),
            aggregationFun = "sum",
            BPPARAM = bpparam
          ),
          error = function(e)
            NULL
        )
      
      if (is.null(tic.plot)) {        
        tic.plot <-
          tryCatch(
            xcms::chromatogram(
              object = filterFile(object = xdata2, file = idx),
              aggregationFun = "sum",
              BPPARAM = bpparam
            ),
            error = function(e)
              NULL
          )
      }
      
      ## Plot all chromatograms.
      save(
        tic.plot,
        file = file.path(intermediate_data_path, "tic.plot"),
        compress = "xz"
      )
      
      plot <-
        plot_chromatogram(object = tic.plot,
                          title = "TIC",
                          group_for_figure = group_for_figure)
      
      if (!is.null(plot)) {
        ggplot2::ggsave(
          filename = file.path(output_path, "TIC.pdf"),
          plot = plot,
          width = 20,
          height = 7
        )
      }
      
      rm(list = c("plot", "tic.plot"))
      message(crayon::red("OK"))
    }
    
    ###BPC
    if (output_bpc) {
      message(crayon::green("Drawing BPC plot..."))
      
      figure_sample_name <-
        pd %>%
        dplyr::filter(sample_group %in% group_for_figure) %>%
        dplyr::pull(sample_name)
      
      idx <-
        which(basename(xdata2@processingData@files) %in% figure_sample_name)
      
      if(length(idx) > 15){
        idx <- sort(sample(idx, 15, replace = FALSE))
      }
      
      bpc.plot <-
        tryCatch(
          xcms::chromatogram(
            object = filterFile(object = xdata2, file = idx),
            aggregationFun = "max",
            BPPARAM = bpparam
          ),
          error = function(e)
            NULL
        )
      
      
      if (is.null(bpc.plot)) {        
        bpc.plot <-
          tryCatch(
            xcms::chromatogram(
              object = filterFile(object = xdata2, file = idx),
              aggregationFun = "max",
              BPPARAM = bpparam
            ),
            error = function(e)
              NULL
          )
      }
      
      ## Plot all chromatograms.
      save(
        bpc.plot,
        file = file.path(intermediate_data_path, "bpc.plot"),
        compress = "xz"
      )
      
      plot <- plot_chromatogram(object = bpc.plot,
                                title = "BPC",
                                group_for_figure = group_for_figure)
      if (!is.null(plot)) {
        ggplot2::ggsave(
          filename = file.path(output_path, "BPC.pdf"),
          plot = plot,
          width = 20,
          height = 7
        )
      }
      
      rm(list = c("plot", "bpc.plot"))
      message(crayon::red("OK"))
    }
    
    #-----------------------------------------------
    ## Perform the correspondence
    message(crayon::green("Grouping peaks across samples..."))
    
    if (any(dir(intermediate_data_path) == "xdata3")) {
      message(crayon::yellow("Use old saved data in Result."))
      load(file.path(intermediate_data_path, "xdata3"))
      xdata3 <-
        change_file_path(object = xdata3, path = path)
      
    } else{
      pdp <- xcms::PeakDensityParam(
        sampleGroups = xdata2$sample_group,
        minFraction = min_fraction,
        bw = bw,
        binSize = binSize,
        minSamples = 1,
        maxFeatures = 100
      )
      
      xdata3 <- xcms::groupChromPeaks(xdata2, param = pdp)
      save(xdata3,
           file = file.path(intermediate_data_path, "xdata3"))
    }
    
    message(crayon::red("OK"))
    rm(list = "xdata2")
    
    if (fill_peaks) {
      ## Filling missing peaks using default settings. Alternatively we could
      ## pass a FillChromPeaksParam object to the method.
      xdata3 <- xcms::fillChromPeaks(xdata3)
      save(xdata3,
           file = file.path(intermediate_data_path, "xdata3"))
    }
    
    message(crayon::green("Outputting peak table..."))
    
    ##output peak table
    values <- xcms::featureValues(xdata3, value = "into")
    definition <- xcms::featureDefinitions(object = xdata3)
    definition <- definition[, -ncol(definition)]
    definition <-
      definition[names(definition) != "peakidx"]
    definition <-
      definition@listData %>%
      do.call(cbind, .) %>%
      as.data.frame()
    
    peak_name <- xcms::groupnames(xdata3)
    peak_name <- paste(peak_name,
                       ifelse(polarity == "positive", "POS", "NEG"), sep = "_")
    
    colnames(values) <-
      stringr::str_replace(
        string = colnames(values),
        pattern = "\\.mz[X|x]{0,1}[M|m][L|l]",
        replacement = ""
      )
    
    peak_table <- data.frame(
      peak.name = peak_name,
      definition,
      values,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    rownames(peak_table) <- NULL
    
    colnames(peak_table) <-
      stringr::str_replace(
        string = colnames(peak_table),
        pattern = "\\.mz[X]{0,1}ML",
        replacement = ""
      )
    
    readr::write_csv(peak_table, file = file.path(output_path, "Peak_table.csv"))
    
    peak_table_for_cleaning <-
      definition %>%
      dplyr::select(-c(
        "mzmin",
        'mzmax',
        'rtmin',
        'rtmax',
        'npeaks',
        unique(sample_group)
      )) %>%
      dplyr::rename(mz = mzmed, rt = rtmed) %>%
      data.frame(
        variable_id = peak_name,
        .,
        values,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    
    readr::write_csv(peak_table_for_cleaning,
                     file = file.path(output_path, "Peak_table_for_cleaning.csv"))
    
    message(crayon::red("OK"))
    
    #####output mass_data class object
    sample_info <-
      pd %>%
      dplyr::rename(sample_id = sample_name,
                    group = sample_group) %>%
      dplyr::mutate(
        class = dplyr::case_when(
          stringr::str_detect(group, "QC") ~ "QC",
          stringr::str_detect(group, "Blank") ~ "Blank",
          stringr::str_detect(group, "Internal_standard") ~ "Internal_standard",
          TRUE ~ "Subject"
        )
      ) %>%
      dplyr::mutate(injection.order = seq_len(nrow(pd))) %>%
      dplyr::mutate(sample_id =
                      stringr::str_replace(sample_id, "\\.mz[M|m][L|l]", "")) %>%
      dplyr::mutate(sample_id = stringr::str_replace(sample_id, "\\.mz[X|x][M|m][L|l]", ""))
    
    expression_data <-
      peak_table_for_cleaning %>%
      dplyr::select(-c(variable_id:rt)) %>%
      as.data.frame()
    
    variable_info <-
      peak_table_for_cleaning %>%
      dplyr::select(variable_id:rt) %>%
      as.data.frame()
    
    rownames(expression_data) <- variable_info$variable_id
    rownames(variable_info) <- NULL
    
    object <-
      massdataset::create_mass_dataset(
        expression_data = expression_data,
        sample_info = sample_info,
        variable_info = variable_info
      )
    
    object@process_info$process_data <-
      massprocesser_parameters
    
    save(object, file = file.path(output_path, "object"))
    
    rm(list = c("peak_table", "peak_table_for_cleaning"))
    message(crayon::red("OK"))
    
    message(crayon::bgRed("All done!"))
  }
