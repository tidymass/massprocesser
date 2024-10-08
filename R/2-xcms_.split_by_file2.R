###remove dot
split_by_file2_sxt <- function(x,
                                msLevel. = unique(msLevel(x)),
                                subsetFeatureData = FALSE,
                                to_class = "OnDiskMSnExp",
                                keep_sample_idx = FALSE) {
  if (is(x, "XCMSnExp") && hasAdjustedRtime(x))
    x@featureData$retentionTime <- adjustedRtime(x)
  if (subsetFeatureData) {
    fcs <- intersect(
      c(
        # MSnbase:::.MSnExpReqFvarLabels,
        MSnExpReqFvarLabels(),
        "centroided",
        "polarity",
        "seqNum"
      ),
      colnames(.fdata(x))
    )
    x <- selectFeatureData(x, fcol = fcs)
  }
  fdl <- split.data.frame(x@featureData, as.factor(fromFile(x)))
  procd <- x@processingData
  expd <- new(
    "MIAPE",
    instrumentManufacturer = x@experimentData@instrumentManufacturer[1],
    instrumentModel = x@experimentData@instrumentModel[1],
    ionSource = x@experimentData@ionSource[1],
    analyser = x@experimentData@analyser[1],
    detectorType = x@experimentData@detectorType[1]
  )
  create_object <- function(i, fd, x, to_class) {
    a <- new(to_class)
    slot(procd, "files", check = FALSE) <- x@processingData@files[i]
    slot(a, "processingData", check = FALSE) <- procd
    slot(a, "featureData", check = FALSE) <-
      extractROWS(fd, which(fd$msLevel %in% msLevel.))
    if (!nrow(a@featureData))
      warning(
        "No MS level ",
        msLevel.,
        " spectra present in file ",
        basename(x@processingData@files[i]),
        call. = FALSE
      )
    else
      a@featureData$fileIdx <- 1L
    slot(a, "experimentData", check = FALSE) <- expd
    slot(a, "spectraProcessingQueue", check = FALSE) <-
      x@spectraProcessingQueue
    slot(a, "phenoData", check = FALSE) <-
      x@phenoData[i, , drop = FALSE]
    a
  }
  if (to_class == "XCMSnExp" &&
      is(x, "XCMSnExp") && hasChromPeaks(x)) {
    if (any(colnames(.chrom_peak_data(x@msFeatureData)) == "ms_level")) {
      pk_idx <- which(.chrom_peak_data(x@msFeatureData)$ms_level %in% msLevel.)
    } else
      pk_idx <- seq_len(nrow(chromPeaks(x@msFeatureData)))
    fct <- as.factor(as.integer(chromPeaks(x@msFeatureData)[pk_idx, "sample"]))
    pksl <- split.data.frame(chromPeaks(x@msFeatureData)[pk_idx, , drop = FALSE], fct)
    pkdl <- split.data.frame(extractROWS(.chrom_peak_data(x@msFeatureData), pk_idx), fct)
    res <- vector("list", length(fileNames(x)))
    for (i in seq_along(res)) {
      a <- create_object(i, fdl[[i]], x, to_class)
      newFd <- new("MsFeatureData")
      pks <- pksl[[as.character(i)]]
      if (!is.null(pks) && nrow(pks)) {
        if (!keep_sample_idx)
          pks[, "sample"] <- 1
        chromPeaks(newFd) <- pks
        chromPeakData(newFd) <- pkdl[[as.character(i)]]
      } else {
        chromPeaks(newFd) <- chromPeaks(x@msFeatureData)[0,]
        chromPeakData(newFd) <-
          .chrom_peak_data(x@msFeatureData)[0,]
      }
      lockEnvironment(newFd, bindings = TRUE)
      slot(a, "msFeatureData", check = FALSE) <- newFd
      res[[i]] <- a
    }
    res
  } else
    mapply(
      seq_along(fileNames(x)),
      fdl,
      FUN = create_object,
      MoreArgs = list(x = x, to_class = to_class)
    )
}



.chrom_peak_data <- function(x) {
  x$chromPeakData
}