findChromPeaks_sxt <-
  function(object,
           param,
           BPPARAM = bpparam(),
           return.type = "XCMSnExp",
           msLevel = 1L,
           ...) {
    return.type <- match.arg(return.type, c("XCMSnExp", "list",
                                            "xcmsSet"))
    # browser()
    startDate <- date()
    ## Restrict to MS X data.
    if (length(msLevel) > 1)
      stop("Currently only peak detection in a single MS level is ",
           "supported",
           call. = FALSE)
    ## Check if the data is centroided
    centroided <-
      all(MSnbase::centroided(object)[msLevel(object) %in% msLevel])
    if (is.na(centroided)) {
      idx <- which(msLevel(object) %in% msLevel)
      idx <- idx[ceiling(length(idx) / 3)]
      suppressWarnings(centroided <-
                         MSnbase::isCentroided(object[[idx]]))
    }
    if (is.na(centroided) || !centroided)
      warning("Your data appears to be not centroided! CentWave",
              " works best on data in centroid mode.")
    resList <-
      bplapply(
        split_by_file2_sxt(object, msLevel. = msLevel),
        FUN = findChromPeaks_OnDiskMSnExp_sxt,
        method = "centWave",
        param = param,
        BPPARAM = BPPARAM
      )
    ## (3) collect the results.
    res <- processResultList(resList,
                             getProcHist = return.type == "xcmsSet",
                             fnames = fileNames(object))
    
    if (return.type == "list")
      return(res$peaks)
    object <-
      peaks_to_result(
        res = res,
        object = object,
        startDate = startDate,
        param = param,
        msLevel = msLevel
      )
    if (return.type == "xcmsSet")
      as(object, "xcmsSet")
    else
      object
  }
