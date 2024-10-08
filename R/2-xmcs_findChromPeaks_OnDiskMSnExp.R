findChromPeaks_OnDiskMSnExp_sxt <-
  function(object, method = "centWave",
           param) {
    # require("xcms", quietly = TRUE, character.only = TRUE)
    if (!requireNamespace("xcms", quietly = TRUE)) {
      stop("Please install 'xcms' package first.")
    }
    if (missing(param))
      stop("'param' has to be specified!")
    file_name <-
      basename(object@processingData@files)
    ## pass the spectra to the _Spectrum_list function
    ## Since we're calling this function already with bplapply ensure that
    ## the spectra call is not firing its own parallel processing!
    findChromPeaks_Spectrum_list_sxt(
      x = spectra(object, BPPARAM = SerialParam()),
      method = method,
      param = param,
      rt = MSnbase::rtime(object),
      file_name = file_name
    )
  }



findChromPeaks_Spectrum_list_sxt <-
  function(x,
           method = "centWave",
           param,
           rt,
           file_name = "file") {
    method <-
      match.arg(
        method,
        c(
          "centWave",
          "massifquant",
          "matchedFilter",
          "MSW",
          "centWaveWithPredIsoROIs"
        )
      )
    method <- paste0("do_findChromPeaks_", method, "_sxt")
    if (method == "do_findChromPeaks_MSW")
      method <- "do_findPeaks_MSW"
    if (method == "do_findChromPeaks_matchedFilter") {
      ## Issue #325: empty spectra is not supported
      x <- lapply(x, function(z) {
        if (!length(z@mz)) {
          z@mz <- 0.0
          z@intensity <- 0.0
        }
        z
      })
    }
    if (missing(param))
      stop("'param' has to be specified!")
    if (missing(rt))
      rt <- unlist(lapply(x, rtime), use.names = FALSE)
    if (is.unsorted(rt))
      stop("Spectra are not ordered by retention time!")
    mzs <- lapply(x, mz)
    vals_per_spect <- lengths(mzs, FALSE)
    if (any(vals_per_spect == 0))
      warning("Found empty spectra. Please run 'filterEmptySpectra' first.",
              call. = FALSE)
    procDat <- date()
    
    ##------------------------------------------------------------------------
    do_findChromPeaks_centWave_sxt <-
      function(mz,
               int,
               scantime,
               valsPerSpect,
               ppm = 25,
               peakwidth = c(20, 50),
               snthresh = 10,
               prefilter = c(3, 100),
               mzCenterFun = "wMean",
               integrate = 1,
               mzdiff = -0.001,
               fitgauss = FALSE,
               noise = 0,
               verboseColumns = FALSE,
               roiList = list(),
               firstBaselineCheck = TRUE,
               roiScales = NULL,
               sleep = 0,
               extendLengthMSW = FALSE,
               file_name = "file") {
        centWave_orig_sxt <- function(mz,
                                      int,
                                      scantime,
                                      valsPerSpect,
                                      ppm = 25,
                                      peakwidth = c(20, 50),
                                      snthresh = 10,
                                      prefilter = c(3, 100),
                                      mzCenterFun = "wMean",
                                      integrate = 1,
                                      mzdiff = -0.001,
                                      fitgauss = FALSE,
                                      noise = 0,
                                      ## noise.local=TRUE,
                                      sleep = 0,
                                      verboseColumns = FALSE,
                                      roiList = list(),
                                      firstBaselineCheck = TRUE,
                                      roiScales = NULL,
                                      extendLengthMSW = FALSE,
                                      file_name = "file") {
          ## Input argument checking.
          if (missing(mz) |
              missing(int) |
              missing(scantime) | missing(valsPerSpect))
            stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
                 " are required!")
          if (length(mz) != length(int) |
              length(valsPerSpect) != length(scantime)
              | length(mz) != sum(valsPerSpect))
            stop(
              "Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
              " have to match. Also, 'length(mz)' should be equal to",
              " 'sum(valsPerSpect)'."
            )
          
          valueCount2ScanIndex <- function(valCount) {
            ## Convert into 0 based.
            valCount <- cumsum(valCount)
            return(as.integer(c(0, valCount[-length(valCount)])))
          }
          
          scanindex <-
            valueCount2ScanIndex(valsPerSpect) ## Get index vector for C calls
          if (!is.double(mz))
            mz <- as.double(mz)
          if (!is.double(int))
            int <- as.double(int)
          ## Fix the mzCenterFun
          mzCenterFun <- paste(
            "mzCenter",
            gsub(
              mzCenterFun,
              pattern = "mzCenter.",
              replacement = "",
              fixed = TRUE
            ),
            sep = "."
          )
          
          mzCenter.wMean <- function(mz, intensity) {
            weighted.mean(mz, intensity)
          }
          
          if (!exists(mzCenterFun, mode = "function"))
            stop("Function '", mzCenterFun, "' not defined !")
          
          if (!is.logical(firstBaselineCheck))
            stop("Parameter 'firstBaselineCheck' should be logical!")
          if (length(firstBaselineCheck) != 1)
            stop("Parameter 'firstBaselineCheck' should be a single logical !")
          if (length(roiScales) > 0)
            if (length(roiScales) != length(roiList) |
                !is.numeric(roiScales))
              stop(
                "If provided, parameter 'roiScales' has to be a numeric with",
                " length equal to the length of 'roiList'!"
              )
          ## if (!is.null(roiScales)) {
          ##     if (!is.numeric(roiScales) | length(roiScales) != length(roiList))
          ##         stop("Parameter 'roiScales' has to be a numeric of length equal to",
          ##              " parameter 'roiList'!")
          ##}
          
          basenames <- c("mz",
                         "mzmin",
                         "mzmax",
                         "rt",
                         "rtmin",
                         "rtmax",
                         "into",
                         "intb",
                         "maxo",
                         "sn")
          verbosenames <-
            c(
              "egauss",
              "mu",
              "sigma",
              "h",
              "f",
              "dppm",
              "scale",
              "scpos",
              "scmin",
              "scmax",
              "lmin",
              "lmax"
            )
          
          ## Peak width: seconds to scales
          scalerange <-
            round((peakwidth / mean(diff(scantime))) / 2)
          
          if (length(z <- which(scalerange == 0)))
            scalerange <- scalerange[-z]
          if (length(scalerange) < 1) {
            warning("No scales? Please check peak width!")
            if (verboseColumns) {
              nopeaks <- matrix(nrow = 0,
                                ncol = length(basenames) +
                                  length(verbosenames))
              colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
              nopeaks <- matrix(nrow = 0, ncol = length(basenames))
              colnames(nopeaks) <- c(basenames)
            }
            return(invisible(nopeaks))
          }
          
          if (length(scalerange) > 1) {
            scales <- seq(from = scalerange[1],
                          to = scalerange[2],
                          by = 2)
          } else{
            scales <- scalerange
          }
          
          minPeakWidth <-  scales[1]
          noiserange <- c(minPeakWidth * 3, max(scales) * 3)
          maxGaussOverlap <- 0.5
          minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
          minCentroids <- minPtsAboveBaseLine
          scRangeTol <-  maxDescOutlier <- floor(minPeakWidth / 2)
          scanrange <- c(1, length(scantime))
          
          ## If no ROIs are supplied then search for them.
          if (length(roiList) == 0) {
            message("Detecting mass traces at ",
                    ppm,
                    " ppm ... ",
                    file_name,
                    " ",
                    appendLF = FALSE)
            ## flush.console();
            ## We're including the findmzROI code in this function to reduce
            ## the need to copy objects etc.
            ## We could also sort the data by m/z anyway; wouldn't need that
            ## much time. Once we're using classes from MSnbase we can be
            ## sure that values are correctly sorted.
            withRestarts(
              tryCatch({
                tmp <- capture.output(
                  roiList <- .Call(
                    "findmzROI",
                    mz,
                    int,
                    scanindex,
                    as.double(c(0.0, 0.0)),
                    as.integer(scanrange),
                    as.integer(length(scantime)),
                    as.double(ppm * 1e-6),
                    as.integer(minCentroids),
                    as.integer(prefilter),
                    as.integer(noise),
                    PACKAGE = 'xcms'
                  )
                )
              },
              error = function(e) {
                if (grepl("m/z sort assumption violated !", e$message)) {
                  invokeRestart("fixSort")
                } else {
                  simpleError(e)
                }
              }),
              fixSort = function() {
                ## Force ordering of values within spectrum by mz:
                ##  o split values into a list -> mz per spectrum, intensity per
                ##    spectrum.
                ##  o define the ordering.
                ##  o re-order the mz and intensity and unlist again.
                ## Note: the Rle split is faster than the "conventional" factor split.
                splitF <- Rle(1:length(valsPerSpect), valsPerSpect)
                mzl <- as.list(S4Vectors::split(mz, f = splitF))
                oidx <- lapply(mzl, order)
                mz <<- unlist(
                  mapply(
                    mzl,
                    oidx,
                    FUN = function(y, z) {
                      return(y[z])
                    },
                    SIMPLIFY = FALSE,
                    USE.NAMES = FALSE
                  ),
                  use.names = FALSE
                )
                int <<- unlist(
                  mapply(
                    as.list(split(int, f = splitF)),
                    oidx,
                    FUN = function(y, z) {
                      return(y[z])
                    },
                    SIMPLIFY = FALSE,
                    USE.NAMES = FALSE
                  ),
                  use.names = FALSE
                )
                rm(mzl)
                rm(splitF)
                tmp <- capture.output(
                  roiList <<- .Call(
                    "findmzROI",
                    mz,
                    int,
                    scanindex,
                    as.double(c(0.0, 0.0)),
                    as.integer(scanrange),
                    as.integer(length(scantime)),
                    as.double(ppm * 1e-6),
                    as.integer(minCentroids),
                    as.integer(prefilter),
                    as.integer(noise),
                    PACKAGE = 'xcms'
                  )
                )
              }
            )
            message("OK")
            ## ROI.list <- findmzROI(object,scanrange=scanrange,dev=ppm * 1e-6,minCentroids=minCentroids, prefilter=prefilter, noise=noise)
            if (length(roiList) == 0) {
              warning("No ROIs found! \n")
              if (verboseColumns) {
                nopeaks <- matrix(nrow = 0,
                                  ncol = length(basenames) +
                                    length(verbosenames))
                colnames(nopeaks) <- c(basenames, verbosenames)
              } else {
                nopeaks <- matrix(nrow = 0, ncol = length(basenames))
                colnames(nopeaks) <- c(basenames)
              }
              return(invisible(nopeaks))
            }
          }
          
          ## Second stage: process the ROIs
          peaklist <- list()
          Nscantime <- length(scantime)
          lf <- length(roiList)
          
          ## cat('\n Detecting chromatographic peaks ... \n % finished: ')
          ## lp <- -1
          message(
            "Detecting chromatographic peaks in ",
            length(roiList),
            " regions of interest ...",
            appendLF = FALSE
          )
          
          for (f in  1:lf) {
            ## ## Show progress
            ## perc <- round((f/lf) * 100)
            ## if ((perc %% 10 == 0) && (perc != lp))
            ## {
            ##     cat(perc," ",sep="");
            ##     lp <- perc;
            ## }
            ## flush.console()
            
            feat <- roiList[[f]]
            N <- feat$scmax - feat$scmin + 1
            peaks <- peakinfo <- NULL
            mzrange <- c(feat$mzmin, feat$mzmax)
            sccenter <- feat$scmin[1] + floor(N / 2) - 1
            scrange <- c(feat$scmin, feat$scmax)
            ## scrange + noiserange, used for baseline detection and wavelet analysis
            sr <- c(max(scanrange[1], scrange[1] - max(noiserange)),
                    min(scanrange[2], scrange[2] + max(noiserange)))
            eic <- .Call(
              "getEIC",
              mz,
              int,
              scanindex,
              as.double(mzrange),
              as.integer(sr),
              as.integer(length(scanindex)),
              PACKAGE = "xcms"
            )
            ## eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
            d <- eic$intensity
            td <- sr[1]:sr[2]
            scan.range <- c(sr[1], sr[2])
            ## original mzROI range
            idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
            mzROI.EIC <-
              list(scan = eic$scan[idxs],
                   intensity = eic$intensity[idxs])
            ## mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
            omz <-
              .Call(
                "getWeightedMZ",
                mz,
                int,
                scanindex,
                as.double(mzrange),
                as.integer(scrange),
                as.integer(length(scantime)),
                PACKAGE = 'xcms'
              )
            ## omz <- rawMZ(object,mzrange=mzrange,scanrange=scrange)
            if (all(omz == 0)) {
              warning("centWave: no peaks found in ROI.")
              next
            }
            od  <- mzROI.EIC$intensity
            otd <- mzROI.EIC$scan
            if (all(od == 0)) {
              warning("centWave: no peaks found in ROI.")
              next
            }
            
            ## scrange + scRangeTol, used for gauss fitting and continuous
            ## data above 1st baseline detection
            ftd <-
              max(td[1], scrange[1] - scRangeTol):min(td[length(td)],
                                                      scrange[2] + scRangeTol)
            fd <- d[match(ftd, td)]
            
            ## 1st type of baseline: statistic approach
            if (N >= 10 * minPeakWidth) {
              ## in case of very long mass trace use full scan range
              ## for baseline detection
              noised <-
                .Call(
                  "getEIC",
                  mz,
                  int,
                  scanindex,
                  as.double(mzrange),
                  as.integer(scanrange),
                  as.integer(length(scanindex)),
                  PACKAGE = "xcms"
                )$intensity
              ## noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity
            } else {
              noised <- d
            }
            ## 90% trimmed mean as first baseline guess
            estimateChromNoise <-
              function(x,
                       trim = 0.05,
                       minPts = 20) {
                gz <- which(x > 0)
                if (length(gz) < minPts)
                  return(mean(x))
                
                mean(x[gz], trim = trim)
              }
            noise <- estimateChromNoise(noised,
                                        trim = 0.05,
                                        minPts = 3 * minPeakWidth)
            ## any continuous data above 1st baseline ?
            continuousPtsAboveThreshold <-
              function(y, threshold, num, istart = 1) {
                if (!is.double(y))
                  y <- as.double(y)
                if (.C(
                  "continuousPtsAboveThreshold",
                  y,
                  as.integer(istart - 1),
                  length(y),
                  threshold = as.double(threshold),
                  num = as.integer(num),
                  n = integer(1),
                  PACKAGE = "xcms"
                )$n > 0)
                  TRUE
                else
                  FALSE
              }
            if (firstBaselineCheck &&
                !continuousPtsAboveThreshold(fd, threshold = noise,
                                             num = minPtsAboveBaseLine))
              next
            ## 2nd baseline estimate using not-peak-range
            getLocalNoiseEstimate <-
              function(d,
                       td,
                       ftd,
                       noiserange,
                       Nscantime,
                       threshold,
                       num) {
                if (length(d) < Nscantime) {
                  ## noiserange[2] is full d-range
                  drange <- which(td %in% ftd)
                  n1 <-
                    d[-drange] ## region outside the detected ROI (wide)
                  continuousPtsAboveThresholdIdx <-
                    function(y, threshold, num, istart = 1) {
                      if (!is.double(y))
                        y <- as.double(y)
                      as.logical(
                        .C(
                          "continuousPtsAboveThresholdIdx",
                          y,
                          as.integer(istart - 1),
                          length(y),
                          threshold = as.double(threshold),
                          num = as.integer(num),
                          n = integer(length(y)),
                          PACKAGE = "xcms"
                        )$n
                      )
                    }
                  n1.cp <-
                    continuousPtsAboveThresholdIdx(n1, threshold = threshold, num = num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
                  n1 <- n1[!n1.cp]
                  if (length(n1) > 1)  {
                    baseline1 <- mean(n1)
                    sdnoise1 <- sd(n1)
                  } else
                    baseline1 <- sdnoise1 <- 1
                  
                  ## noiserange[1]
                  d1 <- drange[1]
                  d2 <- drange[length(drange)]
                  nrange2 <-
                    c(max(1, d1 - noiserange[1]):d1,
                      d2:min(length(d), d2 + noiserange[1]))
                  n2 <-
                    d[nrange2] ## region outside the detected ROI (narrow)
                  n2.cp <-
                    continuousPtsAboveThresholdIdx(n2, threshold = threshold, num = num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
                  n2 <- n2[!n2.cp]
                  if (length(n2) > 1)  {
                    baseline2 <- mean(n2)
                    sdnoise2 <- sd(n2)
                  } else
                    baseline2 <- sdnoise2 <- 1
                  
                } else {
                  trimm <- function(x, trim = c(0.05, 0.95)) {
                    a <- sort(x[x > 0])
                    Na <- length(a)
                    quant <- round((Na * trim[1]) + 1):round(Na * trim[2])
                    a[quant]
                  }
                  trimmed <- trimm(d, c(0.05, 0.95))
                  baseline1 <- baseline2 <- mean(trimmed)
                  sdnoise1 <- sdnoise2 <- sd(trimmed)
                }
                
                c(min(baseline1, baseline2), min(sdnoise1, sdnoise2))
              }
            lnoise <-
              getLocalNoiseEstimate(d,
                                    td,
                                    ftd,
                                    noiserange,
                                    Nscantime,
                                    threshold = noise,
                                    num = minPtsAboveBaseLine)
            ## Final baseline & Noise estimate
            baseline <- max(1, min(lnoise[1], noise))
            sdnoise <- max(1, lnoise[2])
            sdthr <-  sdnoise * snthresh
            ## is there any data above S/N * threshold ?
            if (!(any(fd - baseline >= sdthr)))
              next
            
            MSW.cwt <-
              function (ms,
                        scales = 1,
                        wavelet = "mexh",
                        extendLengthMSW = FALSE) {
                if (wavelet == "mexh") {
                  psi_xval <- seq(-6, 6, length = 256)
                  psi <-
                    (2 / sqrt(3) * pi ^ (-0.25)) * (1 - psi_xval ^ 2) *
                    exp(-psi_xval ^ 2 / 2)
                }
                else if (is.matrix(wavelet)) {
                  if (nrow(wavelet) == 2) {
                    psi_xval <- wavelet[1,]
                    psi <- wavelet[2,]
                  }
                  else if (ncol(wavelet) == 2) {
                    psi_xval <- wavelet[, 1]
                    psi <- wavelet[, 2]
                  }
                  else {
                    stop("Unsupported wavelet format!")
                  }
                }
                else {
                  stop("Unsupported wavelet!")
                }
                oldLen <- length(ms)
                # IF extendLengthMSW is TRUE:
                # The new length is determined by the scales argument, so a larger peakwidth
                # will ensure more scales are run, but may slow it down. See
                # https://github.com/sneumann/xcms/issues/445 for more information about
                # a change from using extendNBase to extendLength.
                if (extendLengthMSW) {
                  newLen <- 2 ^ (ceiling(log2(max(scales) * 12)))
                  ms <-
                    MSW.extendLength(
                      x = ms,
                      addLength = (newLen - length(ms)),
                      method = "open"
                    )
                } else {
                  MSW.extendNBase <- function(x,
                                              nLevel = 1,
                                              base = 2,
                                              ...)
                  {
                    ## from package MassSpecWavelet
                    if (!is.matrix(x))
                      x <- matrix(x, ncol = 1)
                    
                    nR <- nrow(x)
                    if (is.null(nLevel)) {
                      nR1 <- nextn(nR, base)
                    } else {
                      nR1 <- ceiling(nR / base ^ nLevel) * base ^ nLevel
                    }
                    if (nR != nR1) {
                      MSW.extendLength <-
                        function(x,
                                 addLength = NULL,
                                 method = c('reflection', 'open', 'circular'),
                                 direction = c('right', 'left', 'both'))
                        {
                          ## from package MassSpecWavelet
                          if (is.null(addLength))
                            stop('Please provide the length to be added!')
                          if (!is.matrix(x))
                            x <- matrix(x, ncol = 1)
                          method <- match.arg(method)
                          direction <- match.arg(direction)
                          
                          nR <- nrow(x)
                          nR1 <- nR + addLength
                          if (direction == 'both') {
                            left <- right <- addLength
                          } else if (direction == 'right') {
                            left <- 0
                            right <- addLength
                          } else if (direction == 'left') {
                            left <- addLength
                            right <- 0
                          }
                          
                          if (right > 0) {
                            x <- switch(
                              method,
                              reflection = rbind(x, x[nR:(2 * nR - nR1 + 1), , drop =
                                                        FALSE]),
                              open = rbind(x, matrix(
                                rep(x[nR, ], addLength),
                                ncol = ncol(x),
                                byrow = TRUE
                              )),
                              circular = rbind(x, x[1:(nR1 - nR), , drop = FALSE])
                            )
                          }
                          
                          if (left > 0) {
                            x <- switch(
                              method,
                              reflection = rbind(x[addLength:1, , drop = FALSE], x),
                              open = rbind(matrix(
                                rep(x[1, ], addLength),
                                ncol = ncol(x),
                                byrow = TRUE
                              ), x),
                              circular = rbind(x[(2 * nR - nR1 + 1):nR, , drop = FALSE], x)
                            )
                          }
                          if (ncol(x) == 1)
                            x <- as.vector(x)
                          
                          x
                        }
                      x <-
                        MSW.extendLength(x, addLength = nR1 - nR, ...)
                    }
                    x
                  }
                  ms <- MSW.extendNBase(ms, nLevel = NULL, base = 2)
                }
                
                
                len <- length(ms)
                nbscales <- length(scales)
                wCoefs <- NULL
                psi_xval <- psi_xval - psi_xval[1]
                dxval <- psi_xval[2]
                xmax <- psi_xval[length(psi_xval)]
                for (i in 1:length(scales)) {
                  scale.i <- scales[i]
                  f <- rep(0, len)
                  j <-
                    1 + floor((0:(scale.i * xmax)) / (scale.i * dxval))
                  if (length(j) == 1)
                    j <- c(1, 1)
                  lenWave <- length(j)
                  f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
                  if (length(f) > len)
                  {
                    i <- i - 1
                    break
                  }   ##  stop(paste("scale", scale.i, "is too large!"))
                  wCoefs.i <- 1 / sqrt(scale.i) * convolve(ms, f)
                  wCoefs.i <-
                    c(wCoefs.i[(len - floor(lenWave / 2) + 1):len],
                      wCoefs.i[1:(len - floor(lenWave / 2))])
                  wCoefs <- cbind(wCoefs, wCoefs.i)
                }
                if (i < 1)
                  return(NA)
                scales <- scales[1:i]
                if (length(scales) == 1)
                  wCoefs <- matrix(wCoefs, ncol = 1)
                colnames(wCoefs) <- scales
                wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
                wCoefs
              }
            
            wCoefs <- MSW.cwt(
              d,
              scales = scales,
              wavelet = 'mexh',
              extendLengthMSW = extendLengthMSW
            )
            if (!(!is.null(dim(wCoefs)) &&
                  any(wCoefs - baseline >= sdthr)))
              next
            if (td[length(td)] == Nscantime)
              ## workaround, localMax fails otherwise
              wCoefs[nrow(wCoefs), ] <-
              wCoefs[nrow(wCoefs) - 1,] * 0.99
            MSW.getLocalMaximumCWT <-
              function(wCoefs,
                       minWinSize = 5,
                       amp.Th = 0)
              {
                ## from package MassSpecWavelet
                localMax <- NULL
                scales <- as.numeric(colnames(wCoefs))
                
                for (i in 1:length(scales)) {
                  scale.i <- scales[i]
                  winSize.i <- scale.i * 2 + 1
                  if (winSize.i < minWinSize) {
                    winSize.i <- minWinSize
                  }
                  MSW.localMaximum <-
                    function (x, winSize = 5)
                    {
                      ## from package MassSpecWavelet
                      len <- length(x)
                      rNum <- ceiling(len / winSize)
                      
                      ## Transform the vector as a matrix with column length equals winSize
                      ##		and find the maximum position at each row.
                      y <-
                        matrix(c(x, rep(x[len], rNum * winSize - len)), nrow = winSize)
                      y.maxInd <- apply(y, 2, which.max)
                      ## Only keep the maximum value larger than the boundary values
                      selInd <-
                        which(apply(y, 2, function(x)
                          max(x) > x[1] & max(x) > x[winSize]))
                      
                      ## keep the result
                      localMax <- rep(0, len)
                      localMax[(selInd - 1) * winSize + y.maxInd[selInd]] <- 1
                      
                      ## Shift the vector with winSize/2 and do the same operation
                      shift <- floor(winSize / 2)
                      rNum <- ceiling((len + shift) / winSize)
                      y <-
                        matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow =
                                 winSize)
                      y.maxInd <- apply(y, 2, which.max)
                      ## Only keep the maximum value larger than the boundary values
                      selInd <-
                        which(apply(y, 2, function(x)
                          max(x) > x[1] & max(x) > x[winSize]))
                      localMax[(selInd - 1) * winSize + y.maxInd[selInd] - shift] <- 1
                      
                      ## Check whether there is some local maxima have in between distance less than winSize
                      maxInd <- which(localMax > 0)
                      selInd <- which(diff(maxInd) < winSize)
                      if (length(selInd) > 0) {
                        selMaxInd1 <- maxInd[selInd]
                        selMaxInd2 <- maxInd[selInd + 1]
                        temp <- x[selMaxInd1] - x[selMaxInd2]
                        localMax[selMaxInd1[temp <= 0]] <- 0
                        localMax[selMaxInd2[temp > 0]] <- 0
                      }
                      
                      localMax
                    }
                  temp <- MSW.localMaximum(wCoefs[, i], winSize.i)
                  localMax <- cbind(localMax, temp)
                }
                ## Set the values less than peak threshold as 0
                localMax[wCoefs < amp.Th] <- 0
                colnames(localMax) <- colnames(wCoefs)
                rownames(localMax) <- rownames(wCoefs)
                localMax
              }
            localMax <- MSW.getLocalMaximumCWT(wCoefs)
            MSW.getRidge <-
              function(localMax,
                       iInit = ncol(localMax),
                       step = -1,
                       iFinal = 1,
                       minWinSize = 3,
                       gapTh = 3,
                       skip = NULL)
              {
                ## modified from package MassSpecWavelet
                
                scales <- as.numeric(colnames(localMax))
                if (is.null(scales))
                  scales <- 1:ncol(localMax)
                
                maxInd_curr <- which(localMax[, iInit] > 0)
                nMz <- nrow(localMax)
                
                if (is.null(skip))	{
                  skip <- iInit + 1
                }
                
                ## Identify all the peak pathes from the coarse level to detail levels (high column to low column)
                ## Only consider the shortest path
                if (ncol(localMax) > 1)
                  colInd <- seq(iInit + step, iFinal, step)
                else
                  colInd <- 1
                ridgeList <- as.list(maxInd_curr)
                names(ridgeList) <- maxInd_curr
                peakStatus <- as.list(rep(0, length(maxInd_curr)))
                names(peakStatus) <- maxInd_curr
                
                ## orphanRidgeList keep the ridges disconnected at certain scale level
                ## Changed by Pan Du 05/11/06
                orphanRidgeList <- NULL
                orphanRidgeName <- NULL
                nLevel <- length(colInd)
                
                for (j in 1:nLevel) {
                  col.j <- colInd[j]
                  scale.j <- scales[col.j]
                  
                  if (colInd[j] == skip) {
                    oldname <- names(ridgeList)
                    ridgeList <-
                      lapply(ridgeList, function(x)
                        c(x, x[length(x)]))
                    ##peakStatus <- lapply(peakStatus, function(x) c(x, x[length(x)]))
                    names(ridgeList) <- oldname
                    ##names(peakStatus) <- oldname
                    next
                  }
                  
                  if (length(maxInd_curr) == 0) {
                    maxInd_curr <- which(localMax[, col.j] > 0)
                    next
                  }
                  
                  ## The slide window size is proportional to the CWT scale
                  ## winSize.j <- scale.j / 2 + 1
                  winSize.j <- floor(scale.j / 2)
                  if (winSize.j < minWinSize) {
                    winSize.j <- minWinSize
                  }
                  
                  selPeak.j <- NULL
                  remove.j <- NULL
                  for (k in 1:length(maxInd_curr)) {
                    ind.k <- maxInd_curr[k]
                    start.k <-
                      ifelse(ind.k - winSize.j < 1, 1, ind.k - winSize.j)
                    end.k <-
                      ifelse(ind.k + winSize.j > nMz, nMz, ind.k + winSize.j)
                    ind.curr <-
                      which(localMax[start.k:end.k, col.j] > 0) + start.k - 1
                    ##ind.curr <- which(localMax[, col.j] > 0)
                    if (length(ind.curr) == 0) {
                      status.k <- peakStatus[[as.character(ind.k)]]
                      ## bug  work-around
                      if (is.null(status.k))
                        status.k <- gapTh + 1
                      ##
                      if (status.k > gapTh & scale.j >= 2) {
                        temp <- ridgeList[[as.character(ind.k)]]
                        orphanRidgeList <-
                          c(orphanRidgeList, list(temp[1:(length(temp) - status.k)]))
                        orphanRidgeName <-
                          c(orphanRidgeName,
                            paste(col.j + status.k + 1, ind.k, sep = '_'))
                        remove.j <- c(remove.j, as.character(ind.k))
                        next
                      } else {
                        ind.curr <- ind.k
                        peakStatus[[as.character(ind.k)]] <-
                          status.k + 1
                      }
                    } else {
                      peakStatus[[as.character(ind.k)]] <- 0
                      if (length(ind.curr) >= 2)
                        ind.curr <- ind.curr[which.min(abs(ind.curr - ind.k))]
                    }
                    ridgeList[[as.character(ind.k)]] <-
                      c(ridgeList[[as.character(ind.k)]], ind.curr)
                    selPeak.j <- c(selPeak.j, ind.curr)
                  }
                  ## Remove the disconnected lines from the currrent list
                  if (length(remove.j) > 0) {
                    removeInd <- which(names(ridgeList) %in% remove.j)
                    ridgeList <- ridgeList[-removeInd]
                    peakStatus <- peakStatus[-removeInd]
                  }
                  
                  ## Check for duplicated selected peaks and only keep the one with the longest path.
                  dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
                  if (length(dupPeak.j) > 0) {
                    removeInd <- NULL
                    for (dupPeak.jk in dupPeak.j) {
                      selInd <- which(selPeak.j == dupPeak.jk)
                      selLen <- sapply(ridgeList[selInd], length)
                      removeInd.jk <- which.max(selLen)
                      removeInd <- c(removeInd, selInd[-removeInd.jk])
                      orphanRidgeList <-
                        c(orphanRidgeList, ridgeList[removeInd.jk])
                      orphanRidgeName <-
                        c(orphanRidgeName, paste(col.j, selPeak.j[removeInd.jk], sep = '_'))
                    }
                    selPeak.j <- selPeak.j[-removeInd]
                    ridgeList <- ridgeList[-removeInd]
                    peakStatus <- peakStatus[-removeInd]
                  }
                  
                  ## Update the names of the ridgeList as the new selected peaks
                  ##if (scale.j >= 2) {
                  if (length(ridgeList) > 0)
                    names(ridgeList) <- selPeak.j
                  if (length(peakStatus) > 0)
                    names(peakStatus) <- selPeak.j
                  ##}
                  
                  ## If the level is larger than 3, expand the peak list by including other unselected peaks at that level
                  if (scale.j >= 2) {
                    maxInd_next <- which(localMax[, col.j] > 0)
                    unSelPeak.j <-
                      maxInd_next[!(maxInd_next %in% selPeak.j)]
                    newPeak.j <- as.list(unSelPeak.j)
                    names(newPeak.j) <- unSelPeak.j
                    ## Update ridgeList
                    ridgeList <- c(ridgeList, newPeak.j)
                    maxInd_curr <- c(selPeak.j, unSelPeak.j)
                    ## Update peakStatus
                    newPeakStatus <- as.list(rep(0, length(newPeak.j)))
                    names(newPeakStatus) <- newPeak.j
                    peakStatus <- c(peakStatus, newPeakStatus)
                  } else {
                    maxInd_curr <- selPeak.j
                  }
                }
                
                ## Attach the peak level at the beginning of the ridge names
                if (length(ridgeList) > 0)
                  names(ridgeList) <- paste(1, names(ridgeList), sep = '_')
                if (length(orphanRidgeList) > 0)
                  names(orphanRidgeList) <- orphanRidgeName
                ## Combine ridgeList and orphanRidgeList
                ridgeList <- c(ridgeList, orphanRidgeList)
                if (length(ridgeList) == 0)
                  return(NULL)
                
                ## Reverse the order as from the low level to high level.
                ridgeList <- lapply(ridgeList, rev)
                ## order the ridgeList in increasing order
                ##ord <- order(selPeak.j)
                ##ridgeList <- ridgeList[ord]
                
                ## Remove possible duplicated ridges
                ridgeList <- ridgeList[!duplicated(names(ridgeList))]
                
                attr(ridgeList, 'class') <- 'ridgeList'
                attr(ridgeList, 'scales') <- scales
                return(ridgeList)
              }
            rL <- MSW.getRidge(localMax)
            wpeaks <- sapply(rL,
                             function(x) {
                               w <- min(1:length(x), ncol(wCoefs))
                               any(wCoefs[x, w] - baseline >= sdthr)
                             })
            if (any(wpeaks)) {
              wpeaksidx <- which(wpeaks)
              ## check each peak in ridgeList
              for (p in 1:length(wpeaksidx)) {
                opp <- rL[[wpeaksidx[p]]]
                pp <- unique(opp)
                if (length(pp) >= 1) {
                  dv <- td[pp] %in% ftd
                  if (any(dv)) {
                    ## peaks in orig. data range
                    ## Final S/N check
                    if (any(d[pp[dv]] - baseline >= sdthr)) {
                      ## if(!is.null(roiScales)) {
                      ## allow roiScales to be a numeric of length 0
                      if (length(roiScales) > 0) {
                        ## use given scale
                        best.scale.nr <-
                          which(scales == roiScales[[f]])
                        if (best.scale.nr > length(opp))
                          best.scale.nr <- length(opp)
                      } else {
                        ## try to decide which scale describes the peak best
                        inti <- numeric(length(opp))
                        irange <-
                          rep(ceiling(scales[1] / 2), length(opp))
                        for (k in 1:length(opp)) {
                          kpos <- opp[k]
                          r1 <- ifelse(kpos - irange[k] > 1,
                                       kpos - irange[k], 1)
                          r2 <- ifelse(kpos + irange[k] < length(d),
                                       kpos + irange[k],
                                       length(d))
                          inti[k] <- sum(d[r1:r2])
                        }
                        maxpi <- which.max(inti)
                        if (length(maxpi) > 1) {
                          m <- wCoefs[opp[maxpi], maxpi]
                          bestcol <- which(m == max(m),
                                           arr.ind = TRUE)[2]
                          best.scale.nr <- maxpi[bestcol]
                        } else
                          best.scale.nr <- maxpi
                      }
                      
                      best.scale <-  scales[best.scale.nr]
                      best.scale.pos <- opp[best.scale.nr]
                      
                      pprange <- min(pp):max(pp)
                      ## maxint <- max(d[pprange])
                      lwpos <- max(1, best.scale.pos - best.scale)
                      rwpos <-
                        min(best.scale.pos + best.scale, length(td))
                      p1 <- match(td[lwpos], otd)[1]
                      p2 <- match(td[rwpos], otd)
                      p2 <- p2[length(p2)]
                      if (is.na(p1))
                        p1 <- 1
                      if (is.na(p2))
                        p2 <- N
                      mz.value <- omz[p1:p2]
                      mz.int <- od[p1:p2]
                      maxint <- max(mz.int)
                      
                      ## re-calculate m/z value for peak range
                      mzrange <- range(mz.value)
                      mzmean <- do.call(mzCenterFun,
                                        list(mz = mz.value,
                                             intensity = mz.int))
                      
                      ## Compute dppm only if needed
                      dppm <- NA
                      if (verboseColumns) {
                        if (length(mz.value) >= (minCentroids + 1)) {
                          dppm <- round(min(
                            running(
                              abs(diff(mz.value)) /
                                (mzrange[2] *  1e-6),
                              fun = max,
                              width = minCentroids
                            )
                          ))
                        } else {
                          dppm <- round((mzrange[2] - mzrange[1]) /
                                          (mzrange[2] * 1e-6))
                        }
                      }
                      peaks <- rbind(
                        peaks,
                        c(
                          mzmean,
                          mzrange,
                          ## mz
                          NA,
                          NA,
                          NA,
                          ## rt, rtmin, rtmax,
                          NA,
                          ## intensity (sum)
                          NA,
                          ## intensity (-bl)
                          maxint,
                          ## max intensity
                          round((maxint - baseline) / sdnoise),
                          ##  S/N Ratio
                          NA,
                          ## Gaussian RMSE
                          NA,
                          NA,
                          NA,
                          ## Gaussian Parameters
                          f,
                          ## ROI Position
                          dppm,
                          ## max. difference between the [minCentroids] peaks in ppm
                          best.scale,
                          ## Scale
                          td[best.scale.pos],
                          td[lwpos],
                          td[rwpos],
                          ## Peak positions guessed from the wavelet's (scan nr)
                          NA,
                          NA
                        )
                      )                    ## Peak limits (scan nr)
                      peakinfo <- rbind(
                        peakinfo,
                        c(
                          best.scale,
                          best.scale.nr,
                          best.scale.pos,
                          lwpos,
                          rwpos
                        )
                      )
                      ## Peak positions guessed from the wavelet's
                    }
                  }
                }
              }  ##for
            } ## if
            
            ##  postprocessing
            if (!is.null(peaks)) {
              colnames(peaks) <- c(basenames, verbosenames)
              colnames(peakinfo) <- c("scale", "scaleNr", "scpos",
                                      "scmin", "scmax")
              descendMinTol <- function(d, startpos, maxDescOutlier) {
                l <- startpos[1]
                r <- startpos[2]
                outl <- 0
                N <- length(d)
                ## left
                while ((l > 1) && (d[l] > 0) && outl <= maxDescOutlier) {
                  if (outl > 0)
                    vpos <- opos
                  else
                    vpos <- l
                  if (d[l - 1] > d[vpos])
                    outl <- outl + 1
                  else
                    outl <- 0
                  if (outl == 1)
                    opos <- l
                  l <- l - 1
                }
                if (outl > 0)
                  l <- l + outl
                ## right
                outl <- 0
                
                while ((r < N) && (d[r] > 0) && outl <= maxDescOutlier) {
                  if (outl > 0)
                    vpos <- opos
                  else
                    vpos <- r
                  if (d[r + 1] > d[vpos])
                    outl <- outl + 1
                  else
                    outl <- 0
                  if (outl == 1)
                    opos <- r
                  r <- r + 1
                }
                if (outl > 0)
                  r <- r - outl
                c(l, r)
              }
              for (p in 1:dim(peaks)[1]) {
                ## find minima (peak boundaries), assign rt and intensity values
                if (integrate == 1) {
                  lm <- descendMin(wCoefs[, peakinfo[p, "scaleNr"]],
                                   istart = peakinfo[p, "scpos"])
                  gap <-
                    all(d[lm[1]:lm[2]] == 0) # looks like we got stuck in a gap right in the middle of the peak
                  if ((lm[1] == lm[2]) || gap)
                    # fall-back
                    
                    lm <-
                    descendMinTol(d,
                                  startpos = c(peakinfo[p, "scmin"],
                                               peakinfo[p, "scmax"]),
                                  maxDescOutlier)
                } else {
                  lm <- descendMinTol(d,
                                      startpos = c(peakinfo[p, "scmin"],
                                                   peakinfo[p, "scmax"]),
                                      maxDescOutlier)
                }
                ## Narrow peak rt boundaries by removing values below threshold
                narrow_rt_boundaries <- function(lm, d, thresh = 1) {
                  lm_seq <- lm[1]:lm[2]
                  above_thresh <- d[lm_seq] >= thresh
                  if (any(above_thresh)) {
                    ## Expand by one on each side to be consistent with old code.
                    above_thresh <- above_thresh | c(above_thresh[-1], FALSE) |
                      c(FALSE, above_thresh[-length(above_thresh)])
                    lm <- range(lm_seq[above_thresh], na.rm = TRUE)
                  }
                  lm
                }
                lm <- narrow_rt_boundaries(lm, d)
                lm_seq <- lm[1]:lm[2]
                pd <- d[lm_seq]
                
                peakrange <- td[lm]
                peaks[p, "rtmin"] <- scantime[peakrange[1]]
                peaks[p, "rtmax"] <- scantime[peakrange[2]]
                peaks[p, "maxo"] <- max(pd)
                pwid <-
                  (scantime[peakrange[2]] - scantime[peakrange[1]]) /
                  (peakrange[2] - peakrange[1])
                if (is.na(pwid))
                  pwid <- 1
                peaks[p, "into"] <- pwid * sum(pd)
                db <- pd - baseline
                peaks[p, "intb"] <- pwid * sum(db[db > 0])
                peaks[p, "lmin"] <- lm[1]
                peaks[p, "lmax"] <- lm[2]
                
                if (fitgauss) {
                  ## perform gaussian fits, use wavelets for inital parameters
                  td_lm <- td[lm_seq]
                  md <- max(pd)
                  d1 <-
                    pd / md ## normalize data for gaussian error calc.
                  pgauss <- fitGauss(td_lm,
                                     pd,
                                     pgauss = list(
                                       mu = peaks[p, "scpos"],
                                       sigma = peaks[p, "scmax"] -
                                         peaks[p, "scmin"],
                                       h = peaks[p, "maxo"]
                                     ))
                  rtime <- peaks[p, "scpos"]
                  if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                    gtime <- td[match(round(pgauss$mu), td)]
                    if (!is.na(gtime)) {
                      rtime <- gtime
                      peaks[p, "mu"] <- pgauss$mu
                      peaks[p, "sigma"] <- pgauss$sigma
                      peaks[p, "h"] <- pgauss$h
                      peaks[p, "egauss"] <-
                        sqrt((1 / length(td_lm)) *
                               sum(((
                                 d1 - gauss(td_lm, pgauss$h / md,
                                            pgauss$mu, pgauss$sigma)
                               ) ^ 2)))
                    }
                  }
                  peaks[p, "rt"] <- scantime[rtime]
                  ## avoid fitting side effects
                  if (peaks[p, "rt"] < peaks[p, "rtmin"])
                    peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
                } else
                  peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
              }
              joinOverlappingPeaks <-
                function(td,
                         d,
                         otd,
                         omz,
                         od,
                         scantime,
                         scan.range,
                         peaks,
                         maxGaussOverlap = 0.5,
                         mzCenterFun) {
                  ## Fix issue #284: avoid having identical peaks multiple times in this
                  ## matrix.
                  peaks <- unique(peaks)
                  gausspeaksidx <- which(!is.na(peaks[, "mu"]))
                  Ngp <- length(gausspeaksidx)
                  if (Ngp == 0)
                    return(peaks)
                  
                  newpeaks <- NULL
                  
                  gpeaks <- peaks[gausspeaksidx, , drop = FALSE]
                  if (nrow(peaks) - Ngp > 0)
                    notgausspeaks <- peaks[-gausspeaksidx, , drop = FALSE]
                  
                  if (Ngp > 1) {
                    comb <- which(upper.tri(matrix(0, Ngp, Ngp)), arr.ind = TRUE)
                    overlap <- logical(nrow(comb))
                    overlap <- rep(FALSE, dim(comb)[1])
                    for (k in seq_len(nrow(comb))) {
                      p1 <- comb[k, 1]
                      p2 <- comb[k, 2]
                      overlap[k] <- gaussCoverage(
                        xlim = scan.range,
                        h1 = gpeaks[p1, "h"],
                        mu1 = gpeaks[p1, "mu"],
                        s1 = gpeaks[p1, "sigma"],
                        h2 = gpeaks[p2, "h"],
                        mu2 = gpeaks[p2, "mu"],
                        s2 = gpeaks[p2, "sigma"]
                      ) >=
                        maxGaussOverlap
                    }
                  } else
                    overlap <- FALSE
                  
                  if (any(overlap) && (Ngp > 1)) {
                    jlist <- list()
                    if (length(which(overlap)) > 1) {
                      gm <- comb[overlap, ]
                      ## create list of connected components
                      cc <- list()
                      cc[[1]] <- gm[1,] ## copy first entry
                      for (j in 2:dim(gm)[1]) {
                        ## search for connections
                        ccl <- unlist(cc)
                        nl <- sapply(cc, function(x)
                          length(x))
                        ccidx <- rep(1:length(nl), nl)
                        idx <- match(gm[j,], ccl)
                        if (any(!is.na(idx))) {
                          ## connection found, add to list
                          pos <- ccidx[idx[which(!is.na(idx))[1]]]
                          cc[[pos]] <- c(cc[[pos]], gm[j,])
                        } else
                          ## create new list element
                          cc[[length(cc) + 1]] <- gm[j,]
                        
                      }
                      ccn <- list()
                      lcc <- length(cc)
                      ins <- rep(FALSE, lcc)
                      if (lcc > 1) {
                        jcomb <- which(upper.tri(matrix(0, lcc, lcc)), arr.ind = TRUE)
                        for (j in 1:dim(jcomb)[1]) {
                          j1 <- jcomb[j, 1]
                          j2 <- jcomb[j, 2]
                          if (any(cc[[j1]] %in% cc[[j2]]))
                            ccn[[length(ccn) + 1]] <-
                            unique(c(cc[[j1]], cc[[j2]]))
                          else {
                            if (!ins[j1]) {
                              ccn[[length(ccn) + 1]] <- unique(cc[[j1]])
                              ins[j1] <- TRUE
                            }
                            if (!ins[j2]) {
                              ccn[[length(ccn) + 1]] <- unique(cc[[j2]])
                              ins[j2] <- TRUE
                            }
                          }
                        }
                      } else
                        ccn <- cc
                      
                      
                      size <- sapply(ccn, function(x)
                        length(x))
                      s2idx <- which(size >= 2)
                      
                      if (length(s2idx) > 0) {
                        for (j in 1:length(s2idx)) {
                          pgroup <- unique(ccn[[s2idx[j]]])
                          jlist[[j]] <- pgroup
                        }
                      } else
                        stop('(length(s2idx) = 0) ?!?')
                    } else
                      jlist[[1]] <- comb[overlap, ]
                    
                    ## join all peaks belonging to one cc
                    for (j in seq_along(jlist)) {
                      jidx <- jlist[[j]]
                      newpeak <- gpeaks[jidx[1], , drop = FALSE]
                      newmin <- min(gpeaks[jidx, "lmin"])
                      newmax <- max(gpeaks[jidx, "lmax"])
                      newpeak[1, "scpos"] <- -1 ## not defined after join
                      newpeak[1, "scmin"] <- -1 ##    ..
                      newpeak[1, "scmax"] <- -1 ##    ..
                      newpeak[1, "scale"] <- -1 ##    ..
                      
                      newpeak[1, "maxo"] <- max(gpeaks[jidx, "maxo"])
                      newpeak[1, "sn"]   <- max(gpeaks[jidx, "sn"])
                      newpeak[1, "lmin"] <- newmin
                      newpeak[1, "lmax"] <- newmax
                      newpeak[1, "rtmin"] <- scantime[td[newmin]]
                      newpeak[1, "rtmax"] <- scantime[td[newmax]]
                      newpeak[1, "rt"] <- weighted.mean(gpeaks[jidx, "rt"],
                                                        w = gpeaks[jidx, "maxo"])
                      
                      ## Re-assign m/z values
                      p1 <- match(td[newmin], otd)[1]
                      p2 <- match(td[newmax], otd)
                      p2 <- p2[length(p2)]
                      if (is.na(p1))
                        p1 <- 1
                      if (is.na(p2))
                        p2 <- length(omz)
                      mz.value <- omz[p1:p2]
                      mz.int <- od[p1:p2]
                      
                      ## re-calculate m/z value for peak range
                      mzmean <- do.call(mzCenterFun, list(mz = mz.value,
                                                          intensity = mz.int))
                      mzrange <- range(mz.value)
                      newpeak[1, "mz"] <- mzmean
                      newpeak[1, c("mzmin", "mzmax")] <- mzrange
                      
                      ## re-fit gaussian
                      md <- max(d[newmin:newmax])
                      d1 <- d[newmin:newmax] / md
                      pgauss <- fitGauss(td[newmin:newmax],
                                         d[newmin:newmax],
                                         pgauss = list(
                                           mu = td[newmin] +
                                             (td[newmax] - td[newmin]) /
                                             2,
                                           sigma = td[newmax] - td[newmin],
                                           h = max(gpeaks[jidx, "h"])
                                         ))
                      if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                        newpeak[1, "mu"]    <- pgauss$mu
                        newpeak[1, "sigma"] <- pgauss$sigma
                        newpeak[1, "h"]     <- pgauss$h
                        newpeak[1, "egauss"] <-
                          sqrt((1 / length(td[newmin:newmax])) *
                                 sum(((
                                   d1 - gauss(td[newmin:newmax],
                                              pgauss$h / md,
                                              pgauss$mu,
                                              pgauss$sigma)
                                 ) ^ 2)))
                      } else {
                        ## re-fit after join failed
                        newpeak[1, "mu"]       <- NA
                        newpeak[1, "sigma"]    <- NA
                        newpeak[1, "h"]        <- NA
                        newpeak[1, "egauss"]   <- NA
                      }
                      
                      newpeaks <- rbind(newpeaks, newpeak)
                    }
                    ## add the remaining peaks
                    jp <- unique(unlist(jlist))
                    if (dim(peaks)[1] - length(jp) > 0)
                      newpeaks <- rbind(newpeaks, gpeaks[-jp, ])
                    
                  } else
                    newpeaks <- gpeaks
                  
                  grt.min <- newpeaks[, "rtmin"]
                  grt.max <- newpeaks[, "rtmax"]
                  
                  if (nrow(peaks) - Ngp > 0) {
                    ## notgausspeaks
                    for (k in 1:nrow(notgausspeaks)) {
                      ## here we can only check if they are completely overlapped
                      ## by other peaks
                      if (!any((notgausspeaks[k, "rtmin"] >= grt.min) &
                               (notgausspeaks[k, "rtmax"] <= grt.max)))
                        newpeaks <- rbind(newpeaks, notgausspeaks[k,])
                    }
                  }
                  
                  rownames(newpeaks) <- NULL
                  newpeaks
                }
              peaks <- joinOverlappingPeaks(
                td,
                d,
                otd,
                omz,
                od,
                scantime,
                scan.range,
                peaks,
                maxGaussOverlap,
                mzCenterFun = mzCenterFun
              )
            }
            
            ## BEGIN - plotting/sleep
            if ((sleep > 0) && (!is.null(peaks))) {
              tdp <- scantime[td]
              trange <- range(tdp)
              egauss <-
                paste(round(peaks[, "egauss"], 3), collapse = ", ")
              cdppm <- paste(peaks[, "dppm"], collapse = ", ")
              csn <- paste(peaks[, "sn"], collapse = ", ")
              par(bg = "white")
              l <-
                layout(matrix(
                  c(1, 2, 3),
                  nrow = 3,
                  ncol = 1,
                  byrow = T
                ),
                heights = c(.5, .75, 2))
              
              par(mar = c(2, 4, 4, 2) + 0.1)
              ## plotRaw(object,mzrange=mzrange,rtrange=trange,log=TRUE,title='')
              ## Do plotRaw manually.
              raw_mat <- .rawMat(
                mz = mz,
                int = int,
                scantime = scantime,
                valsPerSpect = valsPerSpect,
                mzrange = mzrange,
                rtrange = rtrange,
                log = TRUE
              )
              if (nrow(raw_mat) > 0) {
                y <- raw_mat[, "intensity"]
                ylim <- range(y)
                y <- y / ylim[2]
                colorlut <- terrain.colors(16)
                col <- colorlut[y * 15 + 1]
                plot(
                  raw_mat[, "time"],
                  raw_mat[, "mz"],
                  pch = 20,
                  cex = .5,
                  main = "",
                  xlab = "Seconds",
                  ylab = "m/z",
                  col = col,
                  xlim = trange
                )
              } else {
                plot(
                  c(NA, NA),
                  main = "",
                  xlab = "Seconds",
                  ylab = "m/z",
                  xlim = trange,
                  ylim = mzrange
                )
              }
              ## done
              title(
                main = paste(
                  f,
                  ': ',
                  round(mzrange[1], 4),
                  ' - ',
                  round(mzrange[2], 4),
                  ' m/z , dppm=',
                  cdppm,
                  ', EGauss=',
                  egauss ,
                  ',  S/N =',
                  csn,
                  sep = ''
                )
              )
              par(mar = c(1, 4, 1, 2) + 0.1)
              image(
                y = scales[1:(dim(wCoefs)[2])],
                z = wCoefs,
                col = terrain.colors(256),
                xaxt = 'n',
                ylab = 'CWT coeff.'
              )
              par(mar = c(4, 4, 1, 2) + 0.1)
              plot(tdp, d, ylab = 'Intensity', xlab = 'Scan Time')
              lines(tdp, d, lty = 2)
              lines(scantime[otd], od, lty = 2, col = 'blue') ## original mzbox range
              abline(h = baseline, col = 'green')
              bwh <- length(sr[1]:sr[2]) - length(baseline)
              if (odd(bwh)) {
                bwh1 <-  floor(bwh / 2)
                bwh2 <- bwh1 + 1
              } else {
                bwh1 <- bwh2 <- bwh / 2
              }
              if (any(!is.na(peaks[, "scpos"])))
              {
                ## plot centers and width found through wavelet analysis
                abline(v = scantime[na.omit(peaks[(peaks[, "scpos"] > 0), "scpos"])], col =
                         'red')
              }
              abline(v = na.omit(c(peaks[, "rtmin"], peaks[, "rtmax"])),
                     col = 'green',
                     lwd = 1)
              if (fitgauss) {
                tdx <- seq(min(td), max(td), length.out = 200)
                tdxp <- seq(trange[1], trange[2], length.out = 200)
                fitted.peaks <- which(!is.na(peaks[, "mu"]))
                for (p in fitted.peaks)
                {
                  ## plot gaussian fits
                  yg <-
                    gauss(tdx, peaks[p, "h"], peaks[p, "mu"], peaks[p, "sigma"])
                  lines(tdxp, yg, col = 'blue')
                }
              }
              Sys.sleep(sleep)
            }
            ## -- END plotting/sleep
            
            if (!is.null(peaks)) {
              peaklist[[length(peaklist) + 1]] <- peaks
            }
          } ## f
          
          if (length(peaklist) == 0) {
            warning("No peaks found!")
            
            if (verboseColumns) {
              nopeaks <- matrix(nrow = 0,
                                ncol = length(basenames) +
                                  length(verbosenames))
              colnames(nopeaks) <- c(basenames, verbosenames)
            } else {
              nopeaks <- matrix(nrow = 0, ncol = length(basenames))
              colnames(nopeaks) <- c(basenames)
            }
            message(" FAIL: none found!")
            return(nopeaks)
          }
          p <- do.call(rbind, peaklist)
          if (!verboseColumns)
            p <- p[, basenames, drop = FALSE]
          
          uorder <- order(p[, "into"], decreasing = TRUE)
          pm <-
            as.matrix(p[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
          rectUnique <-
            function(m,
                     order = seq(length = nrow(m)),
                     xdiff = 0,
                     ydiff = 0) {
              nr <- nrow(m)
              nc <- ncol(m)
              
              if (is.null(nr) || nr == 0) {
                ## empty matrix in first place
                return (m)
              }
              
              if (!is.double(m))
                m <- as.double(m)
              .C(
                "RectUnique",
                m,
                as.integer(order - 1),
                nr,
                nc,
                as.double(xdiff),
                as.double(ydiff),
                logical(nrow(m)),
                PACKAGE = "xcms"
              )[[7]]
            }
          uindex <-
            rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) ## allow adjacent peaks
          pr <- p[uindex, , drop = FALSE]
          message(" OK: ", nrow(pr), " found.")
          
          return(pr)
        }
        
        if (getOption("originalCentWave", default = TRUE)) {
          ## message("DEBUG: using original centWave.")
          centWave_orig_sxt(
            mz = mz,
            int = int,
            scantime = scantime,
            valsPerSpect = valsPerSpect,
            ppm = ppm,
            peakwidth = peakwidth,
            snthresh = snthresh,
            prefilter = prefilter,
            mzCenterFun = mzCenterFun,
            integrate = integrate,
            mzdiff = mzdiff,
            fitgauss = fitgauss,
            noise = noise,
            verboseColumns = verboseColumns,
            roiList = roiList,
            firstBaselineCheck = firstBaselineCheck,
            roiScales = roiScales,
            sleep = sleep,
            extendLengthMSW = extendLengthMSW,
            file_name = file_name
          )
        } else {
          ## message("DEBUG: using modified centWave.")
          centWave_orig_sxt(
            mz = mz,
            int = int,
            scantime = scantime,
            valsPerSpect = valsPerSpect,
            ppm = ppm,
            peakwidth = peakwidth,
            snthresh = snthresh,
            prefilter = prefilter,
            mzCenterFun = mzCenterFun,
            integrate = integrate,
            mzdiff = mzdiff,
            fitgauss = fitgauss,
            noise = noise,
            verboseColumns = verboseColumns,
            roiList = roiList,
            firstBaselineCheck = firstBaselineCheck,
            roiScales = roiScales,
            sleep = sleep,
            file_name = file_name
          )
        }
      }
    
    res <- do.call(method, args = c(
      list(
        mz = unlist(mzs, use.names = FALSE),
        int = unlist(lapply(x, intensity),
                     use.names = FALSE),
        valsPerSpect = vals_per_spect,
        scantime = rt,
        file_name = file_name
      ),
      as(param, "list")
    ))
    ## Ensure that we call the garbage collector to eventually clean unused stuff
    rm(mzs)
    rm(x)
    rm(rt)
    gc()
    list(peaks = res, date = procDat)
  }