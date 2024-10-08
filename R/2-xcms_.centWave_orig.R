# ###remove dot
# centWave_orig_sxt <- function(mz,
#                                int,
#                                scantime,
#                                valsPerSpect,
#                                ppm = 25,
#                                peakwidth = c(20, 50),
#                                snthresh = 10,
#                                prefilter = c(3, 100),
#                                mzCenterFun = "wMean",
#                                integrate = 1,
#                                mzdiff = -0.001,
#                                fitgauss = FALSE,
#                                noise = 0,
#                                ## noise.local=TRUE,
#                                sleep = 0,
#                                verboseColumns = FALSE,
#                                roiList = list(),
#                                firstBaselineCheck = TRUE,
#                                roiScales = NULL,
#                                extendLengthMSW = FALSE,
#                                file_name = "file") {
#   ## Input argument checking.
#   if (missing(mz) |
#       missing(int) | missing(scantime) | missing(valsPerSpect))
#     stop("Arguments 'mz', 'int', 'scantime' and 'valsPerSpect'",
#          " are required!")
#   if (length(mz) != length(int) |
#       length(valsPerSpect) != length(scantime)
#       | length(mz) != sum(valsPerSpect))
#     stop(
#       "Lengths of 'mz', 'int' and of 'scantime','valsPerSpect'",
#       " have to match. Also, 'length(mz)' should be equal to",
#       " 'sum(valsPerSpect)'."
#     )
#   scanindex <-
#     valueCount2ScanIndex(valsPerSpect) ## Get index vector for C calls
#   if (!is.double(mz))
#     mz <- as.double(mz)
#   if (!is.double(int))
#     int <- as.double(int)
#   ## Fix the mzCenterFun
#   mzCenterFun <- paste(
#     "mzCenter",
#     gsub(
#       mzCenterFun,
#       pattern = "mzCenter.",
#       replacement = "",
#       fixed = TRUE
#     ),
#     sep = "."
#   )
#   if (!exists(mzCenterFun, mode = "function"))
#     stop("Function '", mzCenterFun, "' not defined !")
#   
#   if (!is.logical(firstBaselineCheck))
#     stop("Parameter 'firstBaselineCheck' should be logical!")
#   if (length(firstBaselineCheck) != 1)
#     stop("Parameter 'firstBaselineCheck' should be a single logical !")
#   if (length(roiScales) > 0)
#     if (length(roiScales) != length(roiList) |
#         !is.numeric(roiScales))
#       stop(
#         "If provided, parameter 'roiScales' has to be a numeric with",
#         " length equal to the length of 'roiList'!"
#       )
#   ## if (!is.null(roiScales)) {
#   ##     if (!is.numeric(roiScales) | length(roiScales) != length(roiList))
#   ##         stop("Parameter 'roiScales' has to be a numeric of length equal to",
#   ##              " parameter 'roiList'!")
#   ##}
#   
#   basenames <- c("mz",
#                  "mzmin",
#                  "mzmax",
#                  "rt",
#                  "rtmin",
#                  "rtmax",
#                  "into",
#                  "intb",
#                  "maxo",
#                  "sn")
#   verbosenames <-
#     c(
#       "egauss",
#       "mu",
#       "sigma",
#       "h",
#       "f",
#       "dppm",
#       "scale",
#       "scpos",
#       "scmin",
#       "scmax",
#       "lmin",
#       "lmax"
#     )
#   
#   ## Peak width: seconds to scales
#   scalerange <- round((peakwidth / mean(diff(scantime))) / 2)
#   
#   if (length(z <- which(scalerange == 0)))
#     scalerange <- scalerange[-z]
#   if (length(scalerange) < 1) {
#     warning("No scales? Please check peak width!")
#     if (verboseColumns) {
#       nopeaks <- matrix(nrow = 0,
#                         ncol = length(basenames) +
#                           length(verbosenames))
#       colnames(nopeaks) <- c(basenames, verbosenames)
#     } else {
#       nopeaks <- matrix(nrow = 0, ncol = length(basenames))
#       colnames(nopeaks) <- c(basenames)
#     }
#     return(invisible(nopeaks))
#   }
#   
#   if (length(scalerange) > 1) {
#     scales <- seq(from = scalerange[1],
#                   to = scalerange[2],
#                   by = 2)
#   } else{
#     scales <- scalerange
#   }
#   
#   minPeakWidth <-  scales[1]
#   noiserange <- c(minPeakWidth * 3, max(scales) * 3)
#   maxGaussOverlap <- 0.5
#   minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
#   minCentroids <- minPtsAboveBaseLine
#   scRangeTol <-  maxDescOutlier <- floor(minPeakWidth / 2)
#   scanrange <- c(1, length(scantime))
#   
#   ## If no ROIs are supplied then search for them.
#   if (length(roiList) == 0) {
#     message("Detecting mass traces at ", ppm, " ppm ... ", file_name, appendLF = FALSE)
#     ## flush.console();
#     ## We're including the findmzROI code in this function to reduce
#     ## the need to copy objects etc.
#     ## We could also sort the data by m/z anyway; wouldn't need that
#     ## much time. Once we're using classes from MSnbase we can be
#     ## sure that values are correctly sorted.
#     withRestarts(
#       tryCatch({
#         tmp <- capture.output(
#           roiList <- .Call(
#             "findmzROI",
#             mz,
#             int,
#             scanindex,
#             as.double(c(0.0, 0.0)),
#             as.integer(scanrange),
#             as.integer(length(scantime)),
#             as.double(ppm * 1e-6),
#             as.integer(minCentroids),
#             as.integer(prefilter),
#             as.integer(noise),
#             PACKAGE = 'xcms'
#           )
#         )
#       },
#       error = function(e) {
#         if (grepl("m/z sort assumption violated !", e$message)) {
#           invokeRestart("fixSort")
#         } else {
#           simpleError(e)
#         }
#       }),
#       fixSort = function() {
#         ## Force ordering of values within spectrum by mz:
#         ##  o split values into a list -> mz per spectrum, intensity per
#         ##    spectrum.
#         ##  o define the ordering.
#         ##  o re-order the mz and intensity and unlist again.
#         ## Note: the Rle split is faster than the "conventional" factor split.
#         splitF <- Rle(1:length(valsPerSpect), valsPerSpect)
#         mzl <- as.list(S4Vectors::split(mz, f = splitF))
#         oidx <- lapply(mzl, order)
#         mz <<- unlist(
#           mapply(
#             mzl,
#             oidx,
#             FUN = function(y, z) {
#               return(y[z])
#             },
#             SIMPLIFY = FALSE,
#             USE.NAMES = FALSE
#           ),
#           use.names = FALSE
#         )
#         int <<- unlist(
#           mapply(
#             as.list(split(int, f = splitF)),
#             oidx,
#             FUN = function(y, z) {
#               return(y[z])
#             },
#             SIMPLIFY = FALSE,
#             USE.NAMES = FALSE
#           ),
#           use.names = FALSE
#         )
#         rm(mzl)
#         rm(splitF)
#         tmp <- capture.output(
#           roiList <<- .Call(
#             "findmzROI",
#             mz,
#             int,
#             scanindex,
#             as.double(c(0.0, 0.0)),
#             as.integer(scanrange),
#             as.integer(length(scantime)),
#             as.double(ppm * 1e-6),
#             as.integer(minCentroids),
#             as.integer(prefilter),
#             as.integer(noise),
#             PACKAGE = 'xcms'
#           )
#         )
#       }
#     )
#     message("OK")
#     ## ROI.list <- findmzROI(object,scanrange=scanrange,dev=ppm * 1e-6,minCentroids=minCentroids, prefilter=prefilter, noise=noise)
#     if (length(roiList) == 0) {
#       warning("No ROIs found! \n")
#       if (verboseColumns) {
#         nopeaks <- matrix(nrow = 0,
#                           ncol = length(basenames) +
#                             length(verbosenames))
#         colnames(nopeaks) <- c(basenames, verbosenames)
#       } else {
#         nopeaks <- matrix(nrow = 0, ncol = length(basenames))
#         colnames(nopeaks) <- c(basenames)
#       }
#       return(invisible(nopeaks))
#     }
#   }
#   
#   ## Second stage: process the ROIs
#   peaklist <- list()
#   Nscantime <- length(scantime)
#   lf <- length(roiList)
#   
#   ## cat('\n Detecting chromatographic peaks ... \n % finished: ')
#   ## lp <- -1
#   message(
#     "Detecting chromatographic peaks in ",
#     length(roiList),
#     " regions of interest ...",
#     appendLF = FALSE
#   )
#   
#   for (f in  1:lf) {
#     ## ## Show progress
#     ## perc <- round((f/lf) * 100)
#     ## if ((perc %% 10 == 0) && (perc != lp))
#     ## {
#     ##     cat(perc," ",sep="");
#     ##     lp <- perc;
#     ## }
#     ## flush.console()
#     
#     feat <- roiList[[f]]
#     N <- feat$scmax - feat$scmin + 1
#     peaks <- peakinfo <- NULL
#     mzrange <- c(feat$mzmin, feat$mzmax)
#     sccenter <- feat$scmin[1] + floor(N / 2) - 1
#     scrange <- c(feat$scmin, feat$scmax)
#     ## scrange + noiserange, used for baseline detection and wavelet analysis
#     sr <- c(max(scanrange[1], scrange[1] - max(noiserange)),
#             min(scanrange[2], scrange[2] + max(noiserange)))
#     eic <- .Call(
#       "getEIC",
#       mz,
#       int,
#       scanindex,
#       as.double(mzrange),
#       as.integer(sr),
#       as.integer(length(scanindex)),
#       PACKAGE = "xcms"
#     )
#     ## eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
#     d <- eic$intensity
#     td <- sr[1]:sr[2]
#     scan.range <- c(sr[1], sr[2])
#     ## original mzROI range
#     idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
#     mzROI.EIC <-
#       list(scan = eic$scan[idxs], intensity = eic$intensity[idxs])
#     ## mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
#     omz <-
#       .Call(
#         "getWeightedMZ",
#         mz,
#         int,
#         scanindex,
#         as.double(mzrange),
#         as.integer(scrange),
#         as.integer(length(scantime)),
#         PACKAGE = 'xcms'
#       )
#     ## omz <- rawMZ(object,mzrange=mzrange,scanrange=scrange)
#     if (all(omz == 0)) {
#       warning("centWave: no peaks found in ROI.")
#       next
#     }
#     od  <- mzROI.EIC$intensity
#     otd <- mzROI.EIC$scan
#     if (all(od == 0)) {
#       warning("centWave: no peaks found in ROI.")
#       next
#     }
#     
#     ## scrange + scRangeTol, used for gauss fitting and continuous
#     ## data above 1st baseline detection
#     ftd <- max(td[1], scrange[1] - scRangeTol):min(td[length(td)],
#                                                    scrange[2] + scRangeTol)
#     fd <- d[match(ftd, td)]
#     
#     ## 1st type of baseline: statistic approach
#     if (N >= 10 * minPeakWidth) {
#       ## in case of very long mass trace use full scan range
#       ## for baseline detection
#       noised <-
#         .Call(
#           "getEIC",
#           mz,
#           int,
#           scanindex,
#           as.double(mzrange),
#           as.integer(scanrange),
#           as.integer(length(scanindex)),
#           PACKAGE = "xcms"
#         )$intensity
#       ## noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity
#     } else {
#       noised <- d
#     }
#     ## 90% trimmed mean as first baseline guess
#     noise <- estimateChromNoise(noised, trim = 0.05,
#                                 minPts = 3 * minPeakWidth)
#     ## any continuous data above 1st baseline ?
#     if (firstBaselineCheck &&
#         !continuousPtsAboveThreshold(fd, threshold = noise,
#                                      num = minPtsAboveBaseLine))
#       next
#     ## 2nd baseline estimate using not-peak-range
#     lnoise <-
#       getLocalNoiseEstimate(d,
#                             td,
#                             ftd,
#                             noiserange,
#                             Nscantime,
#                             threshold = noise,
#                             num = minPtsAboveBaseLine)
#     ## Final baseline & Noise estimate
#     baseline <- max(1, min(lnoise[1], noise))
#     sdnoise <- max(1, lnoise[2])
#     sdthr <-  sdnoise * snthresh
#     ## is there any data above S/N * threshold ?
#     if (!(any(fd - baseline >= sdthr)))
#       next
#     wCoefs <- MSW.cwt(
#       d,
#       scales = scales,
#       wavelet = 'mexh',
#       extendLengthMSW = extendLengthMSW
#     )
#     if (!(!is.null(dim(wCoefs)) && any(wCoefs - baseline >= sdthr)))
#       next
#     if (td[length(td)] == Nscantime)
#       ## workaround, localMax fails otherwise
#       wCoefs[nrow(wCoefs),] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
#     localMax <- MSW.getLocalMaximumCWT(wCoefs)
#     rL <- MSW.getRidge(localMax)
#     wpeaks <- sapply(rL,
#                      function(x) {
#                        w <- min(1:length(x), ncol(wCoefs))
#                        any(wCoefs[x, w] - baseline >= sdthr)
#                      })
#     if (any(wpeaks)) {
#       wpeaksidx <- which(wpeaks)
#       ## check each peak in ridgeList
#       for (p in 1:length(wpeaksidx)) {
#         opp <- rL[[wpeaksidx[p]]]
#         pp <- unique(opp)
#         if (length(pp) >= 1) {
#           dv <- td[pp] %in% ftd
#           if (any(dv)) {
#             ## peaks in orig. data range
#             ## Final S/N check
#             if (any(d[pp[dv]] - baseline >= sdthr)) {
#               ## if(!is.null(roiScales)) {
#               ## allow roiScales to be a numeric of length 0
#               if (length(roiScales) > 0) {
#                 ## use given scale
#                 best.scale.nr <- which(scales == roiScales[[f]])
#                 if (best.scale.nr > length(opp))
#                   best.scale.nr <- length(opp)
#               } else {
#                 ## try to decide which scale describes the peak best
#                 inti <- numeric(length(opp))
#                 irange <- rep(ceiling(scales[1] / 2), length(opp))
#                 for (k in 1:length(opp)) {
#                   kpos <- opp[k]
#                   r1 <- ifelse(kpos - irange[k] > 1,
#                                kpos - irange[k], 1)
#                   r2 <- ifelse(kpos + irange[k] < length(d),
#                                kpos + irange[k],
#                                length(d))
#                   inti[k] <- sum(d[r1:r2])
#                 }
#                 maxpi <- which.max(inti)
#                 if (length(maxpi) > 1) {
#                   m <- wCoefs[opp[maxpi], maxpi]
#                   bestcol <- which(m == max(m),
#                                    arr.ind = TRUE)[2]
#                   best.scale.nr <- maxpi[bestcol]
#                 } else
#                   best.scale.nr <- maxpi
#               }
#               
#               best.scale <-  scales[best.scale.nr]
#               best.scale.pos <- opp[best.scale.nr]
#               
#               pprange <- min(pp):max(pp)
#               ## maxint <- max(d[pprange])
#               lwpos <- max(1, best.scale.pos - best.scale)
#               rwpos <- min(best.scale.pos + best.scale, length(td))
#               p1 <- match(td[lwpos], otd)[1]
#               p2 <- match(td[rwpos], otd)
#               p2 <- p2[length(p2)]
#               if (is.na(p1))
#                 p1 <- 1
#               if (is.na(p2))
#                 p2 <- N
#               mz.value <- omz[p1:p2]
#               mz.int <- od[p1:p2]
#               maxint <- max(mz.int)
#               
#               ## re-calculate m/z value for peak range
#               mzrange <- range(mz.value)
#               mzmean <- do.call(mzCenterFun,
#                                 list(mz = mz.value,
#                                      intensity = mz.int))
#               
#               ## Compute dppm only if needed
#               dppm <- NA
#               if (verboseColumns) {
#                 if (length(mz.value) >= (minCentroids + 1)) {
#                   dppm <- round(min(
#                     running(
#                       abs(diff(mz.value)) /
#                         (mzrange[2] *  1e-6),
#                       fun = max,
#                       width = minCentroids
#                     )
#                   ))
#                 } else {
#                   dppm <- round((mzrange[2] - mzrange[1]) /
#                                   (mzrange[2] * 1e-6))
#                 }
#               }
#               peaks <- rbind(
#                 peaks,
#                 c(
#                   mzmean,
#                   mzrange,
#                   ## mz
#                   NA,
#                   NA,
#                   NA,
#                   ## rt, rtmin, rtmax,
#                   NA,
#                   ## intensity (sum)
#                   NA,
#                   ## intensity (-bl)
#                   maxint,
#                   ## max intensity
#                   round((maxint - baseline) / sdnoise),
#                   ##  S/N Ratio
#                   NA,
#                   ## Gaussian RMSE
#                   NA,
#                   NA,
#                   NA,
#                   ## Gaussian Parameters
#                   f,
#                   ## ROI Position
#                   dppm,
#                   ## max. difference between the [minCentroids] peaks in ppm
#                   best.scale,
#                   ## Scale
#                   td[best.scale.pos],
#                   td[lwpos],
#                   td[rwpos],
#                   ## Peak positions guessed from the wavelet's (scan nr)
#                   NA,
#                   NA
#                 )
#               )                    ## Peak limits (scan nr)
#               peakinfo <- rbind(peakinfo,
#                                 c(
#                                   best.scale,
#                                   best.scale.nr,
#                                   best.scale.pos,
#                                   lwpos,
#                                   rwpos
#                                 ))
#               ## Peak positions guessed from the wavelet's
#             }
#           }
#         }
#       }  ##for
#     } ## if
#     
#     ##  postprocessing
#     if (!is.null(peaks)) {
#       colnames(peaks) <- c(basenames, verbosenames)
#       colnames(peakinfo) <- c("scale", "scaleNr", "scpos",
#                               "scmin", "scmax")
#       for (p in 1:dim(peaks)[1]) {
#         ## find minima (peak boundaries), assign rt and intensity values
#         if (integrate == 1) {
#           lm <- descendMin(wCoefs[, peakinfo[p, "scaleNr"]],
#                            istart = peakinfo[p, "scpos"])
#           gap <-
#             all(d[lm[1]:lm[2]] == 0) # looks like we got stuck in a gap right in the middle of the peak
#           if ((lm[1] == lm[2]) || gap)
#             # fall-back
#             lm <-
#             descendMinTol(d, startpos = c(peakinfo[p, "scmin"],
#                                           peakinfo[p, "scmax"]),
#                           maxDescOutlier)
#         } else {
#           lm <- descendMinTol(d, startpos = c(peakinfo[p, "scmin"],
#                                               peakinfo[p, "scmax"]),
#                               maxDescOutlier)
#         }
#         ## Narrow peak rt boundaries by removing values below threshold
#         lm <- .narrow_rt_boundaries(lm, d)
#         lm_seq <- lm[1]:lm[2]
#         pd <- d[lm_seq]
#         
#         peakrange <- td[lm]
#         peaks[p, "rtmin"] <- scantime[peakrange[1]]
#         peaks[p, "rtmax"] <- scantime[peakrange[2]]
#         peaks[p, "maxo"] <- max(pd)
#         pwid <- (scantime[peakrange[2]] - scantime[peakrange[1]]) /
#           (peakrange[2] - peakrange[1])
#         if (is.na(pwid))
#           pwid <- 1
#         peaks[p, "into"] <- pwid * sum(pd)
#         db <- pd - baseline
#         peaks[p, "intb"] <- pwid * sum(db[db > 0])
#         peaks[p, "lmin"] <- lm[1]
#         peaks[p, "lmax"] <- lm[2]
#         
#         if (fitgauss) {
#           ## perform gaussian fits, use wavelets for inital parameters
#           td_lm <- td[lm_seq]
#           md <- max(pd)
#           d1 <- pd / md ## normalize data for gaussian error calc.
#           pgauss <- fitGauss(td_lm,
#                              pd,
#                              pgauss = list(
#                                mu = peaks[p, "scpos"],
#                                sigma = peaks[p, "scmax"] -
#                                  peaks[p, "scmin"],
#                                h = peaks[p, "maxo"]
#                              ))
#           rtime <- peaks[p, "scpos"]
#           if (!any(is.na(pgauss)) && all(pgauss > 0)) {
#             gtime <- td[match(round(pgauss$mu), td)]
#             if (!is.na(gtime)) {
#               rtime <- gtime
#               peaks[p, "mu"] <- pgauss$mu
#               peaks[p, "sigma"] <- pgauss$sigma
#               peaks[p, "h"] <- pgauss$h
#               peaks[p, "egauss"] <- sqrt((1 / length(td_lm)) *
#                                            sum(((
#                                              d1 - gauss(td_lm, pgauss$h / md,
#                                                         pgauss$mu, pgauss$sigma)
#                                            ) ^ 2)))
#             }
#           }
#           peaks[p, "rt"] <- scantime[rtime]
#           ## avoid fitting side effects
#           if (peaks[p, "rt"] < peaks[p, "rtmin"])
#             peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
#         } else
#           peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
#       }
#       peaks <- joinOverlappingPeaks(td,
#                                     d,
#                                     otd,
#                                     omz,
#                                     od,
#                                     scantime,
#                                     scan.range,
#                                     peaks,
#                                     maxGaussOverlap,
#                                     mzCenterFun = mzCenterFun)
#     }
#     
#     ## BEGIN - plotting/sleep
#     if ((sleep > 0) && (!is.null(peaks))) {
#       tdp <- scantime[td]
#       trange <- range(tdp)
#       egauss <- paste(round(peaks[, "egauss"], 3), collapse = ", ")
#       cdppm <- paste(peaks[, "dppm"], collapse = ", ")
#       csn <- paste(peaks[, "sn"], collapse = ", ")
#       par(bg = "white")
#       l <-
#         layout(matrix(
#           c(1, 2, 3),
#           nrow = 3,
#           ncol = 1,
#           byrow = T
#         ), heights = c(.5, .75, 2))
#       
#       par(mar = c(2, 4, 4, 2) + 0.1)
#       ## plotRaw(object,mzrange=mzrange,rtrange=trange,log=TRUE,title='')
#       ## Do plotRaw manually.
#       raw_mat <- .rawMat(
#         mz = mz,
#         int = int,
#         scantime = scantime,
#         valsPerSpect = valsPerSpect,
#         mzrange = mzrange,
#         rtrange = rtrange,
#         log = TRUE
#       )
#       if (nrow(raw_mat) > 0) {
#         y <- raw_mat[, "intensity"]
#         ylim <- range(y)
#         y <- y / ylim[2]
#         colorlut <- terrain.colors(16)
#         col <- colorlut[y * 15 + 1]
#         plot(
#           raw_mat[, "time"],
#           raw_mat[, "mz"],
#           pch = 20,
#           cex = .5,
#           main = "",
#           xlab = "Seconds",
#           ylab = "m/z",
#           col = col,
#           xlim = trange
#         )
#       } else {
#         plot(
#           c(NA, NA),
#           main = "",
#           xlab = "Seconds",
#           ylab = "m/z",
#           xlim = trange,
#           ylim = mzrange
#         )
#       }
#       ## done
#       title(
#         main = paste(
#           f,
#           ': ',
#           round(mzrange[1], 4),
#           ' - ',
#           round(mzrange[2], 4),
#           ' m/z , dppm=',
#           cdppm,
#           ', EGauss=',
#           egauss ,
#           ',  S/N =',
#           csn,
#           sep = ''
#         )
#       )
#       par(mar = c(1, 4, 1, 2) + 0.1)
#       image(
#         y = scales[1:(dim(wCoefs)[2])],
#         z = wCoefs,
#         col = terrain.colors(256),
#         xaxt = 'n',
#         ylab = 'CWT coeff.'
#       )
#       par(mar = c(4, 4, 1, 2) + 0.1)
#       plot(tdp, d, ylab = 'Intensity', xlab = 'Scan Time')
#       lines(tdp, d, lty = 2)
#       lines(scantime[otd], od, lty = 2, col = 'blue') ## original mzbox range
#       abline(h = baseline, col = 'green')
#       bwh <- length(sr[1]:sr[2]) - length(baseline)
#       if (odd(bwh)) {
#         bwh1 <-  floor(bwh / 2)
#         bwh2 <- bwh1 + 1
#       } else {
#         bwh1 <- bwh2 <- bwh / 2
#       }
#       if (any(!is.na(peaks[, "scpos"])))
#       {
#         ## plot centers and width found through wavelet analysis
#         abline(v = scantime[na.omit(peaks[(peaks[, "scpos"] > 0), "scpos"])], col =
#                  'red')
#       }
#       abline(v = na.omit(c(peaks[, "rtmin"], peaks[, "rtmax"])),
#              col = 'green',
#              lwd = 1)
#       if (fitgauss) {
#         tdx <- seq(min(td), max(td), length.out = 200)
#         tdxp <- seq(trange[1], trange[2], length.out = 200)
#         fitted.peaks <- which(!is.na(peaks[, "mu"]))
#         for (p in fitted.peaks)
#         {
#           ## plot gaussian fits
#           yg <-
#             gauss(tdx, peaks[p, "h"], peaks[p, "mu"], peaks[p, "sigma"])
#           lines(tdxp, yg, col = 'blue')
#         }
#       }
#       Sys.sleep(sleep)
#     }
#     ## -- END plotting/sleep
#     
#     if (!is.null(peaks)) {
#       peaklist[[length(peaklist) + 1]] <- peaks
#     }
#   } ## f
#   
#   if (length(peaklist) == 0) {
#     warning("No peaks found!")
#     
#     if (verboseColumns) {
#       nopeaks <- matrix(nrow = 0,
#                         ncol = length(basenames) +
#                           length(verbosenames))
#       colnames(nopeaks) <- c(basenames, verbosenames)
#     } else {
#       nopeaks <- matrix(nrow = 0, ncol = length(basenames))
#       colnames(nopeaks) <- c(basenames)
#     }
#     message(" FAIL: none found!")
#     return(nopeaks)
#   }
#   p <- do.call(rbind, peaklist)
#   if (!verboseColumns)
#     p <- p[, basenames, drop = FALSE]
#   
#   uorder <- order(p[, "into"], decreasing = TRUE)
#   pm <-
#     as.matrix(p[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE])
#   uindex <-
#     rectUnique(pm, uorder, mzdiff, ydiff = -0.00001) ## allow adjacent peaks
#   pr <- p[uindex, , drop = FALSE]
#   message(" OK: ", nrow(pr), " found.")
#   
#   return(pr)
# }
# 
