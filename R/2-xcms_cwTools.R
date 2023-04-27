# MSW.cwt <-
#   function (ms,
#             scales = 1,
#             wavelet = "mexh",
#             extendLengthMSW = FALSE) {
#     if (wavelet == "mexh") {
#       psi_xval <- seq(-6, 6, length = 256)
#       psi <- (2 / sqrt(3) * pi ^ (-0.25)) * (1 - psi_xval ^ 2) *
#         exp(-psi_xval ^ 2 / 2)
#     }
#     else if (is.matrix(wavelet)) {
#       if (nrow(wavelet) == 2) {
#         psi_xval <- wavelet[1, ]
#         psi <- wavelet[2, ]
#       }
#       else if (ncol(wavelet) == 2) {
#         psi_xval <- wavelet[, 1]
#         psi <- wavelet[, 2]
#       }
#       else {
#         stop("Unsupported wavelet format!")
#       }
#     }
#     else {
#       stop("Unsupported wavelet!")
#     }
#     oldLen <- length(ms)
#     # IF extendLengthMSW is TRUE:
#     # The new length is determined by the scales argument, so a larger peakwidth
#     # will ensure more scales are run, but may slow it down. See
#     # https://github.com/sneumann/xcms/issues/445 for more information about
#     # a change from using extendNBase to extendLength.
#     if (extendLengthMSW) {
#       newLen <- 2 ^ (ceiling(log2(max(scales) * 12)))
#       ms <-
#         MSW.extendLength(
#           x = ms,
#           addLength = (newLen - length(ms)),
#           method = "open"
#         )
#     } else {
#       ms <- MSW.extendNBase(ms, nLevel = NULL, base = 2)
#     }
#     
#     
#     len <- length(ms)
#     nbscales <- length(scales)
#     wCoefs <- NULL
#     psi_xval <- psi_xval - psi_xval[1]
#     dxval <- psi_xval[2]
#     xmax <- psi_xval[length(psi_xval)]
#     for (i in 1:length(scales)) {
#       scale.i <- scales[i]
#       f <- rep(0, len)
#       j <- 1 + floor((0:(scale.i * xmax)) / (scale.i * dxval))
#       if (length(j) == 1)
#         j <- c(1, 1)
#       lenWave <- length(j)
#       f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
#       if (length(f) > len)
#       {
#         i <- i - 1
#         break
#       }   ##  stop(paste("scale", scale.i, "is too large!"))
#       wCoefs.i <- 1 / sqrt(scale.i) * convolve(ms, f)
#       wCoefs.i <- c(wCoefs.i[(len - floor(lenWave / 2) + 1):len],
#                     wCoefs.i[1:(len - floor(lenWave / 2))])
#       wCoefs <- cbind(wCoefs, wCoefs.i)
#     }
#     if (i < 1)
#       return(NA)
#     scales <- scales[1:i]
#     if (length(scales) == 1)
#       wCoefs <- matrix(wCoefs, ncol = 1)
#     colnames(wCoefs) <- scales
#     wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
#     wCoefs
#   }

# This function is no longer used by MSW.cwt(): see above note about the
# switch from extendNBase to calling extendLength directly.
# Possibly now unecessary?
# MSW.extendNBase <- function(x,
#                             nLevel = 1,
#                             base = 2,
#                             ...)
# {
#   ## from package MassSpecWavelet
#   if (!is.matrix(x))
#     x <- matrix(x, ncol = 1)
#   
#   nR <- nrow(x)
#   if (is.null(nLevel)) {
#     nR1 <- nextn(nR, base)
#   } else {
#     nR1 <- ceiling(nR / base ^ nLevel) * base ^ nLevel
#   }
#   if (nR != nR1) {
#     x <- MSW.extendLength(x, addLength = nR1 - nR, ...)
#   }
#   x
# }

# MSW.extendLength <-
#   function(x,
#            addLength = NULL,
#            method = c('reflection', 'open', 'circular'),
#            direction = c('right', 'left', 'both'))
#   {
#     ## from package MassSpecWavelet
#     if (is.null(addLength))
#       stop('Please provide the length to be added!')
#     if (!is.matrix(x))
#       x <- matrix(x, ncol = 1)
#     method <- match.arg(method)
#     direction <- match.arg(direction)
#     
#     nR <- nrow(x)
#     nR1 <- nR + addLength
#     if (direction == 'both') {
#       left <- right <- addLength
#     } else if (direction == 'right') {
#       left <- 0
#       right <- addLength
#     } else if (direction == 'left') {
#       left <- addLength
#       right <- 0
#     }
#     
#     if (right > 0) {
#       x <- switch(
#         method,
#         reflection = rbind(x, x[nR:(2 * nR - nR1 + 1), , drop =
#                                   FALSE]),
#         open = rbind(x, matrix(
#           rep(x[nR,], addLength), ncol = ncol(x), byrow = TRUE
#         )),
#         circular = rbind(x, x[1:(nR1 - nR), , drop = FALSE])
#       )
#     }
#     
#     if (left > 0) {
#       x <- switch(
#         method,
#         reflection = rbind(x[addLength:1, , drop = FALSE], x),
#         open = rbind(matrix(
#           rep(x[1,], addLength), ncol = ncol(x), byrow = TRUE
#         ), x),
#         circular = rbind(x[(2 * nR - nR1 + 1):nR, , drop = FALSE], x)
#       )
#     }
#     if (ncol(x) == 1)
#       x <- as.vector(x)
#     
#     x
#   }

# MSW.getLocalMaximumCWT <-
#   function(wCoefs,
#            minWinSize = 5,
#            amp.Th = 0)
#   {
#     ## from package MassSpecWavelet
#     localMax <- NULL
#     scales <- as.numeric(colnames(wCoefs))
#     
#     for (i in 1:length(scales)) {
#       scale.i <- scales[i]
#       winSize.i <- scale.i * 2 + 1
#       if (winSize.i < minWinSize) {
#         winSize.i <- minWinSize
#       }
#       temp <- MSW.localMaximum(wCoefs[, i], winSize.i)
#       localMax <- cbind(localMax, temp)
#     }
#     ## Set the values less than peak threshold as 0
#     localMax[wCoefs < amp.Th] <- 0
#     colnames(localMax) <- colnames(wCoefs)
#     rownames(localMax) <- rownames(wCoefs)
#     localMax
#   }

# MSW.localMaximum <-
#   function (x, winSize = 5)
#   {
#     ## from package MassSpecWavelet
#     len <- length(x)
#     rNum <- ceiling(len / winSize)
#     
#     ## Transform the vector as a matrix with column length equals winSize
#     ##		and find the maximum position at each row.
#     y <-
#       matrix(c(x, rep(x[len], rNum * winSize - len)), nrow = winSize)
#     y.maxInd <- apply(y, 2, which.max)
#     ## Only keep the maximum value larger than the boundary values
#     selInd <-
#       which(apply(y, 2, function(x)
#         max(x) > x[1] & max(x) > x[winSize]))
#     
#     ## keep the result
#     localMax <- rep(0, len)
#     localMax[(selInd - 1) * winSize + y.maxInd[selInd]] <- 1
#     
#     ## Shift the vector with winSize/2 and do the same operation
#     shift <- floor(winSize / 2)
#     rNum <- ceiling((len + shift) / winSize)
#     y <-
#       matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow =
#                winSize)
#     y.maxInd <- apply(y, 2, which.max)
#     ## Only keep the maximum value larger than the boundary values
#     selInd <-
#       which(apply(y, 2, function(x)
#         max(x) > x[1] & max(x) > x[winSize]))
#     localMax[(selInd - 1) * winSize + y.maxInd[selInd] - shift] <- 1
#     
#     ## Check whether there is some local maxima have in between distance less than winSize
#     maxInd <- which(localMax > 0)
#     selInd <- which(diff(maxInd) < winSize)
#     if (length(selInd) > 0) {
#       selMaxInd1 <- maxInd[selInd]
#       selMaxInd2 <- maxInd[selInd + 1]
#       temp <- x[selMaxInd1] - x[selMaxInd2]
#       localMax[selMaxInd1[temp <= 0]] <- 0
#       localMax[selMaxInd2[temp > 0]] <- 0
#     }
#     
#     localMax
#   }

# MSW.getRidge <-
#   function(localMax,
#            iInit = ncol(localMax),
#            step = -1,
#            iFinal = 1,
#            minWinSize = 3,
#            gapTh = 3,
#            skip = NULL)
#   {
#     ## modified from package MassSpecWavelet
#     
#     scales <- as.numeric(colnames(localMax))
#     if (is.null(scales))
#       scales <- 1:ncol(localMax)
#     
#     maxInd_curr <- which(localMax[, iInit] > 0)
#     nMz <- nrow(localMax)
#     
#     if (is.null(skip))	{
#       skip <- iInit + 1
#     }
#     
#     ## Identify all the peak pathes from the coarse level to detail levels (high column to low column)
#     ## Only consider the shortest path
#     if (ncol(localMax) > 1)
#       colInd <- seq(iInit + step, iFinal, step)
#     else
#       colInd <- 1
#     ridgeList <- as.list(maxInd_curr)
#     names(ridgeList) <- maxInd_curr
#     peakStatus <- as.list(rep(0, length(maxInd_curr)))
#     names(peakStatus) <- maxInd_curr
#     
#     ## orphanRidgeList keep the ridges disconnected at certain scale level
#     ## Changed by Pan Du 05/11/06
#     orphanRidgeList <- NULL
#     orphanRidgeName <- NULL
#     nLevel <- length(colInd)
#     
#     for (j in 1:nLevel) {
#       col.j <- colInd[j]
#       scale.j <- scales[col.j]
#       
#       if (colInd[j] == skip) {
#         oldname <- names(ridgeList)
#         ridgeList <-
#           lapply(ridgeList, function(x)
#             c(x, x[length(x)]))
#         ##peakStatus <- lapply(peakStatus, function(x) c(x, x[length(x)]))
#         names(ridgeList) <- oldname
#         ##names(peakStatus) <- oldname
#         next
#       }
#       
#       if (length(maxInd_curr) == 0) {
#         maxInd_curr <- which(localMax[, col.j] > 0)
#         next
#       }
#       
#       ## The slide window size is proportional to the CWT scale
#       ## winSize.j <- scale.j / 2 + 1
#       winSize.j <- floor(scale.j / 2)
#       if (winSize.j < minWinSize) {
#         winSize.j <- minWinSize
#       }
#       
#       selPeak.j <- NULL
#       remove.j <- NULL
#       for (k in 1:length(maxInd_curr)) {
#         ind.k <- maxInd_curr[k]
#         start.k <-
#           ifelse(ind.k - winSize.j < 1, 1, ind.k - winSize.j)
#         end.k <-
#           ifelse(ind.k + winSize.j > nMz, nMz, ind.k + winSize.j)
#         ind.curr <-
#           which(localMax[start.k:end.k, col.j] > 0) + start.k - 1
#         ##ind.curr <- which(localMax[, col.j] > 0)
#         if (length(ind.curr) == 0) {
#           status.k <- peakStatus[[as.character(ind.k)]]
#           ## bug  work-around
#           if (is.null(status.k))
#             status.k <- gapTh + 1
#           ##
#           if (status.k > gapTh & scale.j >= 2) {
#             temp <- ridgeList[[as.character(ind.k)]]
#             orphanRidgeList <-
#               c(orphanRidgeList, list(temp[1:(length(temp) - status.k)]))
#             orphanRidgeName <-
#               c(orphanRidgeName,
#                 paste(col.j + status.k + 1, ind.k, sep = '_'))
#             remove.j <- c(remove.j, as.character(ind.k))
#             next
#           } else {
#             ind.curr <- ind.k
#             peakStatus[[as.character(ind.k)]] <-
#               status.k + 1
#           }
#         } else {
#           peakStatus[[as.character(ind.k)]] <- 0
#           if (length(ind.curr) >= 2)
#             ind.curr <- ind.curr[which.min(abs(ind.curr - ind.k))]
#         }
#         ridgeList[[as.character(ind.k)]] <-
#           c(ridgeList[[as.character(ind.k)]], ind.curr)
#         selPeak.j <- c(selPeak.j, ind.curr)
#       }
#       ## Remove the disconnected lines from the currrent list
#       if (length(remove.j) > 0) {
#         removeInd <- which(names(ridgeList) %in% remove.j)
#         ridgeList <- ridgeList[-removeInd]
#         peakStatus <- peakStatus[-removeInd]
#       }
#       
#       ## Check for duplicated selected peaks and only keep the one with the longest path.
#       dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
#       if (length(dupPeak.j) > 0) {
#         removeInd <- NULL
#         for (dupPeak.jk in dupPeak.j) {
#           selInd <- which(selPeak.j == dupPeak.jk)
#           selLen <- sapply(ridgeList[selInd], length)
#           removeInd.jk <- which.max(selLen)
#           removeInd <- c(removeInd, selInd[-removeInd.jk])
#           orphanRidgeList <-
#             c(orphanRidgeList, ridgeList[removeInd.jk])
#           orphanRidgeName <-
#             c(orphanRidgeName, paste(col.j, selPeak.j[removeInd.jk], sep = '_'))
#         }
#         selPeak.j <- selPeak.j[-removeInd]
#         ridgeList <- ridgeList[-removeInd]
#         peakStatus <- peakStatus[-removeInd]
#       }
#       
#       ## Update the names of the ridgeList as the new selected peaks
#       ##if (scale.j >= 2) {
#       if (length(ridgeList) > 0)
#         names(ridgeList) <- selPeak.j
#       if (length(peakStatus) > 0)
#         names(peakStatus) <- selPeak.j
#       ##}
#       
#       ## If the level is larger than 3, expand the peak list by including other unselected peaks at that level
#       if (scale.j >= 2) {
#         maxInd_next <- which(localMax[, col.j] > 0)
#         unSelPeak.j <-
#           maxInd_next[!(maxInd_next %in% selPeak.j)]
#         newPeak.j <- as.list(unSelPeak.j)
#         names(newPeak.j) <- unSelPeak.j
#         ## Update ridgeList
#         ridgeList <- c(ridgeList, newPeak.j)
#         maxInd_curr <- c(selPeak.j, unSelPeak.j)
#         ## Update peakStatus
#         newPeakStatus <- as.list(rep(0, length(newPeak.j)))
#         names(newPeakStatus) <- newPeak.j
#         peakStatus <- c(peakStatus, newPeakStatus)
#       } else {
#         maxInd_curr <- selPeak.j
#       }
#     }
#     
#     ## Attach the peak level at the beginning of the ridge names
#     if (length(ridgeList) > 0)
#       names(ridgeList) <- paste(1, names(ridgeList), sep = '_')
#     if (length(orphanRidgeList) > 0)
#       names(orphanRidgeList) <- orphanRidgeName
#     ## Combine ridgeList and orphanRidgeList
#     ridgeList <- c(ridgeList, orphanRidgeList)
#     if (length(ridgeList) == 0)
#       return(NULL)
#     
#     ## Reverse the order as from the low level to high level.
#     ridgeList <- lapply(ridgeList, rev)
#     ## order the ridgeList in increasing order
#     ##ord <- order(selPeak.j)
#     ##ridgeList <- ridgeList[ord]
#     
#     ## Remove possible duplicated ridges
#     ridgeList <- ridgeList[!duplicated(names(ridgeList))]
#     
#     attr(ridgeList, 'class') <- 'ridgeList'
#     attr(ridgeList, 'scales') <- scales
#     return(ridgeList)
#   }

running <-
  function (X,
            Y = NULL,
            fun = mean,
            width = min(length(X), 20),
            allow.fewer = FALSE,
            pad = FALSE,
            align = c("right", "center",
                      "left"),
            simplify = TRUE,
            by,
            ...)
  {
    ## from package gtools
    align = match.arg(align)
    n <- length(X)
    if (align == "left") {
      from <- 1:n
      to <- pmin((1:n) + width - 1, n)
    }
    else if (align == "right") {
      from <- pmax((1:n) - width + 1, 1)
      to <- 1:n
    }
    else {
      from <- pmax((2 - width):n, 1)
      to <- pmin(1:(n + width - 1), n)
      if (!odd(width))
        stop("width must be odd for center alignment")
    }
    elements <- apply(cbind(from, to), 1, function(x)
      seq(x[1],
          x[2]))
    if (is.matrix(elements))
      elements <- as.data.frame(elements)
    names(elements) <- paste(from, to, sep = ":")
    if (!allow.fewer) {
      len <- sapply(elements, length)
      skip <- (len < width)
    }
    else {
      skip <- 0
    }
    run.elements <- elements[!skip]
    if (!invalid(by))
      run.elements <-
      run.elements[seq(from = 1,
                       to = length(run.elements),
                       by = by)]
    if (is.null(Y)) {
      funct <- function(which, what, fun, ...)
        fun(what[which],
            ...)
      if (simplify)
        Xvar <- sapply(run.elements, funct, what = X, fun = fun,
                       ...)
      else
        Xvar <- lapply(run.elements, funct, what = X, fun = fun,
                       ...)
    } else {
      funct <- function(which, XX, YY, fun, ...)
        fun(XX[which],
            YY[which], ...)
      if (simplify)
        Xvar <- sapply(
          run.elements,
          funct,
          XX = X,
          YY = Y,
          fun = fun,
          ...
        )
      else
        Xvar <- lapply(
          run.elements,
          funct,
          XX = X,
          YY = Y,
          fun = fun,
          ...
        )
    }
    if (allow.fewer || !pad)
      return(Xvar)
    if (simplify)
      if (is.matrix(Xvar)) {
        wholemat <- matrix(new(class(Xvar[1, 1]), NA),
                           ncol = length(to),
                           nrow = nrow(Xvar))
        colnames(wholemat) <- paste(from, to, sep = ":")
        wholemat[, -skip] <- Xvar
        Xvar <- wholemat
      }
    else {
      wholelist <- rep(new(class(Xvar[1]), NA), length(from))
      names(wholelist) <- names(elements)
      wholelist[names(Xvar)] <- Xvar
      Xvar <- wholelist
    }
    return(Xvar)
  }

invalid <- function (x)
{
  ## from package gtools
  if (missing(x) || is.null(x) || length(x) == 0)
    return(TRUE)
  if (is.list(x))
    return(all(sapply(x, invalid)))
  else if (is.vector(x))
    return(all(is.na(x)))
  else
    return(FALSE)
}

odd <- function (x)
  x != as.integer(x / 2) * 2


gauss <- function(x, h, mu, sigma) {
  h * exp(-(x - mu) ^ 2 / (2 * sigma ^ 2))
}

fitGauss <- function(td, d, pgauss = NA) {
  if (length(d) < 3)
    return(rep(NA, 3))
  if (!any(is.na(pgauss))) {
    mu <- pgauss$mu
    sigma <- pgauss$sigma
    h <- pgauss$h
  }
  fit <- try(nls(d ~ SSgauss(td, mu, sigma, h)), silent = TRUE)
  if (is(fit, "try-error"))
    fit <-
    try(nls(d ~ SSgauss(td, mu, sigma, h), algorithm = 'port'),
        silent = TRUE)
  if (is(fit, "try-error"))
    return(rep(NA, 3))
  
  as.data.frame(t(fit$m$getPars()))
}

## ' @param
## ' @param d numeric vector with intensities of centroids within the peak.
## ' @param otd
## ' @param omz
## ' @param od
## ' @param scantime
## ' @param scan.range
## ' @param peaks
## ' @noRd
# joinOverlappingPeaks <-
#   function(td,
#            d,
#            otd,
#            omz,
#            od,
#            scantime,
#            scan.range,
#            peaks,
#            maxGaussOverlap = 0.5,
#            mzCenterFun) {
#     ## Fix issue #284: avoid having identical peaks multiple times in this
#     ## matrix.
#     peaks <- unique(peaks)
#     gausspeaksidx <- which(!is.na(peaks[, "mu"]))
#     Ngp <- length(gausspeaksidx)
#     if (Ngp == 0)
#       return(peaks)
#     
#     newpeaks <- NULL
#     
#     gpeaks <- peaks[gausspeaksidx, , drop = FALSE]
#     if (nrow(peaks) - Ngp > 0)
#       notgausspeaks <- peaks[-gausspeaksidx, , drop = FALSE]
#     
#     if (Ngp > 1) {
#       comb <- which(upper.tri(matrix(0, Ngp, Ngp)), arr.ind = TRUE)
#       overlap <- logical(nrow(comb))
#       overlap <- rep(FALSE, dim(comb)[1])
#       for (k in seq_len(nrow(comb))) {
#         p1 <- comb[k, 1]
#         p2 <- comb[k, 2]
#         overlap[k] <- gaussCoverage(
#           xlim = scan.range,
#           h1 = gpeaks[p1, "h"],
#           mu1 = gpeaks[p1, "mu"],
#           s1 = gpeaks[p1, "sigma"],
#           h2 = gpeaks[p2, "h"],
#           mu2 = gpeaks[p2, "mu"],
#           s2 = gpeaks[p2, "sigma"]
#         ) >=
#           maxGaussOverlap
#       }
#     } else
#       overlap <- FALSE
#     
#     if (any(overlap) && (Ngp > 1)) {
#       jlist <- list()
#       if (length(which(overlap)) > 1) {
#         gm <- comb[overlap, ]
#         ## create list of connected components
#         cc <- list()
#         cc[[1]] <- gm[1,] ## copy first entry
#         for (j in 2:dim(gm)[1]) {
#           ## search for connections
#           ccl <- unlist(cc)
#           nl <- sapply(cc, function(x)
#             length(x))
#           ccidx <- rep(1:length(nl), nl)
#           idx <- match(gm[j,], ccl)
#           if (any(!is.na(idx))) {
#             ## connection found, add to list
#             pos <- ccidx[idx[which(!is.na(idx))[1]]]
#             cc[[pos]] <- c(cc[[pos]], gm[j,])
#           } else
#             ## create new list element
#             cc[[length(cc) + 1]] <- gm[j,]
#           
#         }
#         ccn <- list()
#         lcc <- length(cc)
#         ins <- rep(FALSE, lcc)
#         if (lcc > 1) {
#           jcomb <- which(upper.tri(matrix(0, lcc, lcc)), arr.ind = TRUE)
#           for (j in 1:dim(jcomb)[1]) {
#             j1 <- jcomb[j, 1]
#             j2 <- jcomb[j, 2]
#             if (any(cc[[j1]] %in% cc[[j2]]))
#               ccn[[length(ccn) + 1]] <-
#               unique(c(cc[[j1]], cc[[j2]]))
#             else {
#               if (!ins[j1]) {
#                 ccn[[length(ccn) + 1]] <- unique(cc[[j1]])
#                 ins[j1] <- TRUE
#               }
#               if (!ins[j2]) {
#                 ccn[[length(ccn) + 1]] <- unique(cc[[j2]])
#                 ins[j2] <- TRUE
#               }
#             }
#           }
#         } else
#           ccn <- cc
#         
#         
#         size <- sapply(ccn, function(x)
#           length(x))
#         s2idx <- which(size >= 2)
#         
#         if (length(s2idx) > 0) {
#           for (j in 1:length(s2idx)) {
#             pgroup <- unique(ccn[[s2idx[j]]])
#             jlist[[j]] <- pgroup
#           }
#         } else
#           stop('(length(s2idx) = 0) ?!?')
#       } else
#         jlist[[1]] <- comb[overlap, ]
#       
#       ## join all peaks belonging to one cc
#       for (j in seq_along(jlist)) {
#         jidx <- jlist[[j]]
#         newpeak <- gpeaks[jidx[1], , drop = FALSE]
#         newmin <- min(gpeaks[jidx, "lmin"])
#         newmax <- max(gpeaks[jidx, "lmax"])
#         newpeak[1, "scpos"] <- -1 ## not defined after join
#         newpeak[1, "scmin"] <- -1 ##    ..
#         newpeak[1, "scmax"] <- -1 ##    ..
#         newpeak[1, "scale"] <- -1 ##    ..
#         
#         newpeak[1, "maxo"] <- max(gpeaks[jidx, "maxo"])
#         newpeak[1, "sn"]   <- max(gpeaks[jidx, "sn"])
#         newpeak[1, "lmin"] <- newmin
#         newpeak[1, "lmax"] <- newmax
#         newpeak[1, "rtmin"] <- scantime[td[newmin]]
#         newpeak[1, "rtmax"] <- scantime[td[newmax]]
#         newpeak[1, "rt"] <- weighted.mean(gpeaks[jidx, "rt"],
#                                           w = gpeaks[jidx, "maxo"])
#         
#         ## Re-assign m/z values
#         p1 <- match(td[newmin], otd)[1]
#         p2 <- match(td[newmax], otd)
#         p2 <- p2[length(p2)]
#         if (is.na(p1))
#           p1 <- 1
#         if (is.na(p2))
#           p2 <- length(omz)
#         mz.value <- omz[p1:p2]
#         mz.int <- od[p1:p2]
#         
#         ## re-calculate m/z value for peak range
#         mzmean <- do.call(mzCenterFun, list(mz = mz.value,
#                                             intensity = mz.int))
#         mzrange <- range(mz.value)
#         newpeak[1, "mz"] <- mzmean
#         newpeak[1, c("mzmin", "mzmax")] <- mzrange
#         
#         ## re-fit gaussian
#         md <- max(d[newmin:newmax])
#         d1 <- d[newmin:newmax] / md
#         pgauss <- fitGauss(td[newmin:newmax],
#                            d[newmin:newmax],
#                            pgauss = list(
#                              mu = td[newmin] +
#                                (td[newmax] - td[newmin]) /
#                                2,
#                              sigma = td[newmax] - td[newmin],
#                              h = max(gpeaks[jidx, "h"])
#                            ))
#         if (!any(is.na(pgauss)) && all(pgauss > 0)) {
#           newpeak[1, "mu"]    <- pgauss$mu
#           newpeak[1, "sigma"] <- pgauss$sigma
#           newpeak[1, "h"]     <- pgauss$h
#           newpeak[1, "egauss"] <-
#             sqrt((1 / length(td[newmin:newmax])) *
#                    sum(((
#                      d1 - gauss(td[newmin:newmax],
#                                 pgauss$h / md,
#                                 pgauss$mu,
#                                 pgauss$sigma)
#                    ) ^ 2)))
#         } else {
#           ## re-fit after join failed
#           newpeak[1, "mu"]       <- NA
#           newpeak[1, "sigma"]    <- NA
#           newpeak[1, "h"]        <- NA
#           newpeak[1, "egauss"]   <- NA
#         }
#         
#         newpeaks <- rbind(newpeaks, newpeak)
#       }
#       ## add the remaining peaks
#       jp <- unique(unlist(jlist))
#       if (dim(peaks)[1] - length(jp) > 0)
#         newpeaks <- rbind(newpeaks, gpeaks[-jp, ])
#       
#     } else
#       newpeaks <- gpeaks
#     
#     grt.min <- newpeaks[, "rtmin"]
#     grt.max <- newpeaks[, "rtmax"]
#     
#     if (nrow(peaks) - Ngp > 0) {
#       ## notgausspeaks
#       for (k in 1:nrow(notgausspeaks)) {
#         ## here we can only check if they are completely overlapped
#         ## by other peaks
#         if (!any((notgausspeaks[k, "rtmin"] >= grt.min) &
#                  (notgausspeaks[k, "rtmax"] <= grt.max)))
#           newpeaks <- rbind(newpeaks, notgausspeaks[k,])
#       }
#     }
#     
#     rownames(newpeaks) <- NULL
#     newpeaks
#   }

# descendMinTol <- function(d, startpos, maxDescOutlier) {
#   l <- startpos[1]
#   r <- startpos[2]
#   outl <- 0
#   N <- length(d)
#   ## left
#   while ((l > 1) && (d[l] > 0) && outl <= maxDescOutlier) {
#     if (outl > 0)
#       vpos <- opos
#     else
#       vpos <- l
#     if (d[l - 1] > d[vpos])
#       outl <- outl + 1
#     else
#       outl <- 0
#     if (outl == 1)
#       opos <- l
#     l <- l - 1
#   }
#   if (outl > 0)
#     l <- l + outl
#   ## right
#   outl <- 0
#   
#   while ((r < N) && (d[r] > 0) && outl <= maxDescOutlier) {
#     if (outl > 0)
#       vpos <- opos
#     else
#       vpos <- r
#     if (d[r + 1] > d[vpos])
#       outl <- outl + 1
#     else
#       outl <- 0
#     if (outl == 1)
#       opos <- r
#     r <- r + 1
#   }
#   if (outl > 0)
#     r <- r - outl
#   c(l, r)
# }

cent <- function(x) {
  N <- length(x)
  if (N == 1)
    return(1)
  floor(N / 2)
}

gaussCoverage <- function(xlim, h1, mu1, s1, h2, mu2, s2) {
  overlap <- NA
  by = 0.05
  ## Calculate points of intersection
  a <- s2 ^ 2 - s1 ^ 2
  cc <-
    -(2 * s1 ^ 2 * s2 ^ 2 * (log(h1) - log(h2)) + (s1 ^ 2 * mu2 ^ 2) - (s2 ^
                                                                          2 * mu1 ^ 2))
  b <- ((2 * s1 ^ 2 * mu2) - (2 * s2 ^ 2 * mu1))
  D <- b ^ 2 - (a * cc)
  if (a == 0) {
    S1 <- -cc / b
    S2 <- NA
  } else if ((D < 0) || ((b ^ 2 - (4 * a * cc)) < 0)) {
    S1 <- S2 <- NA
  } else {
    S1 <- (-b + sqrt(b ^ 2 - (4 * a * cc))) / (2 * a)
    S2 <- (-b - sqrt(b ^ 2 - (4 * a * cc))) / (2 * a)
    if (S2 < S1)
    {
      tmp <- S1
      S1 <- S2
      S2 <- tmp
    }
  }
  if (!is.na(S1))
    if (S1 < xlim[1] || S1 > xlim[2])
      S1 <- NA
  if (!is.na(S2))
    if (S2 < xlim[1] || S2 > xlim[2])
      S2 <- NA
  
  x <- seq(xlim[1], xlim[2], by = by)
  vsmall <-
    min(sum(gauss(x, h1, mu1, s1)), sum(gauss(x, h2, mu2, s2)))
  
  if (!is.na(S1) && !is.na(S2)) {
    x0 <- seq(xlim[1], S1, by = by)
    xo <- seq(S1, S2, by = by)
    x1 <- seq(S2, xlim[2], by = by)
    if (gauss(x0[cent(x0)], h1, mu1, s1) < gauss(x0[cent(x0)], h2, mu2, s2)) {
      ov1 <- sum(gauss(x0, h1, mu1, s1))
    } else {
      ov1 <- sum(gauss(x0, h2, mu2, s2))
    }
    if (gauss(xo[cent(xo)], h1, mu1, s1) < gauss(xo[cent(xo)], h2, mu2, s2)) {
      ov <- sum(gauss(xo, h1, mu1, s1))
    } else {
      ov <- sum(gauss(xo, h2, mu2, s2))
    }
    if (gauss(x1[cent(x1)], h1, mu1, s1) < gauss(x1[cent(x1)], h2, mu2, s2)) {
      ov2 <- sum(gauss(x1, h1, mu1, s1))
    } else {
      ov2 <- sum(gauss(x1, h2, mu2, s2))
    }
    overlap <- ov1 + ov + ov2
  } else
    if (is.na(S1) &&
        is.na(S2)) {
      ## no overlap -> intergrate smaller function
      if (gauss(x[cent(x)], h1, mu1, s1) < gauss(x[cent(x)], h2, mu2, s2)) {
        overlap <- sum(gauss(x, h1, mu1, s1))
      } else {
        overlap <- sum(gauss(x, h2, mu2, s2))
      }
    } else
      if (!is.na(S1) || !is.na(S2)) {
        if (is.na(S1))
          S0 <- S2
        else
          S0 <- S1
        x0 <- seq(xlim[1], S0, by = by)
        x1 <- seq(S0, xlim[2], by = by)
        g01 <- gauss(x0[cent(x0)], h1, mu1, s1)
        g02 <- gauss(x0[cent(x0)], h2, mu2, s2)
        g11 <- gauss(x1[cent(x1)], h1, mu1, s1)
        g12 <- gauss(x1[cent(x1)], h2, mu2, s2)
        if (g01 < g02)
          ov1 <-
          sum(gauss(x0, h1, mu1, s1))
        else
          ov1 <- sum(gauss(x0, h2, mu2, s2))
        if (g11 < g12)
          ov2 <-
          sum(gauss(x1, h1, mu1, s1))
        else
          ov2 <- sum(gauss(x1, h2, mu2, s2))
        if ((g01 == g02) && (g01 == 0))
          ov1 <- 0
        if ((g11 == g12) && (g11 == 0))
          ov2 <- 0
        overlap <- ov1 + ov2
      }
  
  overlap / vsmall
}

mzCenter.wMean <- function(mz, intensity) {
  weighted.mean(mz, intensity)
}

mzCenter.mean <- function(mz, intensity) {
  mean(mz)
}

mzCenter.apex <- function(mz, intensity) {
  mz[which.max(intensity)]
}

mzCenter.wMeanApex3 <- function(mz, intensity) {
  iap <- which.max(intensity)
  st <- max(1, iap - 1)
  en <- min(iap + 1, length(mz))
  weighted.mean(mz[st:en], intensity[st:en])
}

mzCenter.meanApex3 <- function(mz, intensity) {
  iap <- which.max(intensity)
  st <- max(1, iap - 1)
  en <- min(iap + 1, length(mz))
  mean(mz[st:en])
}

# trimm <- function(x, trim = c(0.05, 0.95)) {
#   a <- sort(x[x > 0])
#   Na <- length(a)
#   quant <- round((Na * trim[1]) + 1):round(Na * trim[2])
#   a[quant]
# }

# estimateChromNoise <- function(x, trim = 0.05, minPts = 20) {
#   gz <- which(x > 0)
#   if (length(gz) < minPts)
#     return(mean(x))
#   
#   mean(x[gz], trim = trim)
# }

# getLocalNoiseEstimate <-
#   function(d,
#            td,
#            ftd,
#            noiserange,
#            Nscantime,
#            threshold,
#            num) {
#     if (length(d) < Nscantime) {
#       ## noiserange[2] is full d-range
#       drange <- which(td %in% ftd)
#       n1 <- d[-drange] ## region outside the detected ROI (wide)
#       n1.cp <-
#         continuousPtsAboveThresholdIdx(n1, threshold = threshold, num = num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
#       n1 <- n1[!n1.cp]
#       if (length(n1) > 1)  {
#         baseline1 <- mean(n1)
#         sdnoise1 <- sd(n1)
#       } else
#         baseline1 <- sdnoise1 <- 1
#       
#       ## noiserange[1]
#       d1 <- drange[1]
#       d2 <- drange[length(drange)]
#       nrange2 <-
#         c(max(1, d1 - noiserange[1]):d1,
#           d2:min(length(d), d2 + noiserange[1]))
#       n2 <- d[nrange2] ## region outside the detected ROI (narrow)
#       n2.cp <-
#         continuousPtsAboveThresholdIdx(n2, threshold = threshold, num = num) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
#       n2 <- n2[!n2.cp]
#       if (length(n2) > 1)  {
#         baseline2 <- mean(n2)
#         sdnoise2 <- sd(n2)
#       } else
#         baseline2 <- sdnoise2 <- 1
#       
#     } else {
#       trimmed <- trimm(d, c(0.05, 0.95))
#       baseline1 <- baseline2 <- mean(trimmed)
#       sdnoise1 <- sdnoise2 <- sd(trimmed)
#     }
#     
#     c(min(baseline1, baseline2), min(sdnoise1, sdnoise2))
#   }


###remove dot
# narrow_rt_boundaries <- function(lm, d, thresh = 1) {
#   lm_seq <- lm[1]:lm[2]
#   above_thresh <- d[lm_seq] >= thresh
#   if (any(above_thresh)) {
#     ## Expand by one on each side to be consistent with old code.
#     above_thresh <- above_thresh | c(above_thresh[-1], FALSE) |
#       c(FALSE, above_thresh[-length(above_thresh)])
#     lm <- range(lm_seq[above_thresh], na.rm = TRUE)
#   }
#   lm
# }


# rectUnique <-
#   function(m,
#            order = seq(length = nrow(m)),
#            xdiff = 0,
#            ydiff = 0) {
#     nr <- nrow(m)
#     nc <- ncol(m)
#     
#     if (is.null(nr) || nr == 0) {
#       ## empty matrix in first place
#       return (m)
#     }
#     
#     if (!is.double(m))
#       m <- as.double(m)
#     .C(
#       "RectUnique",
#       m,
#       as.integer(order - 1),
#       nr,
#       nc,
#       as.double(xdiff),
#       as.double(ydiff),
#       logical(nrow(m)),
#       PACKAGE = "xcms"
#     )[[7]]
#   }


.fdata <- function(x) {
  x@featureData@data
}





# valueCount2ScanIndex <- function(valCount) {
#   ## Convert into 0 based.
#   valCount <- cumsum(valCount)
#   return(as.integer(c(0, valCount[-length(valCount)])))
# }


# continuousPtsAboveThreshold <-
#   function(y, threshold, num, istart = 1) {
#     if (!is.double(y))
#       y <- as.double(y)
#     if (.C(
#       "continuousPtsAboveThreshold",
#       y,
#       as.integer(istart - 1),
#       length(y),
#       threshold = as.double(threshold),
#       num = as.integer(num),
#       n = integer(1),
#       PACKAGE = "xcms"
#     )$n > 0)
#       TRUE
#     else
#       FALSE
#   }



# continuousPtsAboveThresholdIdx <-
#   function(y, threshold, num, istart = 1) {
#     if (!is.double(y))
#       y <- as.double(y)
#     as.logical(
#       .C(
#         "continuousPtsAboveThresholdIdx",
#         y,
#         as.integer(istart - 1),
#         length(y),
#         threshold = as.double(threshold),
#         num = as.integer(num),
#         n = integer(length(y)),
#         PACKAGE = "xcms"
#       )$n
#     )
#   }


###remove the dot
peaks_to_result <-
  function(res, object, startDate, param, msLevel) {
    xph <- XProcessHistory(
      param = param,
      date. = startDate,
      type. = .PROCSTEP.PEAK.DETECTION,
      fileIndex = 1:length(fileNames(object)),
      msLevel = msLevel
    )
    object <- as(object, "XCMSnExp")
    phist <- object@.processHistory
    ## if (hasAdjustedRtime(object) | hasFeatures(object))
    ##     object@msFeatureData <- new("MsFeatureData")
    pks <- do.call(rbind, res$peaks)
    if (length(pks) > 0) {
      chromPeaks(object) <- pks
      chromPeakData(object)$ms_level <- as.integer(msLevel)
      chromPeakData(object)$is_filled <- FALSE
    }
    object@.processHistory <- c(phist, list(xph))
    validObject(object)
    object
  }



###remove dot
processResultList <- function(x, getProcHist = TRUE, fnames) {
  nSamps <- length(x)
  pks <- vector("list", nSamps)
  phList <- vector("list", nSamps)
  for (i in 1:nSamps) {
    n_pks <- nrow(x[[i]]$peaks)
    if (is.null(n_pks))
      n_pks <- 0
    if (n_pks == 0) {
      pks[[i]] <- NULL
      warning("No peaks found in sample number ", i, ".")
    } else {
      pks[[i]] <- cbind(x[[i]]$peaks, sample = rep.int(i, n_pks))
    }
    if (getProcHist)
      phList[[i]] <- ProcessHistory(
        info. = paste0(
          "Chromatographic peak detection in '",
          basename(fnames[i]),
          "': ",
          n_pks,
          " peaks identified."
        ),
        date. = x[[i]]$date,
        type. = .PROCSTEP.PEAK.DETECTION,
        fileIndex. = i
      )
  }
  list(peaks = pks, procHist = phList)
}


XProcessHistory <-
  function(param = NULL,
           msLevel = NA_integer_,
           ...) {
    obj <- ProcessHistory(...)
    obj <- as(obj, "XProcessHistory")
    obj@param <- param
    obj@msLevel <- as.integer(msLevel)
    OK <- validObject(obj)
    if (is.character(OK))
      stop(OK)
    return(obj)
  }



descendMin <- function(y, istart = which.max(y)) {
  if (!is.double(y))
    y <- as.double(y)
  unlist(
    .C(
      "DescendMin",
      y,
      length(y),
      as.integer(istart - 1),
      ilower = integer(1),
      iupper = integer(1),
      PACKAGE = "xcms"
    )[4:5]
  ) + 1
}


.rawMat <-
  function(mz,
           int,
           scantime,
           valsPerSpect,
           mzrange = numeric(),
           rtrange = numeric(),
           scanrange = numeric(),
           log = FALSE) {
    if (length(rtrange) >= 2) {
      rtrange <- range(rtrange)
      ## Fix for issue #267. rtrange outside scanrange causes scanrange
      ## being c(Inf, -Inf)
      scns <-
        which((scantime >= rtrange[1]) & (scantime <= rtrange[2]))
      if (!length(scns))
        return(matrix(
          nrow = 0,
          ncol = 3,
          dimnames = list(character(), c("time", "mz", "intensity"))
        ))
      scanrange <- range(scns)
    }
    if (length(scanrange) < 2)
      scanrange <- c(1, length(valsPerSpect))
    else
      scanrange <- range(scanrange)
    if (!all(is.finite(scanrange)))
      stop("'scanrange' does not contain finite values")
    if (!all(is.finite(mzrange)))
      stop("'mzrange' does not contain finite values")
    if (!all(is.finite(rtrange)))
      stop("'rtrange' does not contain finite values")
    if (scanrange[1] == 1)
      startidx <- 1
    else
      startidx <- sum(valsPerSpect[1:(scanrange[1] - 1)]) + 1
    endidx <- sum(valsPerSpect[1:scanrange[2]])
    scans <- rep(scanrange[1]:scanrange[2],
                 valsPerSpect[scanrange[1]:scanrange[2]])
    masses <- mz[startidx:endidx]
    massidx <- 1:length(masses)
    if (length(mzrange) >= 2) {
      mzrange <- range(mzrange)
      massidx <-
        massidx[(masses >= mzrange[1] & (masses <= mzrange[2]))]
    }
    int <- int[startidx:endidx][massidx]
    if (log && (length(int) > 0))
      int <- log(int + max(1 - min(int), 0))
    cbind(time = scantime[scans[massidx]],
          mz = masses[massidx],
          intensity = int)
  }




ProcessHistory <-
  function(type., date., info., error., fileIndex.) {
    if (missing(type.))
      type. <- .PROCSTEP.UNKNOWN
    if (missing(info.))
      info. <- character()
    if (missing(error.))
      error. <- NULL
    if (missing(date.))
      date. <- date()
    if (missing(fileIndex.))
      fileIndex. <- integer()
    return(
      new(
        "ProcessHistory",
        type = type.,
        info = info.,
        date = date.,
        error = error.,
        fileIndex = as.integer(fileIndex.)
      )
    )
  }


hasChromPeaks <- function(object) {
  as.logical(nrow(object@chromPeaks))
}


.PROCSTEP.UNKNOWN <- "Unknown"
.PROCSTEP.PEAK.DETECTION <- "Peak detection"
.PROCSTEP.PEAK.REFINEMENT <- "Peak refinement"
.PROCSTEP.PEAK.GROUPING <- "Peak grouping"
.PROCSTEP.RTIME.CORRECTION <- "Retention time correction"
.PROCSTEP.PEAK.FILLING <- "Missing peak filling"
.PROCSTEP.CALIBRATION <- "Calibration"
.PROCSTEP.FEATURE.GROUPING <- "Feature grouping"
.PROCSTEPS <- c(
  .PROCSTEP.UNKNOWN,
  .PROCSTEP.PEAK.DETECTION,
  .PROCSTEP.PEAK.REFINEMENT,
  .PROCSTEP.PEAK.GROUPING,
  .PROCSTEP.RTIME.CORRECTION,
  .PROCSTEP.PEAK.FILLING,
  .PROCSTEP.CALIBRATION,
  .PROCSTEP.FEATURE.GROUPING
)