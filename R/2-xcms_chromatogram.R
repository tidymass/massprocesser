# chromatogram <-
#   function(object,
#            rt,
#            mz,
#            aggregationFun = "sum",
#            missing = NA_real_,
#            msLevel = 1L,
#            BPPARAM = bpparam(),
#            adjustedRtime = hasAdjustedRtime(object),
#            filled = FALSE,
#            include = c("apex_within", "any", "none")) {
#     include <- match.arg(include)
#     if (adjustedRtime)
#       adj_rt <- rtime(object, adjusted = TRUE)
#     object_od <- as(object, "OnDiskMSnExp")
#     fcs <-
#       c(
#         "fileIdx",
#         "spIdx",
#         "seqNum",
#         "acquisitionNum",
#         "msLevel",
#         "polarity",
#         "retentionTime",
#         "precursorScanNum"
#       )
#     fcs <- intersect(fcs, colnames(.fdata(object)))
#     object_od <- MSnbase::selectFeatureData(object_od, fcol = fcs)
#     if (adjustedRtime)
#       object_od@featureData$retentionTime <- adj_rt
#     res <- MSnbase::chromatogram(
#       object_od,
#       rt = rt,
#       mz = mz,
#       aggregationFun = aggregationFun,
#       missing = missing,
#       msLevel = msLevel,
#       BPPARAM = BPPARAM
#     )
#     if (!hasChromPeaks(object) | include == "none")
#       return(res)
#     ## Process peaks
#     lvls <- 1:length(fileNames(object))
#     if (missing(rt))
#       rt <- c(-Inf, Inf)
#     if (missing(mz))
#       mz <- c(-Inf, Inf)
#     if (is.matrix(rt) | is.matrix(mz)) {
#       ## Ensure rt and mz are aligned.
#       if (!is.matrix(rt))
#         rt <- matrix(rt, ncol = 2)
#       if (!is.matrix(mz))
#         mz <- matrix(mz, ncol = 2)
#       if (nrow(rt) == 1)
#         rt <- matrix(rep(rt, nrow(mz)), ncol = 2, byrow = TRUE)
#       if (nrow(mz) == 1)
#         mz <- matrix(rep(mz, nrow(rt)), ncol = 2, byrow = TRUE)
#       pk_list <- vector("list", nrow(mz))
#       pkd_list <- vector("list", nrow(mz))
#       for (i in 1:nrow(mz)) {
#         pks <- chromPeaks(object,
#                           rt = rt[i,],
#                           mz = mz[i,],
#                           type = include)
#         pkd <- chromPeakData(object)[rownames(pks), , drop = FALSE]
#         if (!filled) {
#           pks <- pks[!pkd$is_filled, , drop = FALSE]
#           pkd <- extractROWS(pkd, which(!pkd$is_filled))
#         }
#         smpls <- factor(pks[, "sample"], levels = lvls)
#         pk_list[[i]] <- split.data.frame(pks, smpls)
#         pkd_list[[i]] <- split.data.frame(pkd, smpls)
#       }
#       pks <- do.call(rbind, pk_list)
#       pks <- pks[seq_along(pks)]
#       pkd <- do.call(rbind, pkd_list)
#       pkd <- pkd[seq_along(pkd)]
#     } else {
#       pks <- chromPeaks(object,
#                         rt = rt,
#                         mz = mz,
#                         type = include)
#       pkd <- chromPeakData(object)[rownames(pks), , drop = FALSE]
#       if (!filled) {
#         pks <- pks[!pkd$is_filled, , drop = FALSE]
#         pkd <- extractROWS(pkd, which(!pkd$is_filled))
#       }
#       smpls <- factor(pks[, "sample"], levels = lvls)
#       pks <- split.data.frame(pks, smpls)
#       pkd <- split.data.frame(pkd, smpls)
#       mz <- matrix(mz, ncol = 2)
#       rt <- matrix(rt, ncol = 2)
#     }
#     res <- as(res, "XChromatograms")
#     res@.Data <- matrix(
#       mapply(
#         unlist(res),
#         pks,
#         pkd,
#         FUN = function(chr, pk, pd) {
#           chr@chromPeaks <- pk
#           chr@chromPeakData <- pd
#           chr
#         }
#       ),
#       nrow = nrow(res),
#       dimnames = dimnames(res)
#     )
#     res@.processHistory <- object@.processHistory
#     if (hasFeatures(object)) {
#       pks_sub <- chromPeaks(res)
#       ## Loop through each EIC "row" to ensure all features in
#       ## that EIC are retained.
#       fts <- lapply(seq_len(nrow(res)), function(r) {
#         fdev <- featureDefinitions(object, mz = mz[r,],
#                                    rt = rt[r,])
#         if (nrow(fdev)) {
#           fdev$row <- r
#           .subset_features_on_chrom_peaks(fdev, chromPeaks(object), pks_sub)
#         } else
#           DataFrame()
#       })
#       res@featureDefinitions <- do.call(rbind, fts)
#     }
#     validObject(res)
#     res
#     
#   }
