# # mz = args[["mz"]]
# # int = args[["int"]]
# # scantime = args[["scantime"]]
# # valsPerSpect = args[["valsPerSpect"]]
# # ppm = args[["ppm"]]
# # peakwidth = args[["peakwidth"]]
# # snthresh = args[["snthresh"]]
# # prefilter = args[["prefilter"]]
# # mzCenterFun = args[["mzCenterFun"]]
# # integrate = args[["integrate"]]
# # mzdiff  = args[["mzdiff"]]
# # fitgauss  = args[["fitgauss"]]
# # noise = args[["noise"]]
# # verboseColumns = args[["verboseColumns"]]
# # roiList = args[["roiList"]]
# # firstBaselineCheck  = args[["firstBaselineCheck"]]
# # roiScales  = args[["roiScales"]]
# # sleep = 0
# # extendLengthMSW  = args[["extendLengthMSW"]]
# # file_name = args[["file_name"]]
# 
# 
# do_findChromPeaks_centWave_sxt <-
#   function(mz,
#            int,
#            scantime,
#            valsPerSpect,
#            ppm = 25,
#            peakwidth = c(20, 50),
#            snthresh = 10,
#            prefilter = c(3, 100),
#            mzCenterFun = "wMean",
#            integrate = 1,
#            mzdiff = -0.001,
#            fitgauss = FALSE,
#            noise = 0,
#            verboseColumns = FALSE,
#            roiList = list(),
#            firstBaselineCheck = TRUE,
#            roiScales = NULL,
#            sleep = 0,
#            extendLengthMSW = FALSE,
#            file_name = "file") {
#     if (getOption("originalCentWave", default = TRUE)) {
#       ## message("DEBUG: using original centWave.")
#       centWave_orig_sxt(
#         mz = mz,
#         int = int,
#         scantime = scantime,
#         valsPerSpect = valsPerSpect,
#         ppm = ppm,
#         peakwidth = peakwidth,
#         snthresh = snthresh,
#         prefilter = prefilter,
#         mzCenterFun = mzCenterFun,
#         integrate = integrate,
#         mzdiff = mzdiff,
#         fitgauss = fitgauss,
#         noise = noise,
#         verboseColumns = verboseColumns,
#         roiList = roiList,
#         firstBaselineCheck = firstBaselineCheck,
#         roiScales = roiScales,
#         sleep = sleep,
#         extendLengthMSW = extendLengthMSW,
#         file_name = file_name
#       )
#     } else {
#       ## message("DEBUG: using modified centWave.")
#       centWave_orig_sxt(
#         mz = mz,
#         int = int,
#         scantime = scantime,
#         valsPerSpect = valsPerSpect,
#         ppm = ppm,
#         peakwidth = peakwidth,
#         snthresh = snthresh,
#         prefilter = prefilter,
#         mzCenterFun = mzCenterFun,
#         integrate = integrate,
#         mzdiff = mzdiff,
#         fitgauss = fitgauss,
#         noise = noise,
#         verboseColumns = verboseColumns,
#         roiList = roiList,
#         firstBaselineCheck = firstBaselineCheck,
#         roiScales = roiScales,
#         sleep = sleep,
#         file_name = file_name
#       )
#     }
#   }