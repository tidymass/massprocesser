#' @title massprocesser_logo
#' @description massprocesser_logo
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @importFrom xcms CentWaveParam findChromPeaks adjustRtime ObiwarpParam chromatogram PeakDensityParam groupChromPeaks
#' @importFrom xcms featureChromatograms groupnames featureDefinitions featureValues fillChromPeaks
#' @importFrom Biobase featureData
#' @importFrom crayon yellow red green bold bgRed
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom stringr str_split str_replace_all str_trim str_detect str_extract
#' @importFrom dplyr filter select pull everything distinct one_of left_join mutate bind_cols arrange
#' @importFrom tibble as_tibble enframe tibble rownames_to_column
#' @importFrom clisymbols symbol
#' @importFrom cli rule col_cyan tree
#' @importFrom utils packageVersion object.size write.csv tail
#' @importFrom purrr map map2
#' @importFrom plyr dlply .
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Biobase featureData
#' @importFrom MSnbase readMSData
#' @importFrom readr read_csv cols
#' @importFrom readxl read_excel
#' @importFrom masstools get_os mz_rt_match
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom BiocGenerics basename
#' @importFrom stats coefficients lm loess median predict
#' @importFrom stats rgamma rt sd cor p.adjust prcomp t.test wilcox.test
#' @importFrom methods is
#' @import mzR
#' @import tidyr
#' @import tidyverse
#' @importClassesFrom massdataset tidymass_parameter mass_dataset
#' @export
#' @return massprocesser_logo
#' @examples
#' massprocesser_logo()

massprocesser_logo <- function() {
  message(crayon::green("Thank you for using massprocesser!\n"))
  message(crayon::green("Version", massprocesser_version, "(", update_date, ')\n'))
  message(crayon::green("More information: searching 'tidymass massprocesser'.\n"))
  cat(crayon::green(
    c(
      "                 _    __ _              ___  ",
      "                | |  / _| |            |__ \\ ",
      "  _ __ ___   ___| |_| |_| | _____      __ ) |",
      " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
      " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ",
      " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
      "                                             ",
      "                                             "
    )
    
  ), sep = "\n")
}

# massprocesser_version = "0.99.3"
massprocesser_version <-
  as.character(utils::packageVersion(pkg = "massprocesser"))
update_date = as.character(Sys.time())


#' @title get_massprocesser_version
#' @description get_massprocesser_version
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @export
#' @return massprocesser_version
#' @examples
#' get_massprocesser_version()

get_massprocesser_version <- function() {
  return(massprocesser_version)
}

#' .onAttach <- function(libname, pkgname){
#'   packageStartupMessage(crayon::green(
#'     "massprocesser,
#' More information can be found at https://tidymass.github.io/massprocesser/
#' Authors: Xiaotao Shen (shenxt@stanford.edu)
#' Maintainer: Xiaotao Shen.
#' Version 0.1.2 (20201230)
#' Small bug fix."
#'
#'   )
#'   )
#' }
#'
#' # library(cowsay)
#' #https://onlineasciitools.com/convert-text-to-ascii-art
#' # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
#' # art <- readLines("ascii_art.txt")
#' # dput(art)
#' # metflow_logo <-
#' #   c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
#' #     "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
#' #     " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
#' #     "                                             ", "                                             "
#' #   )
#' # cat(metflow_logo, sep = "\n")
#' #
#' #
#' # library(asciify)
#' # bayes_img <- ascii_data("bayes.png")      # path to the bayes image
#' # bayes_map <- ascii_map(file = bayes_img)  # construct ASCII map
#' # bayes_map
#' # ascii_plot(bayes_map, charsize = 2)
#'
#'
#'
#' # cities <- c("S\u00e3o Paulo", "Reykjav\u00edk")
#' # print(cities)
#' # ASCIIfy(cities, 1)
#' # ASCIIfy(cities, 2)
#' #
#' # athens <- "\u0391\u03b8\u03ae\u03bd\u03b1"
#' # print(athens)
#' # ASCIIfy(athens)
#'
#'
#' #' # art <- readLines("metflo2_logo.txt")
#' #' # dput(art)
#' #'
#' #' #https://www.text-image.com/convert/ascii.html
#' #'
