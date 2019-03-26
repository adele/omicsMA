#' Intestinal Metaplasia Dataset
#'
#' This is part of the dataset used in the study by \insertCite{ribeiro2019omicsMA}{omicsMA}.
#'
#' It contains data from 10 two-color microarrays, being 5 from tissues representing
#' type II intestinal metaplasia and 5 from tissues representing the normal condition,
#' obtained from the Tumor Bank at A.C. Camargo Cancer Center / Antonio Prudente Foundation.
#'
#' Each sample is hybridized against a pool of normal tissues using the same orientation of
#' dye labeling. Gene expression levels were measured on Agilent Whole Human Genome
#' Microarrays 4x44K G4112F (design ID 014850).
#'
#' Each slide contais 45015 spots (41093 unique probes) and each spot contains
#'  about 60 foreground pixels.
#'
#' The scanned images of the microarray slides were processed by Agilent Feature Extraction
#' software, version 9.5, where statistics (mean, standard deviation and covariance) of the
#' foreground and local background pixels were computed for each spot, in both test and
#' reference channels.
#'
#' The experiment was conducted with financial support of the
#' Foundation for Research Support of the State of São Paulo (FAPESP), grant 06/03227-2,
#' and in accordance with the recommendations of the international
#' guidelines for investigations involving human beings with written informed consent from all subjects.
#' The protocol was approved by the Ethics Institutional Committee of the A.C. Camargo Cancer Center
#' (process number 1023/07).
#'
#'
#' @format A list with the following elements:
#' \describe{
#' \item{R.mean}{A 45015 x 10 matrix with the means of the red (test) foreground pixels in each spot of each array.}
#' \item{G.mean}{A 45015 x 10 matrix with the means of the green (reference) foreground pixels in each spot of each array.}
#' \item{R.var}{A 45015 x 10 matrix with the variances across red foreground pixels in each spot of each array.}
#' \item{G.var}{A 45015 x 10 matrix with the variances across green foreground pixels in each spot of each array.}
#' \item{RG.cov}{A 45015 x 10 matrix with the covariances between green and red foreground pixels in each spot of each array.}
#' \item{R.bckg}{A 45015 x 10 matrix with the means of the red background pixels in each spot of each array.}
#' \item{G.bckg}{A 45015 x 10 matrix with the means of the green background pixels in each spot of each array.}
#' \item{R.size}{A 45015 x 10 matrix with the number of red foreground pixels in each spot of each array.}
#' \item{G.size}{A 45015 x 10 matrix with the number of green foreground pixels in each spot of each array.}
#' \item{array.ids}{A vector with the ids for the 10 arrays.}
#' \item{gene.ids}{A vector with the ids for the 45015 spots.}
#' \item{gene.names}{A vector with the names of the 45015 genes.}
#' }
#'
#' @references{
#'    \insertAllCited
#' }
#'
"metaplasia"
