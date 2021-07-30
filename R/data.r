#' Simulated differential regulated data
#'
#' A dataset containing simulated data for 300 samples with 1000 genes
#' and 20 miRNAs. The first 10 are true positive.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{predict.target}{The predicted targets of each miRNAs}
#'   \item{group.label}{Group labels for each individuals}
#'   \item{mrna.mat}{Gene expression matrix}
#'   \item{mirna.mat}{miRNA expression matrix}
#'   ...
#' }
"sim.data"