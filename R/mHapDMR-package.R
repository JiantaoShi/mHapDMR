#' mHapDMR: Differential M-score analysis from mHap files
#'
#' @description
#' \pkg{mHapDMR} provides an end-to-end pipeline for differential M-score
#' analysis using mHap-format sequencing data.
#'
#' **M-score definition**
#'
#' For region \eqn{j}, sample \eqn{d}, let \eqn{N_{tjd}} be the number of CpG
#' sites on read \eqn{t} and \eqn{Z_{tjd} \in \{0,1\}} indicate whether the
#' read carries at least one methylated site. The M-score is:
#'
#' \deqn{M_{jd} = \frac{S_{jd}}{T_{jd}}, \quad S_{jd} = \sum_t N_{tjd} Z_{tjd}, \quad T_{jd} = \sum_t N_{tjd}}
#'
#' **DSS substitution**
#'
#' To fit M-score into the DSS Beta-Binomial framework, standard DSS inputs
#' \eqn{(m_{jd}, Y_{jd})} are replaced by \eqn{(\kappa_{jd}, Y'_{jd})}:
#'
#' \deqn{\kappa_{jd} = \frac{T_{jd}^2}{\sum_t N_{tjd}^2}, \qquad Y'_{jd} = S_{jd} \cdot \frac{T_{jd}}{\sum_t N_{tjd}^2}}
#'
#' so that \eqn{Y'_{jd} / \kappa_{jd} = M_{jd}} exactly.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{mhap_dmr}}}{Full pipeline: mHap files → DML results.}
#'   \item{\code{\link{read_mhap_gr}}}{Read mHap records into a \code{GRanges} (mHapGR).}
#'   \item{\code{\link{mscore_region_stats}}}{Compute M-score stats for all regions (rGR).}
#'   \item{\code{\link{build_bsseq_mscore}}}{Build BSseq with kappa/Y' substitution.}
#'   \item{\code{\link{run_dml_test}}}{DSS DML test wrapper.}
#' }
#'
#' @docType package
#' @name mHapDMR-package
"_PACKAGE"
