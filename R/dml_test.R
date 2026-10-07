#' Run a DSS multi-factor differential M-score test
#'
#' @param bsseq     A \code{BSseq} object from \code{\link{build_bsseq_mscore}}.
#' @param design    A \code{data.frame} with one row per sample.
#' @param formula   One-sided formula (default \code{~ group}).
#' @param coef      Coefficient name to test (default \code{"groupTumor"}).
#' @param smoothing Passed to DSS (default FALSE).
#' @return A \code{data.frame} of DML results sorted by ascending p-value.
#' @importFrom DSS DMLfit.multiFactor DMLtest.multiFactor
#' @export
run_dml_test <- function(bsseq, design, formula = ~group, coef = "groupTumor") {
    fit <- DSS::DMLfit.multiFactor(bsseq, design = design, formula = formula, smoothing = FALSE)
    res <- DSS::DMLtest.multiFactor(fit, coef = coef)
    res[order(res$pvals), ]
}


#' Annotate DML results with original region coordinates
#'
#' Joins the \code{data.frame} returned by \code{\link{run_dml_test}} back to
#' the input regions using the \code{chr + pos} key, restoring \code{start},
#' \code{end}, and \code{region_id} columns.  Optionally converts the result
#' to a \code{GRanges}.
#'
#' @param dml_result  A \code{data.frame} from \code{\link{run_dml_test}}.
#' @param region_index A \code{data.frame} from
#'   \code{\link{build_bsseq_mscore}} (the \code{$region_index} element).
#' @param as_granges  If \code{TRUE}, return a \code{GRanges} object with
#'   DML statistics as metadata columns (default \code{FALSE}).
#' @return A \code{data.frame} (or \code{GRanges}) with columns:
#'   \code{region_id}, \code{chr}, \code{start}, \code{end}, \code{pos},
#'   \code{stat}, \code{pvals}, \code{fdrs}.
#' @importFrom GenomicRanges GRanges mcols<-
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' \dontrun{
#' result   <- mhap_dmr(...)
#' annotated <- annotate_dml_result(result$dml_result, result$region_index)
#' head(annotated)
#'
#' # As GRanges
#' dml_gr <- annotate_dml_result(result$dml_result, result$region_index, as_granges = TRUE)
#' dml_gr[dml_gr$fdrs < 0.05]
#' }
annotate_dml_result <- function(dml_result, region_index, as_granges = FALSE) {
    # join on chr + pos
    merged <- merge(dml_result, region_index, by = c("chr", "pos"), all.x = TRUE, sort = FALSE)
    # restore original sort order (by pvals)
    merged <- merged[order(merged$pvals), ]

    # tidy column order
    front <- c("region_id", "chr", "start", "end", "pos")
    rest <- setdiff(names(merged), front)
    merged <- merged[, c(front, rest), drop = FALSE]
    rownames(merged) <- NULL

    if (!as_granges) return(merged)

    # convert to GRanges
    gr <- GenomicRanges::GRanges(seqnames = merged$chr, ranges = IRanges::IRanges(start = merged$start, end = merged$end))
    stat_cols <- setdiff(names(merged), c("chr", "start", "end"))
    GenomicRanges::mcols(gr) <- S4Vectors::DataFrame(merged[, stat_cols, drop = FALSE])
    gr
}


#' Full pipeline: mHap files to annotated differential M-score results
#'
#' @param mhap_files Named character vector of mHap file paths (.mhap.gz).
#' @param cpg_file   Path to the CpG position file (bgzipped, tabix-indexed).
#' @param rGR        A \code{GRanges} of analysis regions.  If the object has
#'   \code{names()}, they are used as \code{region_id}; otherwise ids are
#'   auto-generated as \code{chr:start-end}.
#' @param design     A \code{data.frame} with one row per sample.
#' @param formula    Model formula (default \code{~ group}).
#' @param coef       Coefficient to test (default \code{"groupTumor"}).
#' @param min_reads  Minimum per-region read depth in all samples (default 1).
#' @param margin     Extra bp added to the CpG position lookup (default 150);
#'   does not change the results (see \code{\link{mscore_region_stats}}).
#' @param smoothing  Passed to DSS (default FALSE).
#' @param BPPARAM    \code{BiocParallelParam} for sample-level parallelism.
#' @param keep_mhapgr Retain raw per-region \code{GRanges} in each sample's
#'   result (default \code{FALSE}). Leave \code{FALSE} for large analyses;
#'   retaining them costs one \code{GRanges} per region per sample.
#' @param verbose    Print per-region progress (default TRUE).
#' @return A named list:
#'   \describe{
#'     \item{dml_result}{Annotated \code{data.frame} with \code{region_id},
#'       \code{start}, \code{end} restored, sorted by p-value.}
#'     \item{bsseq}{\code{BSseq} object used for testing.}
#'     \item{region_index}{\code{data.frame} mapping \code{chr+pos} ↔ original
#'       region coordinates and ids.}
#'     \item{rGR_list}{Named list of per-sample \code{GRanges} with M-score
#'       metadata columns.}
#'   }
#' @importFrom BiocParallel bplapply SerialParam
#' @export
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#'
#' # Optionally name your regions
#' rGR <- GRanges("chr1", IRanges(c(1000L, 5000L), c(1200L, 5300L)))
#' names(rGR) <- c("promoter_GENE1", "exon1_GENE2")
#'
#' # Continuous covariate (e.g. age)
#' design <- data.frame(age = c(25, 45, 55, 70), row.names = names(mhap_files))
#' result <- mhap_dmr(..., design = design, formula = ~ age, coef = "age")
#'
#' # Restore region coordinates
#' head(result$dml_result)   # already annotated
#'
#' # Significant regions as GRanges
#' dml_gr <- annotate_dml_result(result$dml_result, result$region_index, as_granges = TRUE)
#' dml_gr[dml_gr$fdrs < 0.05]
#' }
mhap_dmr <- function(mhap_files, cpg_file, rGR, design, formula = ~group, coef = "groupTumor", min_reads = 1L, margin = 150L, smoothing = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 4L), keep_mhapgr = FALSE, verbose = TRUE) {
    if (is.null(names(mhap_files))) stop("mhap_files must be a named character vector.")
    stopifnot(inherits(rGR, "GRanges"))

    # Step 1: per-sample M-score extraction
    message("=== Step 1: computing M-score statistics per sample ===")
    rGR_list <- BiocParallel::bplapply(names(mhap_files), function(sname) mscore_region_stats(mhap_files[[sname]], cpg_file, rGR, margin = margin, verbose = verbose, keep_mhapgr = keep_mhapgr), BPPARAM = BPPARAM)
    names(rGR_list) <- names(mhap_files)

    # Step 2: build BSseq + region index
    message("=== Step 2: building BSseq object ===")
    built <- build_bsseq_mscore(rGR_list, min_reads = min_reads)
    bsseq <- built$bsseq
    region_index <- built$region_index

    # Step 3: differential M-score test
    message("=== Step 3: differential M-score test (DSS) ===")
    dml_raw <- run_dml_test(bsseq, design, formula = formula, coef = coef)

    # Step 4: annotate results with original region coordinates
    dml_result <- annotate_dml_result(dml_raw, region_index)
    m_summary <- summary_m_score(rGR_list, design)
    list(dml_result = dml_result, bsseq = bsseq, region_index = region_index, rGR_list = rGR_list, m_summary = m_summary)
}
