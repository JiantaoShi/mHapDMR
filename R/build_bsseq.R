#' Build a DSS BSseq object from per-sample rGR objects
#'
#' Takes a named list of \code{GRanges} objects produced by
#' \code{\link{mscore_region_stats}} (one per sample) and assembles a
#' \code{BSseq} object using the kappa/Y' substitution, so that DSS operates
#' on M-score rather than raw read fractions.
#'
#' Because every sample shares the \emph{same} regions in the \emph{same} order
#' (they all derive from one \code{rGR}), the coverage (\eqn{\kappa}) and
#' methylation (\eqn{Y'}) matrices are assembled directly by column-binding the
#' per-sample vectors and the \code{BSseq} object is constructed in one shot.
#' This deliberately avoids \code{DSS::makeBSseqData()}, which merges samples
#' one-by-one with \code{merge(..., all = TRUE)} — if any two regions collide on
#' the \code{(chr, pos)} key (region midpoints can, via overlap or integer
#' rounding) that merge forms a Cartesian product that compounds across samples
#' (\eqn{2^{n_{samples}}} rows) and exhausts memory.  We also guarantee the
#' \code{(chr, pos)} key is unique so the downstream coordinate join is correct.
#'
#' @param rGR_list  A \strong{named} list of \code{GRanges} objects, one per
#'   sample, each being the output of \code{mscore_region_stats()}.
#' @param min_reads Minimum read depth required in \emph{every} sample for a
#'   region to be retained (default 1).
#' @return A list with two elements:
#'   \describe{
#'     \item{bsseq}{A \code{BSseq} object ready for \code{\link{run_dml_test}}.}
#'     \item{region_index}{A \code{data.frame} mapping every retained region to
#'       its DSS key (\code{chr}, \code{pos}) plus the original coordinates
#'       (\code{start}, \code{end}) and name (\code{region_id}).  Use this to
#'       join \code{dml_result} back to the input \code{rGR}.}
#'   }
#' @importFrom bsseq BSseq
#' @importFrom GenomicRanges mcols start end seqnames
#' @export
build_bsseq_mscore <- function(rGR_list, min_reads = 1L) {
    if (is.null(names(rGR_list))) stop("rGR_list must be a named list; names are used as sample identifiers.")

    rGR_ref <- rGR_list[[1L]]
    n_regions <- length(rGR_ref)

    # ── coverage filter (all samples must pass) ──────────────────────────────
    keep <- rep(TRUE, n_regions)
    for (rGR in rGR_list) {
        mc <- GenomicRanges::mcols(rGR)
        keep <- keep & !is.na(mc$reads) & (mc$reads >= min_reads) & !is.na(mc$kappa)
    }

    n_keep <- sum(keep)
    if (n_keep == 0L) stop("No regions passed the min_reads filter across all samples.")
    message(sprintf("Regions retained: %d / %d", n_keep, n_regions))

    # ── region coordinates (shared by all samples) ───────────────────────────
    chr_kept <- as.character(GenomicRanges::seqnames(rGR_ref)[keep])
    start_kept <- GenomicRanges::start(rGR_ref)[keep]
    end_kept <- GenomicRanges::end(rGR_ref)[keep]
    pos_kept <- as.integer((start_kept + end_kept) / 2L)

    # ── enforce a unique (chr, pos) DSS key ──────────────────────────────────
    # The midpoint key must be unique: BSseq/DSS index loci by (chr, pos) and the
    # result is joined back to coordinates on that key. Resolve any collisions by
    # nudging duplicate positions by +1 until unique (rare; preserves order).
    key <- paste0(chr_kept, ":", pos_kept)
    if (anyDuplicated(key)) {
        warning(sprintf(
            "%d region(s) share a (chr, midpoint) key; nudging positions to keep the DSS key unique.",
            sum(duplicated(key))))
        repeat {
            dup <- duplicated(key)
            if (!any(dup)) break
            pos_kept[dup] <- pos_kept[dup] + 1L
            key <- paste0(chr_kept, ":", pos_kept)
        }
    }

    # region_id: use existing names(rGR_ref) if set, otherwise auto-generate
    existing_names <- names(rGR_ref)
    if (!is.null(existing_names) && !any(is.na(existing_names)) && !any(existing_names == "")) {
        region_id <- existing_names[keep]
    } else {
        region_id <- paste0(chr_kept, ":", start_kept, "-", end_kept)
    }

    region_index <- data.frame(
        region_id = region_id,
        chr = chr_kept,
        start = start_kept,
        end = end_kept,
        pos = pos_kept, # DSS key — matches dml_result$pos
        stringsAsFactors = FALSE
    )

    # ── assemble Cov (kappa) and M (Y') matrices directly ────────────────────
    # All samples share the same kept regions in the same order, so a column-bind
    # is exact and O(n) — no per-sample merge needed.
    Cov <- vapply(rGR_list, function(rGR) as.double(GenomicRanges::mcols(rGR)$kappa[keep]),   numeric(n_keep))
    M   <- vapply(rGR_list, function(rGR) as.double(GenomicRanges::mcols(rGR)$Y_prime[keep]), numeric(n_keep))
    Cov[is.na(Cov)] <- 0
    M[is.na(M)] <- 0
    colnames(Cov) <- colnames(M) <- names(rGR_list)

    bsseq <- bsseq::BSseq(chr = chr_kept, pos = pos_kept, M = M, Cov = Cov,
                          sampleNames = names(rGR_list))

    list(bsseq = bsseq, region_index = region_index)
}
