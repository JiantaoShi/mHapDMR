#' Accumulate M-score statistics from per-read vectors
#'
#' Core arithmetic shared by the GRanges- and record-based stat functions.
#' Computes the weighted sums that define the M-score and the DSS substitution
#' parameters \eqn{\kappa} and \eqn{Y'} from aligned per-read vectors.
#'
#' @param N   Integer vector of per-read CpG counts (\code{ncpg}).
#' @param Z   Logical/integer vector, methylated state per read
#'   (\code{is_methylated}).
#' @param cnt Integer vector of haplotype multiplicities (\code{count}).
#' @return A named numeric vector with elements:
#'   \describe{
#'     \item{reads}{Total read count (sum of \code{count}).}
#'     \item{Nsum}{\eqn{\sum_t N_t c_t} — weighted sum of per-read CpG counts.}
#'     \item{N2sum}{\eqn{\sum_t N_t^2 c_t} — weighted sum of squared CpG counts.}
#'     \item{Sjd}{\eqn{\sum_t N_t Z_t c_t} — M-score numerator.}
#'     \item{mscore}{\eqn{S_{jd}/T_{jd}} — the M-score (NA if Nsum = 0).}
#'     \item{kappa}{\eqn{T_{jd}^2 / \sum N_t^2 c_t} (NA if N2sum = 0).}
#'     \item{Y_prime}{\eqn{S_{jd} \cdot T_{jd} / \sum N_t^2 c_t} (NA if N2sum = 0).}
#'   }
#' @keywords internal
.mscore_accumulate <- function(N, Z, cnt) {
    out <- c(reads = 0, Nsum = 0, N2sum = 0, Sjd = 0, mscore = NA_real_, kappa = NA_real_, Y_prime = NA_real_)
    if (length(N) == 0L) return(out)

    # Accumulate in double precision: weighted sums for deep regions can exceed
    # the 2^31 integer limit and would otherwise silently overflow to NA.
    N   <- as.double(N)
    Z   <- as.double(Z)
    cnt <- as.double(cnt)

    reads <- sum(cnt)
    Nsum  <- sum(N * cnt)
    N2sum <- sum(N * N * cnt)
    Sjd   <- sum(N * Z * cnt)

    out["reads"] <- reads
    out["Nsum"]  <- Nsum
    out["N2sum"] <- N2sum
    out["Sjd"]   <- Sjd

    if (N2sum > 0) {
        out["kappa"]   <- Nsum^2 / N2sum
        out["Y_prime"] <- Sjd * Nsum / N2sum
    }
    if (Nsum > 0) {
        out["mscore"] <- Sjd / Nsum
    }
    out
}

#' @keywords internal
.mscore_from_mhapgr <- function(mHapGR) {
    if (length(mHapGR) == 0L) return(.mscore_accumulate(integer(0L), integer(0L), integer(0L)))
    mc <- GenomicRanges::mcols(mHapGR)
    .mscore_accumulate(mc$ncpg, mc$is_methylated, mc$count)
}

#' Compute M-score statistics directly from raw mHap lines (streaming, no objects)
#'
#' Fully vectorised, allocation-light computation of the per-region M-score
#' statistics straight from the raw tabix lines.  No per-read list, no
#' methylation strings, and no \code{GRanges} are materialised — only a handful
#' of integer/character vectors of length \emph{number of haplotypes in the
#' region}, which keeps the per-region (and therefore per-worker) memory
#' high-water-mark minimal at scale.
#'
#' Reproduces the exact semantics of \code{\link{.parse_mhap_record}}: a read is
#' used only if the CpG positions spanning \code{[start, end]} match the length
#' of its methylation string, and statistics are restricted to CpGs that fall
#' inside the analysis region.
#'
#' @param raw          Character vector of raw mHap lines for the region
#'   (output of \code{Rsamtools::scanTabix(...)[[1]]}). Columns:
#'   chr, start, end, hapStr, count, and an optional strand.
#' @param cpg_pos      Integer vector of CpG positions (1-based, \strong{sorted
#'   ascending}) covering the region plus margin.
#' @param region_start 1-based inclusive analysis-region start.
#' @param region_end   1-based inclusive analysis-region end.
#' @return Same named numeric vector as \code{.mscore_from_mhapgr}.
#' @keywords internal
.mscore_streaming <- function(raw, cpg_pos, region_start, region_end) {
    empty <- c(reads = 0, Nsum = 0, N2sum = 0, Sjd = 0,
               mscore = NA_real_, kappa = NA_real_, Y_prime = NA_real_)
    if (length(raw) == 0L || length(cpg_pos) == 0L) return(empty)

    # Parse columns 2-5 straight into typed vectors (col 1 skipped, 6+ flushed)
    cols <- tryCatch(
        scan(text = raw, sep = "\t", quiet = TRUE, flush = TRUE,
             what = list(NULL, integer(), integer(), character(), integer())),
        error = function(e) NULL
    )
    if (is.null(cols) || length(cols[[2L]]) == 0L) {
        # Robust fallback for any non-conforming lines
        fields <- strsplit(raw, "\t", fixed = TRUE)
        keepf  <- lengths(fields) >= 5L
        if (!any(keepf)) return(empty)
        fields  <- fields[keepf]
        h_start <- as.integer(vapply(fields, `[`, "", 2L))
        h_end   <- as.integer(vapply(fields, `[`, "", 3L))
        hap_str <- vapply(fields, `[`, "", 4L)
        count   <- as.integer(vapply(fields, `[`, "", 5L))
    } else {
        h_start <- cols[[2L]]; h_end <- cols[[3L]]
        hap_str <- cols[[4L]]; count <- cols[[5L]]
    }

    # CpG index span covered by each read's haplotype string (cpg_pos sorted)
    lo    <- findInterval(h_start - 1L, cpg_pos) + 1L  # first CpG idx >= h_start
    hi    <- findInterval(h_end,        cpg_pos)        # last  CpG idx <= h_end
    nspan <- hi - lo + 1L                               # CpGs the hapStr must cover

    # CpG index span restricted to the analysis region
    rlo  <- findInterval(pmax(h_start, region_start) - 1L, cpg_pos) + 1L
    rhi  <- findInterval(pmin(h_end,   region_end),        cpg_pos)
    ncpg <- rhi - rlo + 1L

    valid <- !is.na(count) & (nspan == nchar(hap_str)) & (ncpg > 0L)
    if (!any(valid)) return(empty)

    # Methylated-in-region <=> any "1" among the region-restricted characters
    cs <- (rlo - lo + 1L)[valid]
    ce <- (rhi - lo + 1L)[valid]
    Z  <- as.integer(grepl("1", substr(hap_str[valid], cs, ce), fixed = TRUE))

    .mscore_accumulate(ncpg[valid], Z, count[valid])
}


#' Compute M-score statistics for all regions in a rGR across one sample
#'
#' Iterates over every interval in \code{rGR} and computes the per-region
#' M-score statistics, storing them back as metadata columns of \code{rGR}.
#'
#' By default (\code{keep_mhapgr = FALSE}) this uses a streaming path that reads
#' each region's raw mHap lines and reduces them straight to summary statistics
#' via \code{\link{.mscore_streaming}} — \strong{no per-read objects, methylation
#' strings, or \code{GRanges} are ever materialised}, so memory stays flat even
#' for tens of thousands of regions across many samples.  Tabix handles are
#' opened once and reused across regions.  Set \code{keep_mhapgr = TRUE} only
#' when you also need the raw haplotypes (builds one \code{GRanges} per region).
#'
#' @param mhap_file Path to mHap file (.mhap.gz, bgzipped + tabix-indexed).
#' @param cpg_file  Path to CpG position file (.gz, bgzipped + tabix-indexed).
#' @param rGR       A \code{GRanges} of analysis regions.
#' @param margin    Extra bp for CpG position lookup (default 150).
#' @param verbose   Print progress every 100 regions (default TRUE).
#' @return The input \code{rGR} with additional metadata columns:
#'   \describe{
#'     \item{reads}{Total read depth.}
#'     \item{Nsum, N2sum, Sjd}{Raw weighted sums (see \code{.mscore_from_mhapgr}).}
#'     \item{mscore}{M-score (\eqn{S_{jd}/T_{jd}}).}
#'     \item{kappa}{Effective sample size for DSS (\eqn{T^2 / \sum N_t^2}).}
#'     \item{Y_prime}{Effective methylation count for DSS.}
#'   }
#'   When \code{keep_mhapgr = TRUE}, a list attribute \code{"mHapGR_list"} is
#'   attached, containing one \code{GRanges} (mHapGR) per region.  This is off
#'   by default because retaining one GRanges per region per sample dominates
#'   memory use at scale (e.g. tens of thousands of regions); only the numeric
#'   statistics above are needed by the downstream DMR pipeline.
#' @param keep_mhapgr Logical. If \code{TRUE}, retain the raw per-region
#'   \code{GRanges} objects as the \code{"mHapGR_list"} attribute. Default
#'   \code{FALSE} for memory efficiency.
#' @importFrom Rsamtools TabixFile scanTabix
#' @export
#' @examples
#' \dontrun{
#' rGR <- GRanges("chr1", IRanges(c(1000L, 5000L), c(1200L, 5300L)))
#' rGR <- mscore_region_stats("sample.mhap.gz", "hg19_CpG.gz", rGR)
#' rGR$mscore
#' }
mscore_region_stats <- function(mhap_file, cpg_file, rGR, margin = 150L, verbose = TRUE, keep_mhapgr = FALSE) {
    stopifnot(inherits(rGR, "GRanges"))
    n <- length(rGR)

    # Pull coordinate vectors once to avoid repeated S4 subsetting in the loop
    seqn  <- as.character(GenomicRanges::seqnames(rGR))
    starts <- GenomicRanges::start(rGR)
    ends   <- GenomicRanges::end(rGR)

    stat_cols <- c("reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime")
    stat_mat  <- matrix(NA_real_, nrow = n, ncol = length(stat_cols), dimnames = list(NULL, stat_cols))
    empty_stat <- c(reads = 0, Nsum = 0, N2sum = 0, Sjd = 0, mscore = NA_real_, kappa = NA_real_, Y_prime = NA_real_)
    mHapGR_list <- if (keep_mhapgr) vector("list", n) else NULL

    # Reuse open tabix handles across all regions (streaming path only)
    if (!keep_mhapgr) {
        cpg_tf  <- Rsamtools::TabixFile(cpg_file)
        mhap_tf <- Rsamtools::TabixFile(mhap_file)
        open(cpg_tf); open(mhap_tf)
        on.exit({ close(cpg_tf); close(mhap_tf) }, add = TRUE)
    }

    for (i in seq_len(n)) {
        if (verbose && i %% 100L == 0L) {
            message(sprintf("[%d/%d] %s:%d-%d", i, n, seqn[i], starts[i], ends[i]))
        }

        if (keep_mhapgr) {
            gr_i <- GenomicRanges::GRanges(seqn[i], IRanges::IRanges(starts[i], ends[i]))
            mHapGR_i <- read_mhap_gr(mhap_file, cpg_file, gr_i, margin)
            mHapGR_list[[i]] <- mHapGR_i
            stat_mat[i, ] <- .mscore_from_mhapgr(mHapGR_i)
        } else {
            # Streaming path: raw lines -> summary stats, nothing retained.
            cpg_pos <- .fetch_cpg_positions(cpg_tf, seqn[i], starts[i], ends[i], margin)
            if (length(cpg_pos) == 0L) { stat_mat[i, ] <- empty_stat; next }
            region_gr <- GenomicRanges::GRanges(seqn[i], IRanges::IRanges(starts[i], ends[i]))
            raw <- tryCatch(Rsamtools::scanTabix(mhap_tf, param = region_gr)[[1L]],
                            error = function(e) character(0L))
            stat_mat[i, ] <- .mscore_streaming(raw, cpg_pos, starts[i], ends[i])
        }
    }

    # Store stats as metadata columns of rGR
    for (col in stat_cols) {
        GenomicRanges::mcols(rGR)[[col]] <- stat_mat[, col]
    }

    if (keep_mhapgr) attr(rGR, "mHapGR_list") <- mHapGR_list
    rGR
}
