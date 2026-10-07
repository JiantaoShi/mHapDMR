# ── internal helpers ──────────────────────────────────────────────────────────

#' Fetch CpG positions within a window from a tabix-indexed CpG file
#'
#' @param cpg_file     Path to bgzipped, tabix-indexed CpG file.
#'                      Expected columns: chr, pos (1-based).
#' @param chrom        Chromosome name.
#' @param region_start 1-based inclusive start.
#' @param region_end   1-based inclusive end.
#' @param margin       Extra bp added on each side (default 150).
#' @return Integer vector of CpG positions (1-based), sorted.
#' @keywords internal
.fetch_cpg_positions <- function(cpg_file, chrom, region_start, region_end, margin = 150L) {
    gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(max(1L, region_start - margin), region_end + margin))
    raw <- tryCatch(Rsamtools::scanTabix(cpg_file, param = gr)[[1L]], error = function(e) character(0L))
    if (length(raw) == 0L) return(integer(0L))
    pos <- vapply(strsplit(raw, "\t", fixed = TRUE), function(x) as.integer(x[2L]), integer(1L))
    sort(pos)
}

#' Parse and trim one mHap line to a region of interest
#'
#' @param fields        Character vector from splitting one mHap line on TAB.
#'                       Expected: chr, start, end, hapStr, count.
#' @param cpg_positions Integer vector of CpG positions in the region (1-based).
#' @param region_start  1-based inclusive start.
#' @param region_end    1-based inclusive end.
#' @return Named list(chrom, read_start, read_end, hap_met, count, ncpg,
#'         is_methylated) or NULL if the record is unusable.
#' @keywords internal
.parse_mhap_record <- function(fields, cpg_positions, region_start, region_end) {
    if (length(fields) < 5L) return(NULL)

    chrom   <- fields[1L]
    h_start <- as.integer(fields[2L])
    h_end   <- as.integer(fields[3L])
    hap_str <- fields[4L]
    count   <- as.integer(fields[5L])
    hap_chars <- strsplit(hap_str, "", fixed = TRUE)[[1L]]
        
    # CpG sites covered by this read
    frag_cpgs <- cpg_positions[cpg_positions >= h_start & cpg_positions <= h_end]
    if (length(frag_cpgs) != length(hap_chars)) return(NULL)

    # Restrict to sites within the analysis region
    keep <- frag_cpgs >= region_start & frag_cpgs <= region_end
    if (!any(keep)) return(NULL)

    hap_kept  <- hap_chars[keep]
    kept_cpgs <- frag_cpgs[keep]
    ncpg      <- length(kept_cpgs)
    mcpg      <- sum(hap_kept == "1")

    list(
        chrom         = chrom,
        read_start    = kept_cpgs[1L],
        read_end      = kept_cpgs[length(kept_cpgs)],
        hap_met       = hap_kept,
        count         = count,
        ncpg          = ncpg,
        mcpg          = mcpg,
        is_methylated = mcpg > 0,
        is_disordered = (ncpg > mcpg) & (mcpg > 0)
    )
}


#' Fetch and parse mHap records for a single interval (no GRanges assembly)
#'
#' Lightweight core shared by \code{\link{read_mhap_gr}} and the M-score
#' statistics path.  Returns a plain list of parsed records so callers that
#' only need summary statistics can avoid constructing a \code{GRanges} per
#' region (which is the dominant memory cost when scanning many regions).
#'
#' @inheritParams read_mhap_gr
#' @return A list of parsed records (see \code{.parse_mhap_record}); empty list
#'   if no usable reads overlap the interval.
#' @keywords internal
.read_mhap_records <- function(mhap_file, cpg_file, gr, margin = 150L) {
    chrom        <- as.character(GenomicRanges::seqnames(gr))
    region_start <- GenomicRanges::start(gr)
    region_end   <- GenomicRanges::end(gr)

    query_gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(region_start, region_end))
    raw <- tryCatch(Rsamtools::scanTabix(mhap_file, param = query_gr)[[1L]], error = function(e) character(0L))
    if (length(raw) == 0L) return(list())
    fields <- strsplit(raw, "\t", fixed = TRUE)

    # CpGs over the region and the full span of every read, so that no read
    # fails the CpG-count check by extending beyond the lookup window
    h_start <- suppressWarnings(as.integer(vapply(fields, function(x) x[2L], "")))
    h_end   <- suppressWarnings(as.integer(vapply(fields, function(x) x[3L], "")))
    cpg_pos <- .fetch_cpg_positions(cpg_file, chrom, min(region_start, h_start, na.rm = TRUE),
                                    max(region_end, h_end, na.rm = TRUE), margin)
    if (length(cpg_pos) == 0L) return(list())

    records <- lapply(fields, function(f) .parse_mhap_record(f, cpg_pos, region_start, region_end))
    Filter(Negate(is.null), records)
}


# ── exported functions ────────────────────────────────────────────────────────

#' Read mHap records for a single genomic interval into a GRanges (mHapGR)
#'
#' Fetches all mHap reads overlapping a genomic interval, trims each read to
#' the interval boundaries, and returns a \code{GRanges} object where each
#' range represents one unique haplotype.  Metadata columns store the
#' information needed to compute M-score.
#'
#' @param mhap_file Path to mHap file (.mhap.gz, bgzipped + tabix-indexed).
#' @param cpg_file  Path to CpG position file (.gz, bgzipped + tabix-indexed).
#' @param gr        A single-interval \code{GRanges} (length 1) defining the
#'                  query region.
#' @param margin    Extra bp added to the CpG position lookup (default 150).
#'   The lookup always covers the region and the full span of every read
#'   overlapping it, so the margin does not change the results.
#' @return A \code{GRanges} object (\strong{mHapGR}) with one range per unique
#'         haplotype and the following metadata columns:
#'         \describe{
#'             \item{count}{Integer. Haplotype multiplicity (mHap column 5).}
#'             \item{ncpg}{Integer. Number of CpG sites in the region on this read.}
#'             \item{mcpg}{Integer. Number of methylated CpG sites on this read.}
#'             \item{hap_str}{Character. Trimmed methylation string ("0"/"1" per CpG).}
#'             \item{is_methylated}{Logical. TRUE if the read carries >=1 methylated site.}
#'             \item{is_disordered}{Logical. TRUE if the read has both 0 and 1 states.}
#'         }
#'         Returns an empty \code{GRanges} with those columns if no reads overlap.
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools scanTabix
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' \dontrun{
#' gr <- GRanges("chr1", IRanges(10000L, 10300L))
#' mHapGR <- read_mhap_gr("sample.mhap.gz", "hg19_CpG.gz", gr)
#' mHapGR
#' }
read_mhap_gr <- function(mhap_file, cpg_file, gr, margin = 150L) {
    stopifnot(inherits(gr, "GRanges"), length(gr) == 1L)

    # Empty mHapGR template (returned on early exit)
    empty_gr <- GenomicRanges::GRanges()
    GenomicRanges::mcols(empty_gr) <- S4Vectors::DataFrame(
        count         = integer(0L),
        ncpg          = integer(0L),
        mcpg          = integer(0L),
        hap_str       = character(0L),
        is_methylated = logical(0L),
        is_disordered = logical(0L)
    )

    records <- .read_mhap_records(mhap_file, cpg_file, gr, margin)
    if (length(records) == 0L) return(empty_gr)

    # Assemble mHapGR
    mHapGR <- GenomicRanges::GRanges(
        seqnames = vapply(records, `[[`, character(1L), "chrom"),
        ranges   = IRanges::IRanges(
            start = vapply(records, `[[`, integer(1L), "read_start"),
            end   = vapply(records, `[[`, integer(1L), "read_end")
        )
    )
    GenomicRanges::mcols(mHapGR) <- S4Vectors::DataFrame(
        count         = vapply(records, `[[`, integer(1L),  "count"),
        ncpg          = vapply(records, `[[`, integer(1L),  "ncpg"),
        mcpg          = vapply(records, `[[`, integer(1L),  "mcpg"),
        hap_str       = vapply(records, function(r) paste(r$hap_met, collapse = ""), character(1L)),
        is_methylated = vapply(records, `[[`, logical(1L),  "is_methylated"),
        is_disordered = vapply(records, `[[`, logical(1L),  "is_disordered")
    )
    mHapGR
}
