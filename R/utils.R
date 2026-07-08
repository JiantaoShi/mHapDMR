#' Check that a bgzipped file has a tabix index
#'
#' @param filepath Path to the bgzipped file (.gz).
#' @return Invisible \code{TRUE}; stops with an informative error if the
#'   companion \code{.tbi} index is absent.
#' @export
#' @examples
#' \dontrun{
#' check_tabix_index("hg19_CpG.gz")
#' check_tabix_index("sample.mhap.gz")
#' }
check_tabix_index <- function(filepath) {
    if (!file.exists(filepath)) stop("File not found: ", filepath)
    tbi <- paste0(filepath, ".tbi")
    if (!file.exists(tbi)) stop("Tabix index not found: ", tbi, "\nGenerate it with: tabix -p bed ", filepath)
    invisible(TRUE)
}

#' Summarize methylation score (mscore) by group
#'
#' Builds a matrix of group-wise mean mscore values from a list of GRanges
#' objects, where each list element corresponds to one sample and shares the
#' same set of intervals.
#'
#' @param rGR_list A named list of GRanges objects. Names must match the row
#'   names of \code{design}. Each element must contain an \code{mscore}
#'   metadata column, and all elements must share identical intervals in the
#'   same order.
#' @param design A data frame with a \code{group} column and row names matching
#'   the names of \code{rGR_list}.
#'
#' @return A numeric matrix with one row per interval (named
#'   \code{seqnames:start-end}) and one column per group, where each value is
#'   the mean \code{mscore} across the samples in that group.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- summary_m_score(rGR_list, design)
#' head(mat)
#' }
summary_m_score <- function(rGR_list, design){
    samples <- rownames(design)
    rGR_list <- rGR_list[samples]

    mscore_mat <- sapply(rGR_list, function(gr) gr$mscore)

    gr1 <- rGR_list[[1]]
    rownames(mscore_mat) <- paste0(seqnames(gr1), ":", start(gr1), "-", end(gr1))

    groups <- unique(as.character(design$group))
    result <- sapply(groups, function(g){
        cols <- samples[as.character(design$group) == g]
        rowMeans(mscore_mat[, cols, drop = FALSE], na.rm = TRUE)
    })

    return(result)
}
