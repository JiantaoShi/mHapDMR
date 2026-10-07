library(testthat)
library(GenomicRanges)

# ── helper: build a mock mHapGR ───────────────────────────────────────────────
make_mock_mhapgr <- function(N_vec, Z_vec, cnt_vec = rep(1L, length(N_vec)),
                              chrom = "chr1") {
  stopifnot(length(N_vec) == length(Z_vec),
            length(N_vec) == length(cnt_vec))
  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(start = seq_along(N_vec) * 100L,
                                width = N_vec * 2L)
  )
  GenomicRanges::mcols(gr) <- S4Vectors::DataFrame(
    count         = as.integer(cnt_vec),
    ncpg          = as.integer(N_vec),
    hap_str       = strrep("1", N_vec),       # dummy string
    is_methylated = as.logical(Z_vec)
  )
  gr
}


# ── kappa / Y' arithmetic ─────────────────────────────────────────────────────

test_that("kappa equals total reads when all reads have the same CpG count", {
  # N_t = k for all reads  →  kappa = (m*k)^2 / (m*k^2) = m
  N   <- rep(4L, 5L)
  cnt <- rep(1L, 5L)
  mHapGR <- make_mock_mhapgr(N, rep(1L, 5L), cnt)
  s <- mHapDMR:::.mscore_from_mhapgr(mHapGR)
  expect_equal(s[["kappa"]], 5, tolerance = 1e-10)
})

test_that("Y_prime / kappa equals M-score (ratio invariance)", {
  N   <- c(3L, 5L, 4L, 4L)
  Z   <- c(1L, 0L, 1L, 0L)
  cnt <- c(2L, 1L, 3L, 2L)

  mHapGR <- make_mock_mhapgr(N, Z, cnt)
  s <- mHapDMR:::.mscore_from_mhapgr(mHapGR)

  expect_equal(s[["Y_prime"]] / s[["kappa"]], s[["mscore"]], tolerance = 1e-12)
})

test_that("kappa satisfies Cauchy-Schwarz bound: 1 <= kappa <= sum(reads)", {
  set.seed(42L)
  for (i in seq_len(30L)) {
    N   <- sample(2L:8L, 5L, replace = TRUE)
    cnt <- sample(1L:5L, 5L, replace = TRUE)
    mHapGR <- make_mock_mhapgr(N, sample(0L:1L, 5L, replace = TRUE), cnt)
    s <- mHapDMR:::.mscore_from_mhapgr(mHapGR)
    expect_gte(s[["kappa"]], 1 - 1e-9)
    expect_lte(s[["kappa"]], sum(cnt) + 1e-9)
  }
})

test_that("M-score is in [0, 1]", {
  mHapGR <- make_mock_mhapgr(c(3L,4L,5L), c(1L,0L,1L), c(2L,2L,2L))
  s <- mHapDMR:::.mscore_from_mhapgr(mHapGR)
  expect_gte(s[["mscore"]], 0)
  expect_lte(s[["mscore"]], 1)
})

test_that("empty mHapGR returns all-zero / NA stats", {
  empty <- GenomicRanges::GRanges()
  GenomicRanges::mcols(empty) <- S4Vectors::DataFrame(
    count = integer(0L), ncpg = integer(0L),
    hap_str = character(0L), is_methylated = logical(0L)
  )
  s <- mHapDMR:::.mscore_from_mhapgr(empty)
  expect_equal(s[["reads"]], 0)
  expect_true(is.na(s[["mscore"]]))
  expect_true(is.na(s[["kappa"]]))
})


# ── streaming reducer (.mscore_streaming) ─────────────────────────────────────

test_that(".mscore_streaming trims reads to the region and matches hand calc", {
  # 5 CpGs at 100,102,104,106,108; analysis region = [102, 106]
  cpg_pos <- c(100L, 102L, 104L, 106L, 108L)
  raw <- c(
    "chr1\t100\t108\t10110\t2\t+",  # region chars "011" -> ncpg=3, meth=TRUE,  cnt=2
    "chr1\t102\t104\t00\t1\t+",     # region chars "00"  -> ncpg=2, meth=FALSE, cnt=1
    "chr1\t104\t108\t111\t3\t+"     # region chars "11"  -> ncpg=2, meth=TRUE,  cnt=3
  )
  s <- mHapDMR:::.mscore_streaming(raw, cpg_pos, 102L, 106L)

  expect_equal(s[["reads"]], 6)
  expect_equal(s[["Nsum"]],  14)   # 3*2 + 2*1 + 2*3
  expect_equal(s[["N2sum"]], 34)   # 9*2 + 4*1 + 4*3
  expect_equal(s[["Sjd"]],   12)   # 3*2 + 0 + 2*3
  expect_equal(s[["mscore"]], 12/14, tolerance = 1e-12)
  expect_equal(s[["kappa"]],  196/34, tolerance = 1e-12)
  expect_equal(s[["Y_prime"]], 168/34, tolerance = 1e-12)
})

test_that(".mscore_streaming returns empty stats when no reads/CpGs", {
  expect_true(is.na(mHapDMR:::.mscore_streaming(character(0L), 1:5, 1L, 5L)[["mscore"]]))
  expect_true(is.na(mHapDMR:::.mscore_streaming("chr1\t1\t2\t1\t1", integer(0L), 1L, 2L)[["mscore"]]))
})

test_that(".mscore_streaming does not overflow on very deep regions", {
  cpg_pos <- as.integer(seq(100L, by = 2L, length.out = 10L))
  raw <- sprintf("chr1\t100\t118\t%s\t300000000\t+", strrep("1", 10L))
  s <- withCallingHandlers(
    mHapDMR:::.mscore_streaming(raw, cpg_pos, 100L, 118L),
    warning = function(w) stop("unexpected warning: ", conditionMessage(w))
  )
  expect_false(anyNA(s[c("reads", "Nsum", "N2sum", "Sjd")]))
  expect_equal(s[["Nsum"]], 3e9)   # 10 CpGs * 3e8 reads, exceeds 2^31
})


# ── build_bsseq_mscore ────────────────────────────────────────────────────────

test_that("build_bsseq_mscore rejects unnamed rGR_list", {
  rGR <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 200L))
  GenomicRanges::mcols(rGR) <- S4Vectors::DataFrame(
    reads=10, Nsum=30, N2sum=100, Sjd=15,
    mscore=0.5, kappa=9, Y_prime=4.5
  )
  expect_error(build_bsseq_mscore(list(rGR)), "named list")
})


# ── region statistics: reads longer than the CpG lookup margin ────────────────

test_that("mscore_region_stats keeps reads with CpGs beyond the margin", {
  tmp <- tempfile("mhapdmr"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  # CpGs every 200 bp; the region [450, 550] holds the CpG at 500
  cpg_txt <- file.path(tmp, "cpg.txt")
  writeLines(sprintf("chr1\t%d", c(100L, 300L, 500L, 700L)), cpg_txt)
  # read 1 starts 350 bp before the region (CpGs 100, 300, 500; "1" at 500),
  # read 2 (count 2) covers 500 and 700 ("0" at 500)
  mhap_txt <- file.path(tmp, "x.mhap")
  writeLines(c("chr1\t100\t500\t001\t1\t+", "chr1\t500\t700\t01\t2\t+"), mhap_txt)
  cpg  <- Rsamtools::bgzip(cpg_txt, file.path(tmp, "cpg.gz"))
  mhap <- Rsamtools::bgzip(mhap_txt, file.path(tmp, "x.mhap.gz"))
  Rsamtools::indexTabix(cpg, seq = 1L, start = 2L, end = 2L)
  Rsamtools::indexTabix(mhap, seq = 1L, start = 2L, end = 3L)
  rGR <- GenomicRanges::GRanges("chr1", IRanges::IRanges(450L, 550L))
  cols <- c("reads", "Nsum", "N2sum", "Sjd", "mscore", "kappa", "Y_prime")

  s <- mscore_region_stats(mhap, cpg, rGR, verbose = FALSE)
  expect_equal(unlist(as.data.frame(GenomicRanges::mcols(s))[1L, cols]),
               c(reads = 3, Nsum = 3, N2sum = 3, Sjd = 1, mscore = 1 / 3, kappa = 3, Y_prime = 1))
  g <- mscore_region_stats(mhap, cpg, rGR, verbose = FALSE, keep_mhapgr = TRUE)
  expect_equal(as.data.frame(GenomicRanges::mcols(g))[, cols], as.data.frame(GenomicRanges::mcols(s))[, cols])
  expect_equal(length(read_mhap_gr(mhap, cpg, rGR)), 2L)
})
