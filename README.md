# mHapDMR

Differential M-score analysis from mHap files using DSS.

## Installation

```r
BiocManager::install(c("DSS", "GenomicRanges", "IRanges",
                       "S4Vectors", "Rsamtools", "BiocParallel"))
devtools::install_local("mHapDMR")
```

## Quick example

```r
library(mHapDMR)
library(GenomicRanges)

rGR <- GRanges("chr1", IRanges(c(1000L, 5000L), c(1200L, 5300L)))

result <- mhap_dmr(
  mhap_files = c(C1="ctrl1.mhap.gz", C2="ctrl2.mhap.gz",
                 T1="tumor1.mhap.gz", T2="tumor2.mhap.gz"),
  cpg_file   = "hg19_CpG.gz",
  rGR        = rGR,
  design     = data.frame(group = c("Control","Control","Tumor","Tumor")),
  formula    = ~ group,
  coef       = "groupTumor"
)

head(result$dml_result)
```

## Key objects

| Object | Class | Description |
|---|---|---|
| `rGR` | `GRanges` | Analysis regions |
| `gr` | `GRanges` | Single analysis region |
| `mHapGR` | `GRanges` | mHap reads for one region (one row per haplotype) |
| `rGR_list` | named list of `GRanges` | Per-sample rGR with M-score metadata |

## Key functions

| Function | Description |
|---|---|
| `mhap_dmr()` | Full pipeline entry point |
| `read_mhap_gr()` | Read mHap records for one region → mHapGR |
| `mscore_region_stats()` | M-score stats for all regions in rGR |
| `build_bsseq_mscore()` | Build BSseq with κ/Y' substitution |
| `run_dml_test()` | DSS DML test wrapper |
| `check_tabix_index()` | Verify tabix index exists |

## M-score definition

$$M_{jd} = \frac{\sum_t N_{tjd} Z_{tjd}}{\sum_t N_{tjd}}$$

where $N_{tjd}$ = CpG count on read $t$, $Z_{tjd}$ = 1 if read has ≥1 methylated site.
