# mHapDMR Tutorial — Differential M-score Analysis (ESCC Tumor vs Normal)

This tutorial walks through a complete, reproducible run of **mHapDMR** using a
public esophageal squamous-cell carcinoma (ESCC) dataset: three tumor and three
normal samples. It follows the exact workflow in `test_example.R`.

**What the pipeline does.** For every genomic region you supply, mHapDMR reads
the mHap haplotype records, computes a per-sample **M-score**, and tests for a
difference between two groups (here Tumor vs Normal) with `DSS`.

The M-score for region *j*, sample *d* is

```
M_jd = S_jd / T_jd ,   S_jd = Σ_t N_tjd · Z_tjd ,   T_jd = Σ_t N_tjd
```

where *N* is the number of CpG sites on read *t* and *Z ∈ {0,1}* indicates
whether the read carries at least one methylated site. It captures read-level
methylation heterogeneity that a bulk average would miss.

---

## Contents

1. [Input data & where to download it](#1-input-data--where-to-download-it)
2. [Building the design table](#2-building-the-design-table)
3. [Running the pipeline & interpreting the results](#3-running-the-pipeline--interpreting-the-results)

---

## Prerequisites

Install the package and its Bioconductor dependencies once:

```r
BiocManager::install(c("DSS", "bsseq", "GenomicRanges", "IRanges",
                       "S4Vectors", "Rsamtools", "BiocParallel"))
install.packages("data.table")          # for reading the region file
devtools::install_local("mHapDMR")       # or install_github("<user>/mHapDMR")
```

You also need **tabix** (from HTSlib / `samtools`) on your system if you ever
re-index your own files. The example files below are already indexed.

---

## 1. Input data & where to download it

Three kinds of input are required. All of them are **bgzip-compressed and
tabix-indexed** — every `*.gz` has a companion `*.gz.tbi` that must sit next to
it, because the pipeline queries the files by genomic coordinate.

### 1.1 mHap files (6 samples + indexes)

An **mHap** file stores read-level methylation as haplotype strings. Each line
is one collapsed haplotype:

| col | name     | meaning                                                              |
|----:|----------|---------------------------------------------------------------------|
| 1   | `chr`    | chromosome                                                          |
| 2   | `start`  | position of the first CpG covered by the read                      |
| 3   | `end`    | position of the last CpG covered by the read                       |
| 4   | `hapStr` | methylation string, one character per CpG: `1` = methylated, `0` = unmethylated |
| 5   | `count`  | number of identical reads collapsed into this haplotype            |
| 6   | `strand` | `+` / `-` (optional; ignored by the M-score computation)           |

```
chr1  10497  10542  011  11  +
chr1  10497  10542  111  65  +
chr1  10497  10542  101   2  +
```

Download the 3 tumor + 3 normal samples (each `*.mhap.gz` **and** its
`*.mhap.gz.tbi`):

| sample       | group  | mHap file                                                                                       |
|--------------|--------|-------------------------------------------------------------------------------------------------|
| SRX8208812   | Tumor  | http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208812.mhap.gz   |
| SRX8208813   | Tumor  | http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208813.mhap.gz   |
| SRX8208814   | Tumor  | http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208814.mhap.gz   |
| SRX8208802   | Normal | http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208802.mhap.gz  |
| SRX8208803   | Normal | http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208803.mhap.gz  |
| SRX8208804   | Normal | http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208804.mhap.gz  |

The index for each is the same URL with `.tbi` appended, e.g.
`.../tumor/SRX8208812.mhap.gz.tbi`.

> **Note:** on the data server the last index is mistakenly named
> `SRX8208804.mhap.gz.tib` (`.tib`). The correct local filename must be
> `SRX8208804.mhap.gz.tbi` — the download script below fixes this for you.

### 1.2 CpG position file (hg19)

A 3-column, BED-like file listing every CpG in the genome. mHapDMR uses it to
map each haplotype character back to a genomic coordinate.

```
chr1  10469  10470
chr1  10471  10472
chr1  10484  10485
```

- `hg19_CpG.gz` — http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz  (~150 MB)
- `hg19_CpG.gz.tbi` — http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz.tbi

Use the CpG build (hg19) that matches the coordinates of your mHap files.

### 1.3 Regions to test

The regions can be any BED file; here we ship a tab-separated table,
**`ESCC_DMR_subset.txt`** (150 regions), included in this repository under
[`inst/extdata/ESCC_DMR_subset.txt`](inst/extdata/ESCC_DMR_subset.txt).

```
Chr     Start       End         Category
chr19   13151273    13151813    Hyper
chr16   87867030    87867838    Hyper
...
```

| column     | meaning                                                                 |
|------------|-------------------------------------------------------------------------|
| `Chr`      | chromosome                                                               |
| `Start`    | region start, **0-based** (BED convention)                              |
| `End`      | region end                                                               |
| `Category` | expected behaviour label — `Hyper` (hypermethylated in tumor), `Hypo` (hypomethylated in tumor), or `NC` (no change / negative control). 50 regions each. |

`Category` is only there so we can *validate* the result at the end — it is not
used by the model. Because `Start` is 0-based, we add 1 when building the
`GRanges` (Section 3).

### 1.4 Download everything (R)

Run this from the directory where you want to work. It creates `mhap/` and
`ref/` folders and skips files that already exist.

```r
dir.create("mhap", showWarnings = FALSE)
dir.create("ref",  showWarnings = FALSE)

options(timeout = 3600)          # the CpG file is large; raise the 60s default

base <- "http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public"
samples <- data.frame(
  id     = c("SRX8208812", "SRX8208813", "SRX8208814",
             "SRX8208802", "SRX8208803", "SRX8208804"),
  subdir = c("tumor", "tumor", "tumor", "normal", "normal", "normal"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(samples))) {
  for (ext in c(".mhap.gz", ".mhap.gz.tbi")) {
    url  <- sprintf("%s/%s/%s%s", base, samples$subdir[i], samples$id[i], ext)
    dest <- file.path("mhap", paste0(samples$id[i], ext))
    if (!file.exists(dest)) download.file(url, dest, mode = "wb")
  }
}

# CpG positions (hg19) + index
if (!file.exists("ref/hg19_CpG.gz"))
  download.file(paste0("http://bioinformatics.sibcb.ac.cn/dataupload/",
                       "iGenome/CpGs/hg19/hg19_CpG.gz"),
                "ref/hg19_CpG.gz", mode = "wb")
if (!file.exists("ref/hg19_CpG.gz.tbi"))
  download.file(paste0("http://bioinformatics.sibcb.ac.cn/dataupload/",
                       "iGenome/CpGs/hg19/hg19_CpG.gz.tbi"),
                "ref/hg19_CpG.gz.tbi", mode = "wb")
```

<details>
<summary>Shell alternative (wget)</summary>

```bash
mkdir -p mhap ref
base=http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public

for s in tumor/SRX8208812 tumor/SRX8208813 tumor/SRX8208814 \
         normal/SRX8208802 normal/SRX8208803 normal/SRX8208804; do
  id=${s#*/}
  wget -O mhap/$id.mhap.gz     $base/$s.mhap.gz
  wget -O mhap/$id.mhap.gz.tbi $base/$s.mhap.gz.tbi   # note: fixes the .tib typo
done

wget -O ref/hg19_CpG.gz     http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz
wget -O ref/hg19_CpG.gz.tbi http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz.tbi
```
</details>

After downloading, your working directory should look like this:

```
.
├── mhap/
│   ├── SRX8208812.mhap.gz   (+ SRX8208812.mhap.gz.tbi)
│   ├── SRX8208813.mhap.gz   (+ .tbi)
│   ├── SRX8208814.mhap.gz   (+ .tbi)
│   ├── SRX8208802.mhap.gz   (+ .tbi)
│   ├── SRX8208803.mhap.gz   (+ .tbi)
│   └── SRX8208804.mhap.gz   (+ .tbi)
├── ref/
│   ├── hg19_CpG.gz
│   └── hg19_CpG.gz.tbi
└── ESCC_DMR_subset.txt
```

### 1.5 Verify the indexes

Before running anything, confirm each `.gz` has its `.tbi` companion:

```r
library(mHapDMR)
mhap_files <- list.files("mhap", pattern = "\\.mhap\\.gz$", full.names = TRUE)
invisible(lapply(mhap_files, check_tabix_index))
check_tabix_index("ref/hg19_CpG.gz")
```

`check_tabix_index()` returns `TRUE` silently if the index is present and stops
with an informative error otherwise.

---

## 2. Building the design table

The design table tells mHapDMR which group each sample belongs to. Two rules:

1. **`mhap_files` must be a *named* vector.** The names are the sample IDs.
2. **`rownames(design)` must equal `names(mhap_files)`** — this is how a column
   of the M-score matrix is linked to its group. Order matters only in that the
   names must line up; keep them identical.

```r
library(mHapDMR)
library(GenomicRanges)
library(data.table)

# 1. Named vector of mHap files.  Names become the sample IDs.
mhap_files <- c(
  "mhap/SRX8208812.mhap.gz",   # tumor
  "mhap/SRX8208813.mhap.gz",   # tumor
  "mhap/SRX8208814.mhap.gz",   # tumor
  "mhap/SRX8208802.mhap.gz",   # normal
  "mhap/SRX8208803.mhap.gz",   # normal
  "mhap/SRX8208804.mhap.gz"    # normal
)
names(mhap_files) <- sub("\\.mhap\\.gz$", "", basename(mhap_files))
names(mhap_files)
#> [1] "SRX8208812" "SRX8208813" "SRX8208814" "SRX8208802" "SRX8208803" "SRX8208804"

# 2. Design table — one row per sample, row names = sample IDs.
design <- data.frame(
  group     = c(rep("Tumor", 3), rep("Normal", 3)),
  row.names = names(mhap_files)
)
design
#>            group
#> SRX8208812  Tumor
#> SRX8208813  Tumor
#> SRX8208814  Tumor
#> SRX8208802 Normal
#> SRX8208803 Normal
#> SRX8208804 Normal
```

**Choosing `coef`.** The model formula is `~ group`. R turns `group` into a
factor with **alphabetically ordered** levels, so `Normal` becomes the
reference and the coefficient that contrasts Tumor against Normal is called
**`groupTumor`**. That is the string you pass to `coef`. A positive test
statistic then means *higher M-score in Tumor*. If your labels were
`case`/`control`, the coefficient would be `groupcase`, and so on.

The design is not limited to two groups or to a single factor. You can add
covariates and test any coefficient — for example a continuous variable:

```r
# Example only — a continuous covariate instead of a two-group contrast
design <- data.frame(age = c(25, 45, 55, 40, 62, 70),
                     row.names = names(mhap_files))
# ... then formula = ~ age, coef = "age"
```

---

## 3. Running the pipeline & interpreting the results

### 3.1 Define the regions

Read `ESCC_DMR_subset.txt` and build a `GRanges`. Convert the 0-based BED
`Start` to a 1-based `GRanges` start by adding 1:

```r
bed <- fread("ESCC_DMR_subset.txt", header = TRUE)   # Chr, Start, End, Category

rGR <- GRanges(
  seqnames = bed$Chr,
  ranges   = IRanges(start = bed$Start + 1L, end = bed$End)
)
rGR$Category <- bed$Category      # carried along for later validation
length(rGR)
#> [1] 150
```

### 3.2 Run `mhap_dmr()`

`mhap_dmr()` is the single entry point that runs all four internal steps:
per-sample M-score extraction → BSseq assembly → DSS test → coordinate
annotation.

```r
results <- mhap_dmr(
  mhap_files = mhap_files,
  cpg_file   = "ref/hg19_CpG.gz",
  rGR        = rGR,
  design     = design,
  formula    = ~ group,
  coef       = "groupTumor",
  min_reads  = 1L,
  margin     = 150L,
  smoothing  = FALSE,
  BPPARAM    = BiocParallel::MulticoreParam(workers = 4L),
  verbose    = TRUE
)
```

Key arguments:

| argument    | value here | meaning                                                                                 |
|-------------|-----------|------------------------------------------------------------------------------------------|
| `formula`   | `~ group` | model formula passed to DSS                                                               |
| `coef`      | `"groupTumor"` | coefficient to test (see Section 2)                                                  |
| `min_reads` | `1L`      | a region is kept only if **every** sample has ≥ this read depth there                     |
| `margin`    | `150L`    | extra bp padded around each region when looking up CpG positions                          |
| `smoothing` | `FALSE`   | passed to DSS; keep `FALSE` for region-level (block) tests                                |
| `BPPARAM`   | 4 workers | sample-level parallelism. On Windows use `SnowParam(workers = 4L)`; `MulticoreParam` falls back to serial there. Reduce workers if memory is tight. |
| `verbose`   | `TRUE`    | print per-region / per-step progress                                                     |

You will see progress messages and finally a note such as
`Regions retained: 150 / 150` (or fewer, if some regions had no coverage in a
sample under `min_reads`).

### 3.3 What `results` contains

`mhap_dmr()` returns a named list:

| element        | class        | description                                                            |
|----------------|--------------|------------------------------------------------------------------------|
| `dml_result`   | `data.frame` | the DSS test result, annotated with region coordinates, sorted by p-value |
| `m_summary`    | `matrix`     | group-mean M-score per region (one column per group)                   |
| `bsseq`        | `BSseq`      | the object handed to DSS (κ as coverage, Y′ as methylation)            |
| `region_index` | `data.frame` | maps the DSS `(chr, pos)` key back to `region_id`/`start`/`end`        |
| `rGR_list`     | named list   | per-sample `GRanges` carrying the raw M-score statistics               |

The main table, `dml_result`, has these columns:

| column      | meaning                                                                 |
|-------------|-------------------------------------------------------------------------|
| `region_id` | region identifier, `chr:start-end` (1-based)                            |
| `chr`, `start`, `end` | original region coordinates                                   |
| `pos`       | region midpoint — the internal DSS locus key                            |
| `stat`      | DSS test statistic; **sign = direction** (positive ⇒ higher in Tumor)   |
| `pvals`     | raw p-value                                                             |
| `fdrs`      | Benjamini–Hochberg FDR (use this for significance)                     |

### 3.4 Assemble a readable results table

Attach the group-mean M-scores (and, for this example, the original `Category`
label) to the test result. Rows of `m_summary` are keyed by `region_id`, so the
join is a direct row lookup — this is exactly the pattern in `test_example.R`:

```r
res <- results$dml_result

# group-mean M-scores (columns: Tumor, Normal)
res <- data.frame(res, results$m_summary[res$region_id, ], row.names = NULL)

# effect size on the M-score scale
res$delta <- res$Tumor - res$Normal

# bring back the expected-behaviour label for validation
bed$region_id <- paste0(bed$Chr, ":", bed$Start + 1L, "-", bed$End)
res$Category  <- bed$Category[match(res$region_id, bed$region_id)]

head(res[, c("region_id", "Category", "Tumor", "Normal", "delta",
             "stat", "pvals", "fdrs")], 10)
```

Each row is one region. Read it as: `Tumor` and `Normal` are the mean M-scores
in each group; `delta = Tumor − Normal` is the effect size; `stat` gives the
signed significance and should share the sign of `delta`; `fdrs` is what you
threshold on.

### 3.5 Call significant regions

```r
sig <- subset(res, fdrs < 0.05)
nrow(sig)

# hypermethylated in tumor (gain of methylation)
subset(sig, delta > 0)

# hypomethylated in tumor (loss of methylation)
subset(sig, delta < 0)
```

To keep the regions as a `GRanges` (e.g. for annotation or export), use the
helper:

```r
dml_gr <- annotate_dml_result(results$dml_result, results$region_index,
                              as_granges = TRUE)
dml_gr[dml_gr$fdrs < 0.05]
```

### 3.6 Validate against the `Category` labels

Because our example regions come pre-labelled, we can confirm the pipeline
behaves as expected. A correct run should show:

- **`Hyper`** regions → `delta > 0` (Tumor > Normal) and mostly significant,
- **`Hypo`** regions → `delta < 0` (Tumor < Normal) and mostly significant,
- **`NC`** regions → `delta ≈ 0` and mostly **non**-significant.

```r
res$sig <- res$fdrs < 0.05

# median effect size per category — sign should match the label
aggregate(delta ~ Category, data = res, FUN = median)

# significant vs. non-significant, per category
table(Category = res$Category, significant = res$sig)
```

If `Hyper` sits at a positive median `delta`, `Hypo` at a negative one, and `NC`
near zero with few significant calls, the M-score DMR analysis is working as
intended. From there, swap in your own regions and design table to run the same
analysis on any two-group (or covariate) comparison.

---

## Appendix — running the steps manually

`mhap_dmr()` is a thin wrapper. If you want the intermediate objects, the same
analysis unrolls to:

```r
# 1. per-sample M-score statistics over all regions
rGR_list <- lapply(mhap_files, mscore_region_stats,
                   cpg_file = "ref/hg19_CpG.gz", rGR = rGR)
names(rGR_list) <- names(mhap_files)

# 2. assemble a BSseq object via the kappa / Y' substitution
built <- build_bsseq_mscore(rGR_list, min_reads = 1L)

# 3. differential test with DSS
dml_raw <- run_dml_test(built$bsseq, design,
                        formula = ~ group, coef = "groupTumor")

# 4. restore region coordinates
dml_result <- annotate_dml_result(dml_raw, built$region_index)
```

See `?mhap_dmr`, `?mscore_region_stats`, and `?build_bsseq_mscore` for full
argument documentation.
