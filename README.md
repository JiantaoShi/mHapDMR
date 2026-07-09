# <center> Differential methylation analysis with mHpDMR <center/>

This tutorial walks through a complete, reproducible run of **mHapDMR** using a public esophageal squamous-cell carcinoma (ESCC) dataset: three tumor and three normal samples. For every genomic region you supply, mHapDMR reads the [mHap](https://jiantaoshi.github.io/mHap/) records, computes a per-sample M-score, and tests for a difference between two groups (here Tumor vs Normal) with [DSS](https://www.bioconductor.org/packages/release/bioc/html/DSS.html).

## Installation

Install the package and its Bioconductor dependencies once:

```R
BiocManager::install(c("DSS", "bsseq", "GenomicRanges", "IRanges", "S4Vectors", "Rsamtools", "BiocParallel"))
remotes::install_github("JiantaoShi/mHapDMR")
```

You also need tabix (from HTSlib / samtools) on your system if you ever re-index your own files. The example files below are already indexed.


```R
library(mHapDMR)
library(GenomicRanges)
library(data.table)
```

## Input

Three kinds of input are required. All of them are bgzip-compressed and tabix-indexed — every *.gz has a companion *.gz.tbi that must sit next to it, because the pipeline queries the files by genomic coordinate.

### mHap files

An [mHap](https://jiantaoshi.github.io/mHap/) file stores read-level methylation as haplotype strings. The mHap files used in this tutorial can be downloaded below:

- [SRX8208812](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208812.mhap.gz)
- [SRX8208813](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208813.mhap.gz)
- [SRX8208814](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/tumor/SRX8208814.mhap.gz)
- [SRX8208802](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208802.mhap.gz)
- [SRX8208803](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208803.mhap.gz)
- [SRX8208804](http://bioinformatics.sibcb.ac.cn/dataupload/cancermhaps/mHap/public/normal/SRX8208804.mhap.gz)

The index for each is the same URL with .tbi appended, e.g. .../tumor/SRX8208812.mhap.gz.tbi.

### CpG position file
A 3-column, BED-like file listing every CpG in the genome. mHapDMR uses it to map each haplotype character back to a genomic coordinate.

- [hg19_CpG.gz](http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz)
- [hg19_CpG.gz.tbi](http://bioinformatics.sibcb.ac.cn/dataupload/iGenome/CpGs/hg19/hg19_CpG.gz.tbi)

Use the [annotation files](https://jiantaoshi.github.io/mHap/AnnotationFiles.html) that matches the coordinates of your mHap files.

### Regions to test

The regions can be any BED file; here we ship a tab-separated table, `ESCC_DMR_subset.txt` (150 regions), included in this repository under `inst/extdata/ESCC_DMR_subset.txt`.


```R
head(read.table('mHapDMR/inst/extdata/ESCC_DMR_subset.txt'))

```


<table class="dataframe">
<caption>A data.frame: 6 × 4</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>Chr  </td><td>Start    </td><td>End      </td><td>Category</td></tr>
	<tr><th scope=row>2</th><td>chr19</td><td>13151273 </td><td>13151813 </td><td>Hyper   </td></tr>
	<tr><th scope=row>3</th><td>chr16</td><td>87867030 </td><td>87867838 </td><td>Hyper   </td></tr>
	<tr><th scope=row>4</th><td>chr19</td><td>39203688 </td><td>39203964 </td><td>Hyper   </td></tr>
	<tr><th scope=row>5</th><td>chr18</td><td>12010537 </td><td>12011054 </td><td>Hyper   </td></tr>
	<tr><th scope=row>6</th><td>chr3 </td><td>128068881</td><td>128069254</td><td>Hyper   </td></tr>
</tbody>
</table>



'Category' is only there so we can validate the result at the end — it is not used by the model.

## Building the design table


```R
# 1. Named vector of mHap files.  Names become the sample IDs.
mhap_files <- c(
  "ESCC/mHap/SRX8208812.mhap.gz",   # tumor
  "ESCC/mHap/SRX8208813.mhap.gz",   # tumor
  "ESCC/mHap/SRX8208814.mhap.gz",   # tumor
  "ESCC/mHap/SRX8208802.mhap.gz",   # normal
  "ESCC/mHap/SRX8208803.mhap.gz",   # normal
  "ESCC/mHap/SRX8208804.mhap.gz"    # normal
)
names(mhap_files) <- sub("\\.mhap\\.gz$", "", basename(mhap_files))

# 2. Design table — one row per sample, row names = sample IDs.
design <- data.frame(
  group     = c(rep("Tumor", 3), rep("Normal", 3)),
  row.names = names(mhap_files)
)
design
```


<table class="dataframe">
<caption>A data.frame: 6 × 1</caption>
<thead>
	<tr><th></th><th scope=col>group</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>SRX8208812</th><td>Tumor </td></tr>
	<tr><th scope=row>SRX8208813</th><td>Tumor </td></tr>
	<tr><th scope=row>SRX8208814</th><td>Tumor </td></tr>
	<tr><th scope=row>SRX8208802</th><td>Normal</td></tr>
	<tr><th scope=row>SRX8208803</th><td>Normal</td></tr>
	<tr><th scope=row>SRX8208804</th><td>Normal</td></tr>
</tbody>
</table>



### Choosing coef

The model formula is `~ group`. R turns group into a factor with alphabetically ordered levels, so Normal becomes the reference and the coefficient that contrasts Tumor against Normal is called groupTumor. That is the string you pass to coef. A positive test statistic then means higher M-score in Tumor. If your labels were case/control, the coefficient would be groupcase, and so on. The design is not limited to two groups or to a single factor. You can add covariates and test any coefficient.

## Differential methylation analysis

### loading test regions


```R
path <- system.file("extdata", "ESCC_DMR_subset.txt", package = "mHapDMR")
region <- fread(path, header = TRUE)
rGR <- GRanges(seqnames = region$Chr, ranges = IRanges(start = region$Start + 1, end = region$End))
rGR$Category <- region$Category
names(rGR) <- paste0('r_', 1:length(rGR))
rGR
```


    GRanges object with 150 ranges and 1 metadata column:
            seqnames              ranges strand |    Category
               <Rle>           <IRanges>  <Rle> | <character>
        r_1    chr19   13151274-13151813      * |       Hyper
        r_2    chr16   87867031-87867838      * |       Hyper
        r_3    chr19   39203689-39203964      * |       Hyper
        r_4    chr18   12010538-12011054      * |       Hyper
        r_5     chr3 128068882-128069254      * |       Hyper
        ...      ...                 ...    ... .         ...
      r_146     chr4 123715987-123717704      * |          NC
      r_147    chr21   43800537-43800923      * |          NC
      r_148     chr7   72843765-72843852      * |          NC
      r_149     chr5 127219488-127220302      * |          NC
      r_150     chr1 186546942-186547461      * |          NC
      -------
      seqinfo: 23 sequences from an unspecified genome; no seqlengths


### main function


```R
results <- mhap_dmr(
  mhap_files = mhap_files,
  cpg_file   = "hg19_CpG.gz",
  rGR        = rGR,
  design     = design,
  formula    = ~ group,
  coef       = "groupTumor",
  min_reads  = 1L,
  margin     = 150L,
  smoothing  = FALSE,
  BPPARAM    = BiocParallel::MulticoreParam(workers = 1L), ## set up workers to speed up
  verbose    = TRUE
)
names(results)
```

    === Step 1: computing M-score statistics per sample ===
    
    [100/150] chr14:22359193-22359478
    
    [100/150] chr14:22359193-22359478
    
    [100/150] chr14:22359193-22359478
    
    [100/150] chr14:22359193-22359478
    
    [100/150] chr14:22359193-22359478
    
    [100/150] chr14:22359193-22359478
    
    === Step 2: building BSseq object ===
    
    Regions retained: 150 / 150
    
    === Step 3: differential M-score test (DSS) ===
    


    Fitting DML model for CpG site: 


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'dml_result'</li><li>'bsseq'</li><li>'region_index'</li><li>'rGR_list'</li><li>'m_summary'</li></ol>



### results


```R
res <- results$dml_result
head(res)
```


<table class="dataframe">
<caption>A data.frame: 6 × 8</caption>
<thead>
	<tr><th></th><th scope=col>region_id</th><th scope=col>chr</th><th scope=col>start</th><th scope=col>end</th><th scope=col>pos</th><th scope=col>stat</th><th scope=col>pvals</th><th scope=col>fdrs</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>r_55</td><td>chr3 </td><td>179716475</td><td>179717697</td><td>179717086</td><td>-12.445317</td><td>1.482951e-35</td><td>2.224426e-33</td></tr>
	<tr><th scope=row>2</th><td>r_53</td><td>chr7 </td><td>153355141</td><td>153355744</td><td>153355442</td><td>-10.445725</td><td>1.532791e-25</td><td>1.149593e-23</td></tr>
	<tr><th scope=row>3</th><td>r_75</td><td>chr3 </td><td>179949384</td><td>179950193</td><td>179949788</td><td>-10.305598</td><td>6.647686e-25</td><td>3.323843e-23</td></tr>
	<tr><th scope=row>4</th><td>r_67</td><td>chr14</td><td> 49875676</td><td> 49877568</td><td> 49876622</td><td>-10.276428</td><td>9.000187e-25</td><td>3.375070e-23</td></tr>
	<tr><th scope=row>5</th><td>r_99</td><td>chr6 </td><td>127833825</td><td>127835103</td><td>127834464</td><td>-10.157850</td><td>3.057457e-24</td><td>9.172370e-23</td></tr>
	<tr><th scope=row>6</th><td>r_64</td><td>chr8 </td><td>138371138</td><td>138371828</td><td>138371483</td><td> -9.891427</td><td>4.535176e-23</td><td>1.133794e-21</td></tr>
</tbody>
</table>




```R
m_summary <- results$m_summary
head(m_summary)
```


<table class="dataframe">
<caption>A matrix: 6 × 2 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Tumor</th><th scope=col>Normal</th></tr>
</thead>
<tbody>
	<tr><th scope=row>chr19:13151274-13151813</th><td>0.9428430</td><td>0.5852567</td></tr>
	<tr><th scope=row>chr16:87867031-87867838</th><td>0.9939556</td><td>0.7107815</td></tr>
	<tr><th scope=row>chr19:39203689-39203964</th><td>0.9383133</td><td>0.5466225</td></tr>
	<tr><th scope=row>chr18:12010538-12011054</th><td>0.9556178</td><td>0.6106748</td></tr>
	<tr><th scope=row>chr3:128068882-128069254</th><td>0.9856147</td><td>0.7172536</td></tr>
	<tr><th scope=row>chr20:35827353-35827866</th><td>0.9903759</td><td>0.6768616</td></tr>
</tbody>
</table>




```R

```
