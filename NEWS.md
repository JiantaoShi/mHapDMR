# mHapDMR 0.1.1

* `mscore_region_stats()`, `read_mhap_gr()` and `mhap_dmr()` no longer drop
  reads with CpGs farther than `margin` bp outside the region. The CpG
  positions were looked up only within `margin` (150 bp) of the region, so
  such a read failed the check that its haplotype covers as many CpGs as the
  CpG file has between its start and end. The lookup now covers the full span
  of every read overlapping the region; `margin` no longer changes the
  results. On 150 ESCC DMRs (reads up to 233 bp), 25 regions had lost 2 to 4
  reads each (M-score change up to 0.02).

# mHapDMR 0.1.0

* First release.
