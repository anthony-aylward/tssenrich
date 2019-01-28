"""Calculate TSS enrichment for ATAC-seq data

refFlat data source:
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/

Functions
---------
    generate_tss
    generate_tss_flanks
    tss_flanks_bed_str
    tss_flanks_bed_tool
    samtools_bedcov
    generate_coverage_values
"""

from tssenrich.tssenrich import (
    generate_tss, generate_tss_flanks, tss_flanks_bed_str, tss_flanks_bed_tool,
    samtools_bedcov, generate_coverage_values
)