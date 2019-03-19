"""Calculate TSS enrichment for ATAC-seq data

definition of TSS enrichment score:
https://www.encodeproject.org/data-standards/terms/

ENCODE standards:
https://www.encodeproject.org/atac-seq/

refFlat data source:
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/

Functions
---------
    generate_tss
    generate_tss_flanks
    tss_flanks_bed_str
    tss_flanks_bed_tool
    samtools_bedcov
    generate_coverage_values
    calculate_enrichment
    tss_enrichment
"""

from tssenrich.tssenrich import (
    generate_tss, generate_tss_flanks, tss_flanks_bed_str, tss_flanks_bed_tool,
    samtools_bedcov, generate_coverage_values, calculate_enrichment,
    tss_enrichment
)