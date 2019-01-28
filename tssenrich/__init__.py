"""refFlat data source:
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
"""

from tssenrich.tssenrich import (
    generate_tss, generate_flanks, flanks_bed_str, flanks_bed_tool,
    samtools_bedcov
)