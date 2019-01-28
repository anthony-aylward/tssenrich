#===============================================================================
# tssenrich.py
#===============================================================================

"""Calculate TSS enrichment for ATAC-seq data"""




# Imports ======================================================================

import gzip
import os.path
import pybedtools
import pysam




# Exceptions ===================================================================

class Error(Exception):
   """Base class for other exceptions"""

   pass


class MissingSAMToolsError(Error):
    """Missing samtools error"""

    pass


class MemoryLimitError(Error):
    """Memory limit error"""
    
    pass




# Functions ====================================================================

def generate_tss(genome='hg38'):
    """A generator yielding the coordinates of transcription start sites

    Parameters
    ----------
    genome : str
        Genome build from which to draw TSS's. Must be either 'hg38' or 'hg19'.
    """

    with gzip.open(
        os.path.join(os.path.dirname(__file__), f'{genome}.refFlat.txt.gz'),
        'rt'
    ) as f:
        for line in f:
            gene_name, name, chrom, strand, tss, *rest = line.split()
            yield chrom, int(tss)


def generate_flanks(tss):
    for chrom, pos in tss:
        if pos >= 1_000:
            yield chrom, pos - 1_000, pos - 900
        yield chrom, pos + 900, pos + 1_000


def flanks_bed_tool(flanks):
    """A BedTool representing the TSS flanks

    Parameters
    ----------
    flanks
        iterable giving the coordinates of the flanks
    
    Returns
    -------
    BedTool
        the TSS flanks
    """

    return pybedtools.BedTool(
        '\n'.join(' '.join(str(flank) for flank in flanks)) + '\n',
        from_string=True
    ).sort()
