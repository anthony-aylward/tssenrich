#===============================================================================
# tssenrich.py
#===============================================================================

"""Calculate TSS enrichment for ATAC-seq data"""




# Imports ======================================================================

import gzip
import os.path
import pybedtools
import pysam
import tempfile




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
    
    Yields
    ------
    tuple
        the chromosome (str) and position (int) of a TSS 
    """

    with gzip.open(
        os.path.join(os.path.dirname(__file__), f'{genome}.refFlat.txt.gz'),
        'rt'
    ) as f:
        for line in f:
            gene_name, name, chrom, strand, tss, *rest = line.split()
            yield chrom, int(tss)


def generate_flanks(tss):
    """Generate coordinates of TSS flanks

    Parameters
    ----------
    tss
        an iterable containing coordinates of TSS's
    
    Yields
    ------
    tuple
        The coordinates of a TSS flank and the corresponding TSS, in the form:
        chrom, flank_start, flank_end, tss_pos
    """
    for chrom, pos in tss:
        if pos >= 1_000:
            yield chrom, pos - 1_000, pos - 900, pos
        yield chrom, pos + 900, pos + 1_000, pos


def flanks_bed_str(flanks):
    """Create a string object containing TSS flanks in BED format

    Parameters
    ----------
    flanks
        coordinates of TSS flanks
    
    Returns
    -------
    str
        a BED file containing TSS flanks
    """
    return '\n'.join(
        '\t'.join(str(coord) for coord in flank) for flank in flanks
    ) + '\n'


def flanks_bed_tool(flanks_str: str):
    """A BedTool representing the TSS flanks

    Parameters
    ----------
    flanks_str
        string giving input BED file
    
    Returns
    -------
    BedTool
        the TSS flanks
    """

    return pybedtools.BedTool(flanks_str, from_string=True).sort()


def bedcov(flanks_path, bam_file_path, mapping_quality: int = 0):
    with tempfile.TemporaryDirectory() as temp_dir_name:
        sorted_path = os.path.join(temp_dir_name, 'sorted.bam')
        pysam.sort('-o', sorted_path, bam_file_path)
        return pysam.bedcov(
            '-Q', str(mapping_quality),
            flanks_path,
            sorted_path
        )
