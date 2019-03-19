#===============================================================================
# tssenrich.py
#===============================================================================

"""Calculate TSS enrichment for ATAC-seq data"""




# Imports ======================================================================

import argparse
import gzip
import itertools
import os
import os.path
import pybedtools
import shutil
import subprocess
import tempfile

from functools import partial
from multiprocessing import Pool



# Constants ====================================================================

SAMTOOLS_PATH = os.environ.get('SAMTOOLS_PATH', shutil.which('samtools'))

ENCODE_STANDARDS = '''ENCODE standards:
| Genome | Concerning | Acceptable | Ideal |
| ------ | ---------- | ---------- | ----- |
| hg19   | < 6        | 6 - 10     | > 10  |
| hg38   | < 5        | 5 - 7      | > 7   |
| mm10   | < 10       | 10 - 15    | > 15  |
'''




# Exceptions ===================================================================

class Error(Exception):
   """Base class for other exceptions"""

   pass


class MissingSAMToolsError(Error):
    """Missing samtools error"""

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


def generate_tss_flanks(
    tss,
    flank_distance: int = 1_000,
    flank_size: int = 100
):
    """Generate coordinates of TSS flanks

    Parameters
    ----------
    tss
        an iterable containing coordinates of TSS's
    flank_distance : int
        distance from tss of outer ends of flanks
    flank_size : int
        size of flanks (for determining average depth)

    Yields
    ------
    tuple
        The coordinates of a TSS flank and the corresponding TSS, in the form:
        chrom, flank_start, flank_end, tss_pos
    """
    for chrom, pos in tss:
        if pos >= flank_distance:
            yield (
                chrom,
                pos - flank_distance,
                pos - flank_distance + flank_size,
                pos
            )
            yield chrom, pos - 1, pos, pos
            yield (
                chrom,
                pos + flank_distance - flank_size,
                pos + flank_distance,
                pos
            )


def tss_flanks_bed_str(flanks):
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


def tss_flanks_bed_tool(flanks_str: str, temp_file_dir=None):
    """A BedTool representing the TSS flanks

    Parameters
    ----------
    flanks_str
        string giving input BED file
    temp_file_dir
        directory to use for temporary files

    Returns
    -------
    BedTool
        the TSS flanks
    """

    pybedtools.set_tempdir(temp_file_dir if temp_file_dir else '/tmp')
    return pybedtools.BedTool(flanks_str, from_string=True).sort()


def samtools_bedcov(
    bed_file_path,
    bam_file_path,
    memory_gb: float = 1.0,
    threads: int = 1,
    mapping_quality: int = 0,
    samtools_path: str = SAMTOOLS_PATH,
    log_file_path=None,
    temp_file_dir=None
):
    """Apply samtools bedcov to a bed file & a bam file

    Parameters
    ----------
    bed_file_path : str
        path to a BED file
    bam_file_path : str
        path to a BAM file
    memory_gb : float
        memory limit in gigabytes
    threads : int
        number of threads to use for sorting
    mapping_quality : int
        ignore reads with mapping quality below this value [0]
    samtools_path : str
        path to the samtools executable
    log_file_path
        path to a log file
    temp_file_dir
        directory to use for temporary files
    
    Returns
    -------
    bytes
        the output of samtools bedcov
    """

    if not samtools_path:
        raise MissingSAMToolsError(
            '''samtools was not found! Please provide the `samtools_path`
            parameter, or set the `SAMTOOLS_PATH` environment variable, or make
            sure `samtools` is installed and can be found via the `PATH`
            environment variable.
            '''
        )
    
    log_file = open(log_file_path, 'wb') if log_file_path else None
    with tempfile.TemporaryDirectory(dir=temp_file_dir) as temp_dir_name:
        sorted_path = os.path.join(temp_dir_name, 'sorted.bam')
        with open(sorted_path, 'wb') as temp_sorted:
            subprocess.run(
                (
                    samtools_path, 'sort',
                    '-m', '{}M'.format(int(1024 / threads * memory_gb)),
                    '-@', str(threads),
                    bam_file_path
                ),
                stdout=temp_sorted,
                stderr=log_file
            )
            subprocess.run((samtools_path, 'index', sorted_path))
        with subprocess.Popen(
            (
                samtools_path, 'bedcov',
                '-Q', str(mapping_quality),
                bed_file_path,
                sorted_path
            ),
            stdout=subprocess.PIPE,
            stderr=log_file
        ) as bedcov:
            return bedcov.communicate()[0]
        if log_file_path:
            log_file.close()


def generate_coverage_values(bedcov: bytes):
    """Generate coverage values from the output of samtools bedcov

    Parameters
    ----------
    bedcov : bytes
        output from samtools bedcov
    
    Yields
    ------
    tuple
        (tss_center_depth, flank_depth) per TSS
    """

    for tss, intervals in itertools.groupby(
        sorted(
            (
                (chrom, int(start), int(end), int(tss), int(cov))
                for chrom, start, end, tss, cov in (
                    line.split() for line in bedcov.decode().splitlines()
                )
            ),
            key=lambda interval: (interval[0], interval[3])
        ),
        key=lambda interval: (interval[0], interval[3])
    ):
        lower_flank_cov, tss_cov, upper_flank_cov = (
            interval[4] for interval in sorted(set(intervals))
        )
        yield tss_cov, (lower_flank_cov + upper_flank_cov) / 200


def calculate_enrichment(coverage_values):
    """Calculate TSS enrichment value for a dataset

    Parameters
    ----------
    coverage_values
        iterable of tuples (tss_center_depth, flank_depth) per TSS
    
    Returns
    -------
    float
        the TSS enrichment value
    """

    tss_depth, flank_depth = (sum(z) for z in zip(*coverage_values))
    return tss_depth / flank_depth


def tss_enrichment(
    bam_file_path,
    genome='hg38',
    memory_gb: float = 1.0,
    threads: int = 1,
    mapping_quality: int = 0,
    samtools_path: str = SAMTOOLS_PATH,
    log_file_path=None,
    temp_file_dir=None,
    flank_distance: int = 1_000,
    flank_size: int = 100
):
    """Calculate TSS enrichment from ATAC-seq data
    
    Parameters
    ----------
    bam_file_path : str
        path to input BAM file
    genome : str
        genome build to use (must be 'hg38' or 'hg19') [hg38]
    memory_gb : float
        memory limit in gigabytes [1.0]
    threads : int
        number of threads to use for sorting [1]
    mapping_quality : int
        ignore reads with mapping quality below this value [0]
    samtools_path : str
        path to the samtools executable
    log_file_path
        path to a log file
    
    Returns
    -------
    float
        the TSS enrichment value
    """
    
    tss_flanks = tss_flanks_bed_tool(
        tss_flanks_bed_str(
            generate_tss_flanks(
                generate_tss(genome=genome),
                flank_distance=flank_distance,
                flank_size=flank_size
            )
        ),
        temp_file_dir=temp_file_dir
    )
    return calculate_enrichment(
        generate_coverage_values(
            samtools_bedcov(
                tss_flanks.fn,
                bam_file_path,
                memory_gb=memory_gb,
                threads=threads,
                mapping_quality=mapping_quality,
                samtools_path=samtools_path,
                log_file_path=log_file_path,
                temp_file_dir=temp_file_dir
            )
        )
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='calculate TSS enrichment for ATAC-seq data',
        epilog=ENCODE_STANDARDS,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'bam',
        metavar='<path/to/file.bam>',
        nargs='+',
        help='Path to input BAM file'
    )
    parser.add_argument(
        '--genome',
        choices=('hg38', 'hg19', 'mm10'),
        default='hg38',
        help='genome build [hg38]'
    )
    parser.add_argument(
        '--memory',
        metavar='<float>',
        type=float,
        default=1.0,
        help='memory limit in GB [2]'
    )
    parser.add_argument(
        '--processes',
        metavar='<int>',
        type=int,
        default=1,
        help='number of processes/threads to use [1]'
    )
    parser.add_argument(
        '--mapping-quality',
        metavar='<int>',
        type=int,
        default=0,
        help='ignore reads with mapping quality below the given value [0]'
    )
    parser.add_argument(
        '--samtools-path',
        metavar='<path/to/samtools>',
        default=SAMTOOLS_PATH,
        help=f'path to an alternate samtools executable [{SAMTOOLS_PATH}]'
    )
    parser.add_argument(
        '--log',
        metavar='<path/to/log.txt>',
        help='path to log file'
    )
    parser.add_argument(
        '--tmp-dir',
        metavar='<temp/file/dir/>',
        help='directory to use for temporary files'
    )
    parser.add_argument(
        '--flank-distance',
        metavar='<int>',
        default=1_000,
        help='distance from tss of outer ends of flanks [1000]'
    )
    parser.add_argument(
        '--flank-size',
        metavar='<int>',
        default=100,
        help='size of flanks (for determining average depth) [100]'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    n_bam = len(args.bam)
    with Pool(processes=min(args.processes, n_bam)) as pool:
        values = pool.map(
            partial(
                tss_enrichment,
                genome=args.genome,
                memory_gb=args.memory / n_bam,
                threads=int(args.processes / n_bam),
                mapping_quality=args.mapping_quality,
                samtools_path=args.samtools_path,
                log_file_path=args.log,
                temp_file_dir=args.tmp_dir,
                flank_distance=args.flank_distance,
                flank_size=args.flank_size
            ),
            args.bam
        )
    for value in values:
        print(value)
