# tssenrich

Calculate [TSS enrichment](https://www.encodeproject.org/data-standards/terms/#enrichment) for ATAC-seq data. (see also the [ENCODE standards](https://www.encodeproject.org/atac-seq/))

## Installation
```
pip3 install tssenrich
```
or
```
pip3 install --user tssenrich
```

## Usage
```
usage: tssenrich [-h] [--genome {hg38,hg19,mm10}] [--memory <float>]
                 [--processes <int>] [--mapping-quality <int>]
                 [--samtools-path <path/to/samtools>]
                 [--log <path/to/log.txt>] [--tmp-dir <temp/file/dir/>]
                 [--flank-distance <int>] [--flank-size <int>]
                 <path/to/file.bam> [<path/to/file.bam> ...]

calculate TSS enrichment for ATAC-seq data

positional arguments:
  <path/to/file.bam>    Path to input BAM file

optional arguments:
  -h, --help            show this help message and exit
  --genome {hg38,hg19,mm10}
                        genome build [hg38]
  --memory <float>      memory limit in GB [2]
  --processes <int>     number of processes/threads to use [1]
  --mapping-quality <int>
                        ignore reads with mapping quality below the given value [0]
  --samtools-path <path/to/samtools>
                        path to an alternate samtools executable [/usr/local/bin/samtools]
  --log <path/to/log.txt>
                        path to log file
  --tmp-dir <temp/file/dir/>
                        directory to use for temporary files
  --flank-distance <int>
                        distance from tss of outer ends of flanks [1000]
  --flank-size <int>    size of flanks (for determining average depth) [100]

ENCODE standards:
| Genome | Concerning | Acceptable | Ideal |
| ------ | ---------- | ---------- | ----- |
| hg19   | < 6        | 6 - 10     | > 10  |
| hg38   | < 5        | 5 - 7      | > 7   |
| mm10   | < 10       | 10 - 15    | > 15  |
```

## Examples
```sh
tssenrich example.bam > score.txt
```

```sh
tssenrich --genome hg19 --log log.txt --memory 2 --processes 2 example1.bam example2.bam > scores.txt
```
