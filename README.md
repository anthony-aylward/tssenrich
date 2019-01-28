# tssenrich

Calculate TSS enrichment for ATAC-seq data

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
usage: tssenrich [-h] [--genome {hg38,hg19}] [--memory <float>]
                 [--processes <int>] [--mapping-quality <int>]
                 [--samtools-path <path/to/samtools>]
                 [--log <path/to/log.txt>]
                 <path/to/file.bam>

calculate TSS enrichment for ATAC-seq data

positional arguments:
  <path/to/file.bam>    Path to input BAM file

optional arguments:
  -h, --help            show this help message and exit
  --genome {hg38,hg19}  genome build [hg38]
  --memory <float>      memory limit in GB [5]
  --processes <int>     number of processes/threads to use [1]
  --mapping-quality <int>
                        ignore reads with mapping quality below the given value [0]
  --samtools-path <path/to/samtools>
                        path to an alternate samtools executable
  --log <path/to/log.txt>
                        path to log file

ENCODE standards:
| Genome | Concerning | Acceptable | Ideal |
| ------ | ---------- | ---------- | ----- |
| hg19   | < 6        | 6 - 10     | > 10  |
| hg38   | < 5        | 5 - 7      | > 7   |
```

## Example
```
tssenrich --genome hg19 --log log.txt --memory 2 --processes 2 example.bam > score.txt
```
