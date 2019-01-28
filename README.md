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

## Example
```
tssenrich --genome hg19 --log log.txt --memory 2 --processes 2 example.bam > score.txt
```