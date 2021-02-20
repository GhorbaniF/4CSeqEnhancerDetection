# 4C-seq-enhancer-detection

Discriotion of pipeline should be added here 

- Pipe4C 
- 

reefer to pipe4C and other repositoris that we have used for this project

## Citation


## Prerequisites

- A Unix like shell (we used Ubuntu 18.04)
- Bowtie2 v2.3+ available from http://bowtie-bio.sourceforge.net/bowtie2/.  (we used version 2.3.4.1 64-bit)
- SAMtools v1.3+ available from http://www.htslib.org/. 
   - **Note:** The pipeline will produce a sort error when older versions are used. (we used samtools 1.7)
- R v3.5+ available from https://www.r-project.org/.
- The following R packages available from CRAN:
  - optparse
  - caTools
  - config
- The following R packages available from Bioconductor:
  - shortRead
  - genomicRanges
  - genomicAlignments
  - BSgenome of interest
- The peakC package available from https://github.com/deWitLab/peakC/.

## Installation

Download the latest version of the pipeline from this git repository using:

```
    $ wget https://github.com/4CSeqEnhancerDetection/
    $ unzip master.zip
    $ cd ./pipe4C-master
```
**Note:** the pipe4C.R and functions.R files need to be placed in the same folder. 


## Steps

- We used index (e.g., CGATGT) in the reverse primers to seperate diffrent genes in one sequence run. Based on this index we seperated the genes. Since for one index, there is four lanes in the sequencing, we have to combine the fastq files for each index:

```
   cat *CGATGT*.fq.gz > all_CGATGT.fastq.gz
```

