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

- Since our data is already demultiplexed, we have to run the following commands once just for generating /outF/ folder

```
Rscript pipe4C.R --vpFile=./example/VPinfo.txt --fqFolder=./example/ --outFolder=./outF/ --cores 8 --plot --wig
```

After couple of seconds, we need to stop the process using ctrl+c, and copy our fastq files generated in the previous step to the /outF/FASTQ/ folder. This way, the pipe4C program will use the existing fastq files and does not perform demultiplexing. 

- We finally run the script with our data using the following command:  

```
Rscript pipe4C.R --vpFile=./example/VPinfo.txt --fqFolder=./example/ --outFolder=./outF/ --cores 8 --plot --wig
```

**Notes:** 

- If everything works correctly, we should see the following lines in the terminal: 

```
------ Demultiplexing Fastq files based on VPinfo file
      ### WARNING: File ./outF/FASTQ/all_ACTTGA.fastq.gz exists. continuing with exisiting file.
      ### WARNING: File ./outF/FASTQ/all_CAGATC.fastq.gz exists. continuing with exisiting file.
      ### WARNING: File ./outF/FASTQ/all_GATCAG.fastq.gz exists. continuing with exisiting file.

```

- It is importnat to set the position of choromosom correctly in the VPinfo.txt file, otherwise pickC can not visulize the picks. 


|expname     | spacer| primer	            | firstenzyme 	| secondenzyme	| genome	| vpchr	| vppos		| analysis | fastq 	        | 
|all_GATCAG  | 0		| TCCAGACAAATAAACATG	| NlaIII		   | DpnII		   | hg19	| 14		| 92573009	| cis	     | all_GATCAG.fq.gz|
|all_ACTTGA  | 0		| TCCAGACAAATAAACATG	| NlaIII		   | DpnII		   | hg19	| 14		| 92573009	| cis	     | all_ACTTGA.fq.gz|
|all_CAGATC  | 0		| TCCAGACAAATAAACATG	| NlaIII		   | DpnII		   | hg19	| 14		| 92573009	| cis		  | all_CAGATC.fq.gz|

**Table 1.** example of VPinfo.txt file


*** To set the Y axis of the coverage plots, which are generated in /outF/PLOTS/ folder, you need edit the relevant param at the end of the conf.yml file

*** To generate the pickC plot, you need to add a set of addresses in pickC_analysis.r file, and run it in Rstudio.



