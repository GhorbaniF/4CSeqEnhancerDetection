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

  
| Name            | Description                                                                                                                                                                                   |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| fragFolder      | Path to the folder containing the fragment end libraries of the reference genomes                                                                                                             |
| normalizeFactor | Reads mapped to the 4C fragment end library are normalized to account for sequencing depth according to the normalizeFactor                                                                   |
| enzymes         | Enzyme names used in the viewpoint file and their corresponding recognition motifs                                                                                                            |
| genomes         | Genome names used in the viewpoint file plus corresponding BSgenome packages                                                                                                                  |
| bowtie2         | Path to corresponding bowtie2 index of reference genome. The reference genome assembly used to generate the index should match to the reference genome that was used to generate the BSgenome. NOTE: Some genome builds contain haplotype copies of the same chromosome. You should not include multiple haplotypes in the Bowtie2 index because the pipeline will consider reads mapping to two haplotypes as non-uniquely determined alignments. The pipeline removes chromosomes with 'hap' in the name when generating the fragmented genome.|
| maxY            | Maximal Y value in local 4C cis plot                                                                                                                                                          |
| plotView        | Number of bp to plot around viewpoint in local 4C cis plot                                                                                                                                    |
| xaxisUnit       | X-axis unit (Mb, Kb or bp)                                                                                                                                                                    |
| plotType        | Plots will either be saved as PDF or PNG                                                                                                                                                      |
| binSize         | Genome bin size used in the genome plot                                                                                                                                                       |
| qualityCutoff   | Q-score. Trim 3′-end of all sequences using a sliding window as soon as 2 out of 5 nucleotides have quality encoding less than the Q-score. 0 = no trimming                                   |
| trimLength      | Trim reads to defined capture length from 3′-end. 0 = no trimming                                                                                                                             |
| minAmountReads  | Minimum required amount of reads containing the primer sequence. If less reads are identified the experiment will not be further processed                                                    |
| readsQuality    | Bowtie2 minimum required mapping quality score for mapped reads                                                                                                                               |
| mapUnique       | Extract uniquely mapped reads, based on the lack of XS tag                                                                                                                                    |
| cores           | Number of CPU cores for parallelization                                                                                                                                                       |
| wSize           | The running mean window size                                                                                                                                                                  |
| nTop            | Top fragment ends discarded for calculation of normalizeFactor                                                                                                                                |
| nonBlind        | Only keep non-blind fragments                                                                                                                                                                 |
| wig             | Create wig files for all samples                                                                                                                                                              |
| plot            | Create viewpoint coverage plot for all samples                                                                                                                                                |
| genomePlot      | Create genomeplot for all samples (only possible if analysis is “all” in vpFile)                                                                                                              |
| tsv             | Create tab separated value file for all samples                                                                                                                                               |
| bins            | Count reads for binned regions                                                                                                                                                                |
| mismatchMax     | The maximum number of mismatches allowed during demultiplexing                                                                                                                                |
  
  **Table 1.** Description of parameters that need to be defined in the configuration file.
  
  <BR>

*** To set the Y axis of the coverage plots, which are generated in /outF/PLOTS/ folder, you need edit the relevant param at the end of the conf.yml file

*** To generate the pickC plot, you need to add a set of addresses in pickC_analysis.r file, and run it in Rstudio.



