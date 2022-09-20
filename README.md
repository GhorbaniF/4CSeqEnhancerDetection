# 4CSeq Putative Enhancers (PutE) Detection

To identify the Putative Enhancers (PutE) of four genes (ATXN1, ATXN3, TBP and ITPR1) in the human cerebellum, we used the following workflow: 
(I) Circularized chromosome conformation capture sequencing (4C-seq) on human cerebellum to capture genomic regions interacting with the promoters of the four  genes; 
(II) Peak calling with "PeakC" package to identify cis-genomic regions significantly interacting with the gene promoters (4C-contact peaks); 
(III) Annotating the 4C-contact peaks using publicly available data sets to identify Putative Enhancers (PutE); 

An example of the analyzed data for the TBP gene can be seen below.


![4CSeq Enhancer Detection](./imgs/tbp.png)


## 

This project has been built on top of the following packages, projects, and repositories:

- pipe4C: https://github.com/deLaatLab/pipe4C
- pyGenomeTracks: https://github.com/deeptools/pyGenomeTracks


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
- pyGenomeTracks package available from https://github.com/deeptools/pyGenomeTracks

## Installation Steps

### Step1: download the pipeline 

- First, Download all required files and functionalities:

```
    $ git clone https://github.com/GhorbaniF/4CSeqEnhancerDetection.git
    $ cd 4CSeqEnhancerDetection
    $ unzip pipe4C.zip
    $ cd pipe4C
```

**Note:** the pipe4C.R and functions.R files need to be placed in the same folder. 

### Step2: generating coverage plots and all relevant data to use for data analysis

- We used index (e.g., CGATGT) in the reverse primers to seperate diffrent genes in one sequence run. Based on this index we seperated the genes. Since for one index, there is four lanes in the sequencing, we have to combine the fastq files for each index:

```
cat *CGATGT*.fq.gz > all_CGATGT.fastq.gz
```

- Since our data is already demultiplexed, we have to run the following commands once just for generating /outF/ folder

```
Rscript pipe4C.R --vpFile=./example/VPinfo.txt --fqFolder=./example/ --outFolder=./outF/ --cores 8 --plot --wig
```

After couple of seconds, we need to stop the process using ctrl+c, and copy our fastq files, generated in the previous step, to the /outF/FASTQ/ folder. This way, the pipe4C program will use the existing fastq files and does not perform demultiplexing. 

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

- It is importnat to set the position of choromosom correctly in the VPinfo.txt file, otherwise peakC will not be able to visulize the peaks.

**Table 1.** example of VPinfo.txt file
| expname | spacer | primer | firstenzyme | secondenzyme	| genome	| vpchr	| vppos	| analysis | fastq |
| :---         | :---:  |     :---:  |     :---: |  :---: |  :---: |  :---: |  :---: |  :---: | :---: |
| all_GATCAG   | 0 | TCCAGACAAATAAACATG | NlaIII    | DpnII		   | hg19	| 14		| 92573009	| cis	     | all_GATCAG.fq.gz |
| all_ACTTGA   | 0 | TCCAGACAAATAAACATG | NlaIII    | DpnII		   | hg19	| 14		| 92573009	| cis	     | all_ACTTGA.fq.gz |
| all_CAGATC   | 0 | TCCAGACAAATAAACATG | NlaIII	 | DpnII		   | hg19	| 14		| 92573009	| cis		  | all_CAGATC.fq.gz |

- You can adjust 4C prameters in conf.yml file based on your experimental setup (e.g., genomes can be mm9, mm10, hg19, hg38)
- To set the Y axis of the coverage plots, which are generated in /outF/PLOTS/ folder, we need to edit the relevant parameters at the end of the conf.yml file

### Step3: peakC to call significant peaks:

- To generate the peakC plot, you need to add a set of addresses in **peakC_analysis.r** file, and run it in Rstudio. The parameters are including:
  - pipe4CFunctionsFile -> path to the **../pipe4C/functions.R** file
  - FileDirectory -> path to the **../pipe4C/outF/RDS** folder
  - If you need to change **alphaFDR, qWd, qWr** parameters for your experiments, they can be changed in **../pipe4C/functions.R** file (line ~1570, doPeakC function)
  - 

```
 doPeakC <- function(rdsFiles, vpRegion=2e6, wSize=21,alphaFDR=0.05,qWd=1.5,qWr=1,minDist=15e3) 
```

### Step4: identifing putative enhancers using public data with pygenoumtracks

#### Downloading Public Dataset: 
- We downloaded the **Atac-Seq**, **DNA-Seq** and **CHIP-Seq H3K27AC** bigWig files from **NCBI repository**:
   - Atac-Seq bigWig file (GSE138293_SKNSH.macs2.dnase.broad.SPMR.sorted.bw -- ATAC-seq SK-N-SH hg19 -- 506.13 Mb) from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138293
   
   - DNA-Seq bigWig file (GSM1008585_hg19_wgEncodeOpenChromDnaseSknshSig.bigWig -- 2.94GB) from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1008585

   - CHIP-Seq H3K27AC file (GSM2664335_SHSY5Y_2.K27ac.rep3.wig.bw -- 54.3 Mb) from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2664335
   
   
#### Install pyGenomeTracks and HICexplorer 
- We mainly use pyGenomeTracks and HICexplorer for visualizing and identifing putative enhancers using public data. Please check out the following web pages to obtain information about how to install pyGenomeTracks and HICexplorer and check their main functionalities:

	- PyGenomeTracks:
		- https://github.com/deeptools/pyGenomeTracks
		- https://pythonawesome.com/python-module-to-plot-beautiful-and-highly-customizable-genome-browser-tracks/
		- https://pygenometracks.readthedocs.io/en/latest/content/examples.html#examples-with-hi-c-data

	- HICexplorer (to visualize hic data)

		- https://hicexplorer.readthedocs.io/en/latest/index.html
		- https://hicexplorer.readthedocs.io/en/latest/content/installation.html


#### config *.ini file for visualizing all data using pyGenomeTracks: 

First, we need to create a TBP.ini file containing the following information:

```
[bigwig file test]
file = TBP/output.bw
color = red
height = 4
title = 4C-Seq 
min_value = 0
max_value = 1000

[spacer]

[narrow]
file = TBP/TBP_alphaFDR_0.05_qwr_1.bed
height = 1
type = box
color = black
#title = Peaks

[spacer]
height = 0.05

[bigwig file test]
file = ATAC_Seq/ATAC_Seq_GSE138293_SKNSH.macs2.dnase.broad.SPMR.sorted.bw
height = 2
color = #0000FF80
title = ATAC-Seq
min_value = 0
max_value = 0.5
#grid = 1

[spacer]
[bigwig file test]
file = DNA_Seq/ENCFF093ZBZ_encode.bw
height = 2
color = green
title = DNase-Seq
min_value = 0
max_value = 0.2

[spacer]
[bigwig file test]
file = CHIP_Seq/ENCFF712CGL_encode_H3K27ac.bw
height = 2
color = black
title = ChIP-Seq-H3K27ac
min_value = 0
max_value = 10

[spacer]
[bigwig file test]
file = CHIP_Seq/ENCFF545HWA_encode_H3K4ME1.bw
height = 2
color = black2
title = ChIP-Seq-H3K4me1
min_value = 0
max_value = 10


[spacer]

[vlines]
file = TBP/peacks_TBP.bed
type = vlines
line_width = 1

[vlines]
file = TBP/TBP3.bed
type = vlines
line_width = 1
[spacer]
#### peacks visualization 
[narrow]
file = TBP/intersect_TBP.bed
height = 1
type = box
color = brown
title = Intersect of peaks
[spacer]

[genes]
file = genome_genes_hg19/UCSC_exons_modif_canonical.bed
height = 5
#title = genes
fontsize = 10
file_type = bed
gene_rows = 10

[x-axis]
fontsize = 16
```

Then, open a terminal in the same path and run the following code: 

```
pyGenomeTracks --tracks TBP.ini --trackLabelFraction 0.2 --width 38 --dpi 130 --region chr6:169,000,000-171,000,000 -o TBP.png
```

