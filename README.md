Putative Enhancers (PutE) Detection in human cerebellum tissue

To identify the Putative Enhancers (PutE) of four genes (ATXN1, ATXN3, TBP and ITPR1) in the human cerebellum, we used the following workflow: 
(I) Circularized chromosome conformation capture sequencing (4C-seq) on human cerebellum to capture genomic regions interacting with the promoters of the four  genes; 
(II) Peak calling with "PeakC" package to identify cis-genomic regions significantly interacting with the gene promoters (4C-contact peaks); 
(III) Annotating the 4C-contact peaks using publicly available data sets to identify Putative Enhancers (PutE); 

An example of the analyzed data for the TBP gene can be seen below.
![image](https://user-images.githubusercontent.com/25032978/191241838-03e62dfc-0479-4ea4-a4f2-ca11e53ea84e.png)



## 

This project has been built based on the two main following GitHub repositories:

- pipe4C: https://github.com/deLaatLab/pipe4C (Krijger, Peter HL, et al. "4C-seq from beginning to end: a detailed protocol for sample preparation and data analysis." Methods 170 (2020): 17-32)
- pyGenomeTracks: https://github.com/deeptools/pyGenomeTracks (Lopez-Delisle, Lucille, et al. "pyGenomeTracks: reproducible plots for multivariate genomic data sets." Bioinformatics (2021))


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

- Download all required files and functionalities:

```
    $ git clone https://github.com/GhorbaniF/4CSeqEnhancerDetection.git
    $ cd 4CSeqEnhancerDetection
    $ unzip pipe4C.zip
    $ cd pipe4C
```

**Note:** the pipe4C.R and functions.R files need to be placed in the same folder. 

### Step2: generating coverage plots and all relevant data to use for data analysis

-Since the four genes were pooled in one sequencing run (Illumina Nextseq500 using Mid-Output v2 kit with single-read run and 75bp read length), the four genes were first demultiplexed based on the forward primers to obtain the sequencing data for each gene (we received the data demultipexed from our service center, so it is not addressed here). 
-For each gene we have three replicates with three different indexes in the reverse primer which we use to seperate the three replicates of each gene. In this example we used index (e.g., CGATGT) in the reverse primer to seperate one of the replicates of the TBP gene from the other two. However, for each index there is four lanes which we have to combine the fastq files:

```
cat *CGATGT*.fq.gz > all_CGATGT.fastq.gz
```

-We do the same for the other two replicates (ATCACG and TTAGGC)
- Since our data was already demultiplexed and is not considered in our pipeline, we have to run the following commands once just for generating /outF/ folder

```
Rscript pipe4C.R --vpFile=./example/VPinfo.txt --fqFolder=./example/ --outFolder=./outF/ --cores 8 --plot --wig
```

After couple of seconds, we need to stop the process using ctrl+c, and copy our fastq files, generated in the previous step, to the /outF/FASTQ/ folder. This way, the pipe4C program will use the existing fastq files and does not perform demultiplexing (because it was already demultiplexed). 

- We finally run the script with our data using the following command:  

```
Rscript pipe4C.R --vpFile=./example/VPinfo.txt --fqFolder=./example/ --outFolder=./outF/ --cores 8 --plot --wig
```

**Notes:** 

- If everything works correctly, we should see the following lines in the terminal: 

```
------ Demultiplexing Fastq files based on VPinfo file
      ### WARNING: File ./outF/FASTQ/all_CGATGT.fastq.gz exists. continuing with exisiting file.
      ### WARNING: File ./outF/FASTQ/all_ATCACG.fastq.gz exists. continuing with exisiting file.
      ### WARNING: File ./outF/FASTQ/all_TTAGGC.fastq.gz exists. continuing with exisiting file.

```

- It is important to set the position of the choromosome correctly in the VPinfo.txt file, otherwise peakC will not be able to visualize the peaks.

**Table 1.** example of VPinfo.txt file
| expname | spacer | primer | firstenzyme | secondenzyme	| genome	| vpchr	| vppos	| analysis | fastq |
| :---         | :---:  |     :---:  |     :---: |  :---: |  :---: |  :---: |  :---: |  :---: | :---: |
| all_CGATGT   | 0 | ATCCGCTTCTTCTCCATG  | NlaIII    | DpnII		   | hg19	| 6		| 170863380	| cis	     | all_CGATGT.fq.gz |
| all_ATCACG   | 0 | ATCCGCTTCTTCTCCATG  | NlaIII    | DpnII		   | hg19	| 6		| 170863380	| cis	     | all_ATCACG.fq.gz |
| all_TTAGGC   | 0 | ATCCGCTTCTTCTCCATG  | NlaIII	 | DpnII		   | hg19	| 6		| 170863380	| cis		  | all_TTAGGC.fq.gz |

- You can adjust 4C parameters in conf.yml file based on your experimental setup (e.g., genomes can be mm9, mm10, hg19, hg38)
- To set the Y axis of the coverage plots, which are generated in /outF/PLOTS/ folder, we need to edit the relevant parameters at the end of the conf.yml file

With these steps we will generate the coverage plot and data for each of the replicates of TBP. In the next step, the data generated from the three replicates are used together to generate one plot for TBP. 


### Step3: peakC to call significant peaks:

- To generate the peakC plot, you need to add a set of addresses in **peakC_analysis.r** file, and run it in Rstudio. The parameters include:
  - pipe4CFunctionsFile -> path to the **../pipe4C/functions.R** file
  - FileDirectory -> path to the **../pipe4C/outF/RDS** folder
  - If you need to change **alphaFDR, qWd, qWr** parameters for your experiments, they can be changed in **../pipe4C/functions.R** file (line ~1570, doPeakC function)
  - 

```
 doPeakC <- function(rdsFiles, vpRegion=2e6, wSize=21,alphaFDR=0.05,qWd=1.5,qWr=1,minDist=15e3) 
```
In the case of TBP, the generated plot:

![image](https://user-images.githubusercontent.com/25032978/191950645-ce737028-c55a-41f2-91d5-841061a71f79.png)


### Step4: identifying putative enhancers using public data with pygenometracks

#### Downloading Public Dataset: 
- We downloaded the **ATAC-Seq**, **DNA-Seq**, **CHiP-Seq H3K27AC** and **HiC-Seq** bigWig files for human cerebellar tissue:
ATAC-seq (GSE101912)
ChiP-seq for H3K27ac (GSM1119154)
DNase-seq (GSM736538) 
HiC-seq (ENCFF198DDF) 
   
#### Install pyGenomeTracks and HICexplorer 
- We mainly use pyGenomeTracks and HICexplorer for visualizing and identifying putative enhancers using public data. Please check out the following web pages to obtain information about how to install pyGenomeTracks and HICexplorer and check their main functionalities:

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
file = ATAC_Seq/GSE101912.bw
height = 2
color = #0000FF80
title = ATAC-Seq
min_value = 0
max_value = 0.5
#grid = 1

[spacer]
[bigwig file test]
file = DNA_Seq/GSM736538.bw
height = 2
color = green
title = DNase-Seq
min_value = 0
max_value = 0.2

[spacer]
[bigwig file test]
file = CHIP_Seq/GSM1119154.bw
height = 2
color = black
title = ChIP-Seq-H3K27ac
min_value = 0
max_value = 10


[spacer]

[vlines]
file = TBP/peaks_TBP.bed
type = vlines
line_width = 1

[vlines]
file = TBP/TBP3.bed
type = vlines
line_width = 1
[spacer]
#### peaks visualization 
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

The final plot for TBP will look like this:
![TBP](https://user-images.githubusercontent.com/25032978/191953635-d61c7b5d-7c22-4ba5-97e4-79a506eb661b.png)

