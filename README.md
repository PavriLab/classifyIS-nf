# classifyIS-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

**classifyIS-nf** is a bioinformatics analysis pipeline for quantification of NS-seq data in mapped initiation sites (see [`iniseq-nf`](https://github.com/pavrilab/inisite-nf)) and their classification on the basis of this quantification.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Pipeline summary
The pipeline performs the following steps:

1.  Merging of called IS sites with [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html)
2.  Quantification of SNS-seq reads for either condition in merged IS with [`deeptools`](https://deeptools.readthedocs.io/en/develop/index.html)
3.  Perform classification of IS based on the given cutoff values
4.  Plot the results

## Quick Start
i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install [`BEDTools`](https://bedtools.readthedocs.io/en/latest/) and [`deepTools`](https://deeptools.readthedocs.io/en/develop/) and the [`pandas`](https://pandas.pydata.org/docs/index.html), [`numpy`](https://numpy.org/), [`scipy`](https://www.scipy.org/) and [`matplotlib`](https://matplotlib.org/) Python packages

iii. Clone repository with 
```bash
nextflow pull pavrilab/classifyIS-nf
```

iv. Start running your own analysis!
```bash
nextflow run pavrilab/classifyIS-nf --sitesA A_IS.bed --bamA conditionA.bam --labelA WT --sitesB B_IS.bed  --bamB conditionB.bam --labelB KD
````

## Main arguments
#### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. For example `-profile cbe` invokes the execution of processes using the [`slurm`](https://slurm.schedmd.com/documentation.html) workload manager. If no profile is given the pipeline will be executed locally.

#### `--sitesA`/`--sitesB`
Initiaton sites for conditions A and B (e.g. mapped with [`iniseq-nf`](https://github.com/pavrilab/inisite-nf))

#### `--bamA`/`--bamB`
BAM files containing aligned reads from which the initiation sites were mapped

## Generic arguments
#### `--labelA`/`--labelB`
Optional label for conditions A and B. Use this to customize labels of the conditions in the result files (defaults to A and B respectively)

#### `--FC`
log2(RPM) cutoff to use for assigning upregulation or downregulation (i.e. differential regulation) for a given peak (Default: 0.585 = log2(1.5))

#### `--t`
log2(RPM) cutoff to use for specifying putative peak call threshold (i.e. putatively dormant or absent peaks; Default: 2)

#### `--axMin`
Lower bound for axis values of x- and y-axis (Default: 0)

#### `--axMax`
Upper bound for axis values of x- and y-axis (Default: 8)

#### `--filePrefix`
Prefix for the result files name

#### `--outputDir`
Folder to which results will be written (is created if not existing)

## Results
The pipeline generates five result files:

1.  The `*.master.tsv` file is a tab-separated file with the basic layout of a BEDfile containing the genomic coordinates of the peaks their names and quantification results for either condition and the assigned class
2.  The `*.count.tsv` file contains the results as generated by deepTool multiBamSummary
3.  The two PDF files `*.density.pdf` and `*.class.pdf` are the visualization of the quantification results. The density plot show the distribution of peaks in the quantification space and the class plot visualizes the class assignment results

## Credits
The pipeline was developed by [Daniel Malzl](mailto:daniel.malzl@gmx.at) for use at the [IMP](https://www.imp.ac.at/), Vienna.

Many thanks to others who have helped out along the way too, including (but not limited to): [@t-neumann](https://github.com/t-neumann), [@pditommaso](https://github.com/pditommaso).

## Citations
### Pipeline tools
* [Nextflow](https://www.ncbi.nlm.nih.gov/pubmed/28398311/)
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

* [BEDTools](https://www.ncbi.nlm.nih.gov/pubmed/20110278/)
  > Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PubMed PMID: 20110278. PubMed Central PMCID: PMC2832824.
  
* [deepTools](https://www.ncbi.nlm.nih.gov/pubmed/27079975/)
  > Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dündar F, Manke T. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Res. 2016 Jul 8;44(W1):W160-5. doi: 10.1093/nar/gkw257. Epub 2016 Apr 13. PubMed PMID: 27079975; PubMed Central PMCID: PMC4987876.
  
### Python Packages
* [pandas](https://pandas.pydata.org/docs/index.html)
  > Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
  
* [numpy](https://numpy.org/)
  > Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011). doi: 10.1109/MCSE.2011.37

* [scipy](https://www.scipy.org/)
  > Virtanen, P., Gommers, R., Oliphant, T.E. et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17, 261–272 (2020). doi: 10.1038/s41592-019-0686-2
  
* [matplotlib](https://matplotlib.org/)
  > John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007). doi: 10.1109/MCSE.2007.55
