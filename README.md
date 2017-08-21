# BacPipe

## A rapid, user-friendly whole genome sequencing pipeline for clinical diagnostic bacteriology and outbreak detection


Despite rapid advances in whole genome sequencing (WGS) technologies, their integration into routine microbiological diagnostics and infection control has been hampered by the need for downstream bioinformatics analyses that require considerable expertise. We have developed a comprehensive, rapid, and computationally low-resource bioinformatics pipeline (BacPipe) for the analysis of bacterial whole genome sequences obtained from second and third-generation sequencing technologies. Users can choose to directly analyse raw sequencing reads or contigs or scaffolds in BacPipe. The pipeline is an ensemble of state-of-the-art, open-access bioinformatics tools for quality verification, genome assembly and Annotation, and identification of the bacterial genotype (MLST and emm typing), resistance genes, plasmid(s), virulence genes, and single nucleotide polymorphisms (SNPs). The outbreak module in BacPipe can be used, along with the SNPs and patient metadata, to simultaneously analyse many strains to understand evolutionary relationships and rapidly construct bacterial transmission routes. Importantly, BacPipe is designed to run multiple tools simultaneously which considerably reduces the time-to-result. We validated BacPipe using prior published WGS datasets from confirmed outbreaks of MRSA, carbapenem-resistant Klebsiella pneumoniae, and Salmonella enterica, and from transmission studies of Clostridium difficile and Mycobacterium tuberculosis where BacPipe helped build the same analyses and conclusions within a few hours. We believe this fully automated pipeline will contribute to overcoming one of the primary hurdles faced by microbiologists for analysing and interpreting WGS data, facilitating its direct application for routine patient care in hospitals and public health and infection control monitoring.

OCToPUS software can be downloaded from the release section here (https://github.com/M-Mysara/OCToPUS/releases) 

# Installation Requirement
Both Perl and Java needed to be installed to run OCToPUS. All other software packages that are required to run OCToPUS are included in the downloaded file (OCToPUS_V?.run). In case you are interested in the source code of OCToPUS, this is also included in the downloaded file. Only in case you want to run the source code, you will need to install those software components separately, and adapt the source code referring to those software components accordingly. In all other cases, we encourage the end-user to use the OCToPUS_V?.run executable.

# Included Software
Software listed below is used by the OCToPUS algorithm. However you do NOT need to install it separately as these software modules are included in the OCToPUS software.

    Mothur v.1.33.3:
         Available at http://www.mothur.org/wiki/Download_mothur. 
         Note about changes made in the mothur package integrated in this package:
         The command called "pre.cluster" is modified to be compatible with the IPED algorithm.
         The command called "make.contigs" is modified to produce an additional IPED-formatted quality file.
    WEKA 3.7.11: 
         Available online at http://www.cs.waikato.ac.nz/ml/weka/.
    SPAdes 3.5.0:   
         Available at http://spades.bioinf.spbau.ru/
    IPED v.1
         Available online at https://github.com/M-Mysara/IPED/
         Note about changes made in IPED to remove redundant steps and increase the compatabilities   
    CATCh v.1
         Available online at  http://science.sckcen.be/en/Institutes/EHS/MCB/MIC/Bioinformatics/CATCh
         Note about changes made in CATCh to remove redundant steps and increase the compatabilities
         
All of these software (in addition to OCToPUS) are under the GNU licence.

In addition due to license restrictions, usearch has to be downloaded and installed by the user. It will be parsed to the OCToPUS via _u option (see below)

    usearch 8.1.1861
         Available at http://www.drive5.com/usearch/

# Syntax:
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 ||                      Welcome To LMM WGS pipeline                ||
 || Software for analysing whole genome sequencing data||
 ||     for clinical diagnostics and outbreaks assessment     ||
 ||                                 Copyright (C) 2017                            ||
 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 
 This is free software, under GNU, also it includes other software also under GNU Copyright
 
 ####################################################################################################
 Settings and inputs:
 The pipeline can be run in whole or from intermediate/specific step(s), thus it is important to know the expected input for each step if you want to run specific step(s). It start with the raw reads (forward and reverse *.fastq.gz) and process them as following:
 a) Trim Galore: quality trimming
	INPUT: *_R1*.fastq.gz/*_R2*.fastq.gz
	OUTPUT: *_R1_001_val_1.fq.gz/*_R2_001_val_2.fq.gz
b) SPAdes: assembly and scaffolding
	INPUT: *_R1_001_val_1.fq.gz/*_R2_001_val_2.fq.gz
	OUTPUT: *.fasta
c) MLST typing:
	INPUT: *.fasta (DNA sequences)
d) emm typing (for Streptococcus):
	INPUT: *.fasta (DNA sequences)
e) Plasmids Finder:
	INPUT: *.fasta (DNA sequences)
f) ResFinder (Plasmid mediated Resistance only):
	INPUT: *.fasta (DNA sequences)
g) Virulence Finder:
	INPUT: *.fasta (DNA sequences)
h) Prokka (annotation):
	INPUT: *.fasta (DNA sequences)
i) CARD search (extensive resistance search):
	INPUT: *.faa (poteins sequences)
j) VirDB search:
	INPUT: *.faa (poteins sequences)
k) Outbreak assessment (parSNP)
	INPUT: *.fasta (DNA sequences, > 2 files)
l) ResFams search (extensive resistance search):
	INPUT: *.faa (poteins sequences)
m) Summarise output:
	OUTPUT: excel file for each sample (each selected tools shown in one sheet)
	
	The progress is shown in the Progress tab, where for each sample the percentage of completed steps are shown. When the run is finished the results tab will be accessible and the results will be shown.
	
####################################################################################################
Output:
There are three types of outputs, overall, summarised and details. All will be produced at the end of the analysis.
	A) Overall: A cross sample illustration of the various tools results. It is useful where more than one sample are analysed. It only covers MLST, resFinder, virulenceFinder and PlasmidFinder tools results. In cases where the identified gene is found in more than one location within a sample, their % identify to that gene are shown within the same cell (Comma separated). Additionally, a text format of sap-diversity tree will be illustrated (when selecting the outbreak option).
	
	B) Summarised: An excel excel file for each sample will be created. Within this sheet, the tools results will be shown (one sheet per tool). These excel files will be grouped in the â€œSummaryâ€ folder.
	
	C) Detailed:
	A folder will be created per sample within the output directory. Within this folder a sub-folders for each of the tools output will be created stating the detailed results for each tools.
	
	Additionally, a log-file will be created stating the step-wise update of the tools performed (pipeline.log) and a detailed log-file with their output (log.txt) for any possible errors.
