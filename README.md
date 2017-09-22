# BacPipe

## A rapid, user-friendly whole genome sequencing pipeline for clinical diagnostic bacteriology and outbreak detection


Despite rapid advances in whole genome sequencing (WGS) technologies, their integration into routine microbiological diagnostics and infection control has been hampered by the need for downstream bioinformatics analyses that require considerable expertise. We have developed a comprehensive, rapid, and computationally low-resource bioinformatics pipeline (BacPipe) for the analysis of bacterial whole genome sequences obtained from second and third-generation sequencing technologies. Users can choose to directly analyse raw sequencing reads or contigs or scaffolds in BacPipe. The pipeline is an ensemble of state-of-the-art, open-access bioinformatics tools for quality verification, genome assembly and Annotation, and identification of the bacterial genotype (MLST and emm typing), resistance genes, plasmid(s), virulence genes, and single nucleotide polymorphisms (SNPs). The outbreak module in BacPipe can be used, along with the SNPs and patient metadata, to simultaneously analyse many strains to understand evolutionary relationships and rapidly construct bacterial transmission routes. Importantly, BacPipe is designed to run multiple tools simultaneously which considerably reduces the time-to-result. We validated BacPipe using prior published WGS datasets from confirmed outbreaks of MRSA, carbapenem-resistant Klebsiella pneumoniae, and Salmonella enterica, and from transmission studies of Clostridium difficile and Mycobacterium tuberculosis where BacPipe helped build the same analyses and conclusions within a few hours. We believe this fully automated pipeline will contribute to overcoming one of the primary hurdles faced by microbiologists for analysing and interpreting WGS data, facilitating its direct application for routine patient care in hospitals and public health and infection control monitoring.

BacPipe software can be downloaded from the release section here (https://github.com/wholeGenomeSequencingAnalysisPipeline/BacPipe/releases) 


# Installation 
## Automatic installation
Perl, java and Python (2.7) needed to be installed to run BacPipe, normally pre-installed in most unix/linux and macOS (except java for macOS). If not, Python and Perl can be downloaded and installed via these instructions: https://www.python.org/downloads/, https://www.java.com/en/downloads and https://www.perl.org/get.html respectively.


An automatic installation script is added (but require sudo privilege). To start the installation first you have to give execution permission to the ".run" file:
```
chmod +x BacPipe.v?.?.*.run
```
Then you need to run the installation script, as following:
```
./BacPipe.v?.?.mac.run install /PATH/FOR/Prokka
./BacPipe.v?.?.unix.run install PATH/FOR/Prokka
```
or to the extract the .run file into a folder and install there, please add "--target PATH" after running the tool, as following:
```
./BacPipe.v?.?.run install /PATH/FOR/Prokka --target /PATH/to/extraction/location/
cd /PATH/to/extraction/location/
./Main.command install PATH/FOR/Prokka
```
#### Please note that Prokka PATH has to be absolute PATH i.e. "/home/users/yourname/Output" rather than "./Output".
#### The installation require sudo privilege if you do not have it, you may use the instructions below without "sudo".

## Detailed installation steps (if needed)

The following Perl packages needed to be installed to run BacPipe: Time::Piece, XML::Simple, Bio::Perl, Digest::MD5, Try::Tiny::Retry, Bio::TreeIO, SVG::Graph and Excel::Writer::XLSX. To install them follow these commands: 

You need to install them as following:
```bash
sudo cpan Try::Tiny::Retry
sudo cpan Excel::Writer::XLSX
```
**Centos/Fedora/RHEL (RPM)**
```bash
sudo yum install perl-Time-Piece
sudo yum install perl-XML-Simple
sudo yum install perl-Digest-MD5
sudo yum install git java perl-CPAN perl-Module-Build
sudo cpan -i Bio::Perl
```
**Ubuntu/Debian/Mint (APT)**
```bash
sudo apt-get install libdatetime-perl
sudo apt-get install libxml-simple-perl
sudo apt-get install libdigest-md5-perl
sudo apt-get install git default-jre bioperl
```
**Mac OS X**
```bash
sudo cpan Time::Piece
sudo cpan XML::Simple
sudo cpan Digest::MD5
sudo cpan Bio::Perl
```


The following Python packages needed to be installed to run BacPipe: appjar, Yaml. To install them follow these commands: 
	
	sudo easy_install pip
	sudo pip install appjar (or follow the instructions here: appjar.info/Install/)
	sudo pip install pyyaml, python -m pip install, sudo apt-get install python-yaml, or sudo yum install python-yaml

Additionally, other software are required to run BacPipe such as Prokka (annotation tool) and tbl2asn. Detailed Prokka installation can be found here: 

	http://www.vicbioinformatics.com/software.prokka.shtml

tbl2asn can be obtained and installed as following (mac):

	curl -O ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/mac.tbl2asn.gz
	gunzip mac.tbl2asn.gz
	mv mac.tbl2asn tbl2asn
	chmod +x tbl2asn
	sudo cp tbl2asn /usr/local/bin

or for unix/Linux tbl2asn can be obtained and installed as following:

	curl -O ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
	gunzip linux64.tbl2asn.gz
	mv linux64.tbl2asn tbl2asn
	chmod +x tbl2asn
	sudo cp tbl2asn /usr/local/bin



## For Windows users (or those who prefer Virtualization)

BacPipe was tested on macOS Sierra and Linux (Mint 18 "Sarah") both virtual and physical, thus I recommend using either one of them as your virtual machine.
A very good demonistrations on how to create a virtual machine with macOS Sierra or Mint OS are shown via this video:
```
https://www.youtube.com/watch?v=_cPlKqp8nEA
https://www.youtube.com/watch?v=_95fA2hmGYM
```
To install MacOS Sierra on Windows, you can follow the instructions detailed via this link (please check its legality within MacOS license):
```
http://www.wikigain.com/install-macos-sierra-10-12-virtualbox/
```
Please not that in Step #5 "Add VirtualBox Code to the CMD" the CMD has to be opened with administrator privilege and the Virtual tool has to be closed.
	
# Included software and licence
Software/Database listed below is used by the BacPipe pipeline. However you do NOT need to install them separately as these software modules are included in the BacPipe software.

VFDB:

	Yang J, Chen LH, Sun LL, Yu J and Jin Q, 2008. VFDB 2008 release: an enhanced web-based resource for comparative pathogenomics. Nucleic Acids Res.
	Available at: http://www.mgc.ac.cn/VFs/download.htm
Emm Typing:

	Facklam R, Beall B, Efstratiou A, Fischetti V, Johnson D, Kaplan E, et al. emm Typing and Validation of Provisional M Types for Group A Streptococci. Emerg Infect Dis. 1999;5(2):247-253. 
	Available at: https://www2a.cdc.gov/ncidod/biotech/strepblast.asp
	Licence: GNU
ParSNP (v1.2):

	Treangen TJ, Ondov BD, Koren S, Phillippy AM (2014) Rapid Core-Genome Alignment and Visualization for Thousands of Microbial Genomes.
	Available at: http://harvest.readthedocs.io/en/latest/content/parsnp.html
	Licence owner: Battelle National Biodefense Institute (BNBI)
CARD (v1.1.7):

	Jia et al. 2017. CARD 2017: expansion and model-centric curation of the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 45, D566-573.
	Available at: https://card.mcmaster.ca
	Licence owner: McMaster University
VirulenceFinder (v1.4):

	Joensen KG, Scheutz F, Lund O, Hasman H, Kaas RS, Nielsen EM, Aarestrup FM, Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli. J. Clin. Micobiol. 2014. 52(5): 1501-1510.
	Available at: http://www.genomicepidemiology.org
	Licence: Apache License
ResFinder (v2.1):

	Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV, Identification of acquired antimicrobial resistance genes, J Antimicrob Chemother. 2012 Jul 10.
	Available at: http://www.genomicepidemiology.org
	Licence: Apache License 
PlasmidFinder (v1.3):

	Carattoli A, Zankari E, Garcia-Fernandez A, Voldby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.,PlasmidFinder and pMLST: in silico detection and typing of plasmids. Antimicrob. Agents Chemother. 2014. April 28th.
	Available at: http://www.genomicepidemiology.org
	Licence: Apache License 
MLST (v1.8):

	Larsen MV, Cosentino S, Rasmussen S, Friis C, Hasman H, Marvig RL, Jelsbak L, Sicheritz-Pontén T, Ussery DW, Aarestrup FM and Lund O., Multilocus Sequence Typing of Total Genome Sequenced Bacteria. J. Clin. Micobiol. 2012. 50(4): 1355-1361.
	Available at: http://www.genomicepidemiology.org
	Licence: Apache License 
Resfams (v1.2):
	
	Gibson MK, Forsberg KJ, Dantas G. Improved annotation of antibiotic resistance functions reveals microbial resistomes cluster by ecology. The ISME Journal. 2014, doi:ISMEJ.2014.106')
	Available at: http://www.dantaslab.org/resfams/
	Licence: 
TrimGalore (v0.4.2):

	Available at: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
	Licence: GNU
SPAdes (v3.10.0):

	Bankevich, A., Nurk, S., Antipov, D., Gurevich, A.A., Dvorkin, M., et al.: SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology 19(5), 455–477 (2012)
	Available at http://spades.bioinf.spbau.ru/
	Licence: GNU
NCBI Blast (v2.4.17):

	Johnson M, Zaretskaya I, Raytselis Y, Merezhuk Y, McGinnis S, & Madden T.L. (2008) "NCBI BLAST: a better web interface" Nucleic Acids Res. 36:W5-W9.
	Available at https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/
	Licence: GNU

cutadapt (v1.12):

	MARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may. 2011.
	Available at https://pypi.python.org/pypi/cutadapt 
	Licence: MIT license




	
         
Few of these software (in addition to BacPipe) are under the GNU/GPL licence. Others under Apache license (ResFinder, 

In addition due to license restrictions, usearch has to be downloaded and installed by the user. It will be parsed to the OCToPUS via _u option (see below)

    usearch 8.1.1861
         Available at http://www.drive5.com/usearch/
	 
# Running the pipeline
The pipeline GUI is initiated by double clicking the executable (.run). Alternatively, the tool can be called from the terminal via running:
```
./BacPipe.v?.?.run run
```
If the user wants the single bundel to be extracted in a specific location (rather than the tmp folder which is deleted post analysis), please add "--target PATH" after running the tool, as following:
```
./BacPipe.v?.?.run run --target /PATH/to/extraction/location/
```

The user can start from the raw sequencing data full a comprehensive analysis or a assign specific tools to run. Alternatively, it is possible to start at a specific step within the pipeline and assign specific tools to run. 
## Inputs:
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
## Progress:

The progress is shown in the Progress tab, where for each sample the percentage of completed steps are shown. When the run is finished the results tab will be accessible and the results will be shown.
	
## Output:

There are three types of outputs, overall, summarised and details. All will be produced at the end of the analysis.

### A) Overall: 

A cross sample illustration of the various tools results. It is useful where more than one sample are analysed. It only covers MLST, resFinder, virulenceFinder and PlasmidFinder tools results. In cases where the identified gene is found in more than one location within a sample, their % identify to that gene are shown within the same cell (Comma separated). Additionally, a text format of sap-diversity tree will be illustrated (when selecting the outbreak option).
	
### B) Summarised:

An excel excel file for each sample will be created. Within this sheet, the tools results will be shown (one sheet per tool). These excel files will be grouped in the â€œSummaryâ€ folder.
	
### C) Detailed:

A folder will be created per sample within the output directory. Within this folder a sub-folders for each of the tools output will be created stating the detailed results for each tools.
	
Additionally, a log-file will be created stating the step-wise update of the tools performed (pipeline.log) and a detailed log-file with their output (log.txt) for any possible errors.


# Contact: 

For questions, bugs and suggestions, please refer to XX@X.X
Developed by XXX 2017
