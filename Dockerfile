FROM ubuntu:16.04

MAINTAINER Basil Britto <basilbritto.xavier@uantwerpen.be> and Mohamed Mysara <mohamed.mysara@sckcen.be>
RUN apt update && apt-get -y install python3-pip curl
#RUN pip3 install --upgrade cutadapt

# Install Trim Galore
#RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.tar.gz -o trim_galore.tar.gz
#RUN tar xvzf trim_galore.tar.gz

#ENV PATH="/NGStools/TrimGalore-0.4.5/:$PATH"

# Install spades
RUN apt-get -y install bash python wget

WORKDIR /NGStools
#RUN wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz && tar -xf SPAdes-3.12.0-Linux.tar.gz
#ENV PATH="/NGStools/SPAdes-3.12.0-Linux/bin:$PATH"

# Install BLAST
# Blast+
#RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && tar xf ncbi-blast-2.7.1+-x64-linux.tar.gz && rm ncbi-blast-2.7.1+-x64-linux.tar.gz

#ENV PATH="/NGStools/ncbi-blast-2.7.1+/bin:${PATH}"

# Install Prokka
# Dependencies
RUN apt-get upgrade -y perl && apt-get install -y parallel make wget git python-pip locales && pip install -U setuptools && locale-gen --purge en_US.UTF-8 && DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8

# Barrnap
RUN wget https://github.com/tseemann/barrnap/archive/0.8.tar.gz && tar xf 0.8.tar.gz && rm 0.8.tar.gz

# Infernal
RUN wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz && tar xf infernal-1.1.2-linux-intel-gcc.tar.gz && rm infernal-1.1.2-linux-intel-gcc.tar.gz && cd infernal-1.1.2-linux-intel-gcc && ./configure && make && make install && cd easel && make install && cd /NGStools && rm -rf infernal-1.1.2-linux-intel-gcc/

# GenomeTools
RUN apt-get install -y libpango1.0-dev && wget http://genometools.org/pub/genometools-1.5.9.tar.gz && tar xf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make && make install && cd gtpython && python setup.py install && cd /NGStools && rm -rf genometools-1.5.9/

# gff3toembl
RUN mkdir gff3toembl && cd gff3toembl && git clone https://github.com/sanger-pathogens/gff3toembl.git && python gff3toembl/setup.py install && cd /NGStools

# CD-HIT
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-0621-source.tar.gz && tar xf cd-hit-v4.6.8-2017-0621-source.tar.gz && rm cd-hit-v4.6.8-2017-0621-source.tar.gz && cd cd-hit-v4.6.8-2017-0621 && make && cd cd-hit-auxtools && make && cd /NGStools

# Exonerate
RUN wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz && tar xf exonerate-2.2.0-x86_64.tar.gz --no-same-owner && chmod 755 exonerate-2.2.0-x86_64/ && rm exonerate-2.2.0-x86_64.tar.gz

# Prokka
RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre bioperl && git clone https://github.com/tseemann/prokka.git && prokka/bin/prokka --setupdb

# tbl2asn
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz && \
    gunzip linux.tbl2asn.gz && \
    mv linux.tbl2asn tbl2asn && \
    chmod +x tbl2asn && \
    mv tbl2asn /NGStools/prokka/binaries/linux/

ENV PATH="/NGStools/gff3toembl/gff3toembl/scripts:/NGStools/ncbi-blast-2.7.1+/bin:/NGStools/cd-hit-v4.6.8-2017-0621:/NGStools/cd-hit-v4.6.8-2017-0621/cd-hit-auxtools:/NGStools/exonerate-2.2.0-x86_64/bin:/NGStools/prokka/bin:${PATH}" LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH" PYTHONPATH="/NGStools/gff3toembl/gff3toembl:$PYTHONPATH"

# Install parsnp
WORKDIR /NGStools

#RUN wget https://github.com/marbl/parsnp/releases/download/v1.2/parsnp-Linux64-v1.2.tar.gz
#RUN tar -xvf parsnp-Linux64-v1.2.tar.gz

#ENV PATH="/NGStools/Parsnp-Linux64-v1.2:$PATH"

# Install CARD
#WORKDIR /NGStools/card

#RUN wget https://card.mcmaster.ca/download/6/prevalence-v3.0.1.tar.gz
#RUN tar -xvf prevalence-v3.0.1.tar.gz
#RUN rm prevalence-v3.0.1.tar.gz

#RUN mkdir nucleotide protein
#RUN mv protein_* protein/
#RUN mv nucleotide_* nucleotide/

# Install emmtyping
#WORKDIR /NGStools/emmtyping

#COPY emm_trimmed.fasta emm_trimmed.fasta

# Install mlstfinder
#WORKDIR /NGStools/mlstfinder

#RUN apt update && apt install unzip

#RUN wget https://bitbucket.org/genomicepidemiology/mlst_db/get/035f642871f6.zip
#RUN unzip 035f642871f6.zip

# Install plasmidfinder
#WORKDIR /NGStools/plasmidfinder

#COPY plasmidfinder_20_2_2017.zip plasmidfinder_20_2_2017.zip

#RUN unzip plasmidfinder_20_2_2017.zip

# Install resfinder
#WORKDIR /NGStools/resfinder

#COPY resfinder_2_11_2016.zip resfinder_2_11_2016.zip

#RUN unzip resfinder_2_11_2016.zip

# Install virdb
#WORKDIR /NGStools/virdb

#COPY VFDB_setA_pro.fas VFDB_setA_pro.fas

# Install virulencefinder
#WORKDIR /NGStools/virulencefinder

#COPY virulencefinder_16_3_2016.zip virulencefinder_16_3_2016.zip

#RUN unzip virulencefinder_16_3_2016.zip

# Install bacpipe
WORKDIR /NGStools/BacPipe

#RUN wget https://www.dropbox.com/s/am32gc7u49jigg1/BacPipe.v1.7.unix.run?dl=0
RUN curl -L https://www.dropbox.com/s/brsxvjlmbqjl2bv/BacPipe.zip?dl=0 > BacPipe.v1.7.unix.zip
#RUN mv BacPipe.v1.7.unix.run?dl=0 BacPipe.v1.7.unix.run
RUN unzip BacPipe.v1.7.unix.zip
#RUN chmod +x BacPipe.v1.7.unix.run

WORKDIR /NGStools

RUN apt install -y python-pip python-tk && pip install appJar

RUN \wget -O - https://install.perlbrew.pl | bash
ENV PATH=${PATH}:/root/perl5/perlbrew/bin/
RUN perlbrew install-cpanm
RUN cpanm Bio::Perl
RUN cpanm LWP::UserAgent
RUN cpanm Try::Tiny::Retry
RUN cpanm Excel::Writer::XLSX

WORKDIR /NGStools/BacPipe

#CMD /NGStools/BacPipe/BacPipe.v1.7.unix.run run
