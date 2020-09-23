FROM ubuntu:16.04

MAINTAINER Basil Britto <basilbritto.xavier@uantwerpen.be> and Mohamed Mysara <mohamed.mysara@sckcen.be>

RUN apt update && apt-get -y install python3-pip curl

RUN apt-get update && apt-get -y install \
bash \
python \
python3 \
wget

#RUN echo Y | apt-get install software-properties-common
#RUN add-apt-repository ppa:deadsnakes/ppa
#RUN apt-get update
#RUN echo Y | apt-get install python3.6

RUN apt-get upgrade -y perl && apt-get install -y parallel make wget git python-pip locales && pip install --upgrade pip && pip install -U setuptools && locale-gen --purge en_US.UTF-8 && DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8

RUN git clone https://github.com/wholeGenomeSequencingAnalysisPipeline/BacPipe.git

WORKDIR /BacPipe

#RUN apt-get unzip

RUN chmod 755 /BacPipe/SPAdes-3.13.0-Linux/bin/*

RUN gunzip -c /BacPipe/mlst/blast-2.2.26/bin/makeblastdb.zip > /BacPipe/mlst/blast-2.2.26/bin/makeblastdb

#RUN unzip /BacPipe/mlst/blast-2.2.26/bin/makeblastdb.zip -d /BacPipe/mlst/blast-2.2.26/bin/

RUN chmod +x -R /BacPipe/

RUN apt install -y python-pip python-tk

RUN pip3 install --upgrade cutadapt

# Dependencies

RUN apt-get update && apt-get upgrade -y perl && apt-get install -y parallel make wget git python-pip locales && pip install -U setuptools && locale-gen --purge en_US.UTF-8 && DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8


##prokka#################
RUN apt-get update && apt-get -y install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre bioperl
RUN cpan Test::Most 
RUN cpan Bio::Perl && cpan List::Util 
RUN git clone https://github.com/tseemann/prokka.git --branch v1.13.3
#RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz && tar -zxvf ncbi-blast-2.9.0+-x64-linux.tar.gz
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && tar -zxvf ncbi-blast-2.9.0+-x64-linux.tar.gz
ENV PATH=${PATH}:"ncbi-blast-2.9.0+/bin"
RUN cp ./ncbi-blast-2.9.0+/bin/* /usr/bin/
RUN prokka/bin/prokka --setupdb
WORKDIR /BacPipe
######################

# Barrnap
RUN wget https://github.com/tseemann/barrnap/archive/0.8.tar.gz && tar xf 0.8.tar.gz && rm 0.8.tar.gz


# Infernal
RUN wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz && tar xf infernal-1.1.2-linux-intel-gcc.tar.gz && rm infernal-1.1.2-linux-intel-gcc.tar.gz && cd infernal-1.1.2-linux-intel-gcc && ./configure && make && make install && cd easel && make install && rm -rf infernal-1.1.2-linux-intel-gcc/

# GenomeTools
RUN apt-get install -y libpango1.0-dev && wget http://genometools.org/pub/genometools-1.5.9.tar.gz && tar xf genometools-1.5.9.tar.gz && rm genometools-1.5.9.tar.gz && cd genometools-1.5.9 && make && make install && cd gtpython && python setup.py install && rm -rf genometools-1.5.9/ && cd /BacPipe/
# gff3toembl
RUN mkdir gff3toembl && cd gff3toembl && git clone https://github.com/sanger-pathogens/gff3toembl.git && python gff3toembl/setup.py install && cd /BacPipe/
# CD-HIT
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-0621-source.tar.gz && tar xf cd-hit-v4.6.8-2017-0621-source.tar.gz && rm cd-hit-v4.6.8-2017-0621-source.tar.gz && cd cd-hit-v4.6.8-2017-0621 && make && cd cd-hit-auxtools && make && cd /BacPipe/
# Exonerate
RUN wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz && tar xf exonerate-2.2.0-x86_64.tar.gz --no-same-owner && chmod 755 exonerate-2.2.0-x86_64/ && rm exonerate-2.2.0-x86_64.tar.gz

RUN wget -O - https://install.perlbrew.pl | bash
ENV PATH=${PATH}:/root/perl5/perlbrew/bin/
RUN perlbrew install-cpanm
RUN cpanm Bio::Perl
RUN cpanm -f LWP::UserAgent
RUN cpanm Try::Tiny::Retry
RUN cpanm Excel::Writer::XLSX





#install plasmidfinder
RUN cd /BacPipe; \
apt-get update -qq; \
apt-get install -y -qq git \
apt-utils \
wget \
python3-pip \
ncbi-blast+ \
libz-dev \
; \
rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

# Install python dependencies
RUN pip3 install --upgrade pip
RUN pip3 install -U biopython==1.76 tabulate cgecore==1.3.6;

# Install kma
RUN git clone --branch 1.0.1 --depth 1 https://bitbucket.org/genomicepidemiology/kma.git; \
cd kma && make; \
mv kma* /bin/
WORKDIR /BacPipe

#PLASMIDFINDER update
##rename outdated plasmidfinder
RUN mv /BacPipe/plasmidfinder /BacPipe/plasmidfinder_perl
##download the new version
RUN git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git
WORKDIR /BacPipe/plasmidfinder
RUN chmod 755 ./plasmidfinder.py;
##Download the database
RUN git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
WORKDIR /BacPipe/plasmidfinder/plasmidfinder_db
RUN PLASMID_DB=$(pwd)
RUN python3 /BacPipe/plasmidfinder/plasmidfinder_db/INSTALL.py kma_index
WORKDIR /BacPipe

#RESFINDER updates (no need to upgrade the tool only the database)
#RUN mv /BacPipe/resfinder /BacPipe/resfinder_perl
#RUN git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
WORKDIR /BacPipe/resfinder
#RUN chmod 755 ./resfinder.py;
##Download the database
RUN git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
WORKDIR /BacPipe/resfinder/resfinder_db
RUN RES_DB=$(pwd)
RUN python3 /BacPipe/resfinder/resfinder_db/INSTALL.py kma_index
WORKDIR /BacPipe

#MLST updates
##rename outdated plasmidfinder
#RUN mv /BacPipe/mlst /BacPipe/mlst_perl
##Download the new version
#RUN git clone https://bitbucket.org/genomicepidemiology/mlst.git
#WORKDIR /BacPipe/mlst
#RUN chmod 755 ./mlst.py
##Download the database
#RUN git clone https://bitbucket.org/genomicepidemiology/mlst_db.git
#WORKDIR /BacPipe/mlst/mlst_db
#RUN MLST_DB=$(pwd)
#RUN python3 /BacPipe/mlst/mlst_db/INSTALL.py kma_index
#RUN mv /BacPipe/mlst_perl/blast-2.2.26 /BacPipe/mlst/
WORKDIR /BacPipe

#prokka

##############
#prokka
#WORKDIR /BacPipe

#RUN apt-get install libdatetime-perl
#RUN apt-get upgrade -y perl && apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl



#RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre bioperl && git clone https://github.com/tseemann/prokka.git && prokka/bin/prokka --setupdb

#RUN apt-get install -y libdatetime-perl
#RUN apt-get install -y libxml-simple-perl
#RUN apt-get install -y libdigest-md5-perl
#RUN apt-get install -y default-jre bioperl
#RUN git clone https://github.com/tseemann/prokka.git  && prokka/bin/prokka --setupdb
#####################


# tbl2asn
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz && \
gunzip linux.tbl2asn.gz && \
mv linux.tbl2asn tbl2asn && \
chmod +x tbl2asn && \
mv tbl2asn /BacPipe/prokka/binaries/linux/

ENV PATH=${PATH}:"/BacPipe/gff3toembl/gff3toembl/scripts:/BacPipe/ncbi-blast-2.7.1+/bin:/BacPipe/cd-hit-v4.6.8-2017-0621:/BacPipe/cd-hit-v4.6.8-2017-0621/cd-hit-auxtools:/BacPipe/exonerate-2.2.0-x86_64/bin:/BacPipe/prokka/bin:${PATH}" LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH" PYTHONPATH="/BacPipe/gff3toembl/gff3toembl:$PYTHONPATH"


