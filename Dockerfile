FROM ubuntu:16.04

MAINTAINER Basil Britto <basilbritto.xavier@uantwerpen.be> and Mohamed Mysara <mohamed.mysara@sckcen.be>

RUN apt update && apt-get -y install python3-pip curl

RUN apt-get update && apt-get -y install \
	bash \
	python \
	python3 \
	wget

RUN apt-get upgrade -y perl && apt-get install -y parallel make wget git python-pip locales && pip install -U setuptools && locale-gen --purge en_US.UTF-8 && DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8

RUN git clone https://github.com/wholeGenomeSequencingAnalysisPipeline/BacPipe.git

WORKDIR /BacPipe

#RUN cd BacPipe

RUN chmod +x -R /BacPipe/

#RUN gunzip /BacPipe/mlst/blast-2.2.26/bin/makeblastdb.zip

RUN apt install -y python-pip python-tk 

RUN pip3 install --upgrade cutadapt

# Dependencies

RUN apt-get update && apt-get upgrade -y perl && apt-get install -y parallel make wget git python-pip locales && pip install -U setuptools && locale-gen --purge en_US.UTF-8 && DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8



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
#prokka
RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre bioperl && git clone https://github.com/tseemann/prokka.git && prokka/bin/prokka --setupdb
# tbl2asn
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz && \
    gunzip linux.tbl2asn.gz && \
    mv linux.tbl2asn tbl2asn && \
    chmod +x tbl2asn && \
    mv tbl2asn /BacPipe/prokka/binaries/linux/

ENV PATH="/BacPipe/gff3toembl/gff3toembl/scripts:/BacPipe/ncbi-blast-2.7.1+/bin:/BacPipe/cd-hit-v4.6.8-2017-0621:/BacPipe/cd-hit-v4.6.8-2017-0621/cd-hit-auxtools:/BacPipe/exonerate-2.2.0-x86_64/bin:/BacPipe/prokka/bin:${PATH}" LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH" PYTHONPATH="/BacPipe/gff3toembl/gff3toembl:$PYTHONPATH"

RUN wget -O - https://install.perlbrew.pl | bash
ENV PATH=${PATH}:/root/perl5/perlbrew/bin/
RUN perlbrew install-cpanm
RUN cpanm Bio::Perl
RUN cpanm LWP::UserAgent
RUN cpanm Try::Tiny::Retry
RUN cpanm Excel::Writer::XLSX

