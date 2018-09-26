FROM ummidock/fastqc:0.11.7-1

MAINTAINER Mohamed Mysara <mohamed.mysara@sckcen.be>

FROM ubuntu:16.04

WORKDIR /NGStools

RUN apt-get update && apt-get -y install \
	bash \
	python \
	wget

# Prokka
RUN apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl default-jre bioperl && git clone https://github.com/tseemann/prokka.git && prokka/bin/prokka --setupdb

# tbl2asn
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz && \
    gunzip linux.tbl2asn.gz && \
    mv linux.tbl2asn tbl2asn && \
    chmod +x tbl2asn && \
    mv tbl2asn /NGStools/prokka/binaries/linux/

ENV PATH="/NGStools/gff3toembl/gff3toembl/scripts:/NGStools/ncbi-blast-2.7.1+/bin:/NGStools/cd-hit-v4.6.8-2017-0621:/NGStools/cd-hit-v4.6.8-2017-0621/cd-hit-auxtools:/NGStools/exonerate-2.2.0-x86_64/bin:/NGStools/prokka/bin:${PATH}" LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH" PYTHONPATH="/NGStools/gff3toembl/gff3toembl:$PYTHONPATH"

# Install bacpipe
WORKDIR /NGStools/BacPipe

RUN wget https://www.dropbox.com/s/am32gc7u49jigg1/BacPipe.v1.7.unix.run?dl=0

RUN mv BacPipe.v1.7.unix.run?dl=0 BacPipe.v1.7.unix.run
RUN chmod +x BacPipe.v1.7.unix.run

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

