#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd "${DIR}"
if [ $1 == "mac" ]
then
echo Installing dependencies \for $1
	echo Checking \if \"pip\" \command is installed	
	if ! type pip 2>/dev/null; then
		sudo easy_install pip
	else
		echo pip is \installed
		echo upgrading pip
		sudo pip install --upgrade pip
	fi

	echo Checking \if python setuptools package is installed
	python -c "import setuptools" 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo pip install -U pip setuptools
	else
		echo python setuptools package is installed	
	fi
	
	echo Checking \if python tkinter package is installed	
	python -c "import Tkinter" 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		conda install -c anaconda tk
	else
		echo python tkinter package is installed	
	fi

	echo Checking \if python yaml package is installed
	python -c "import yaml" 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo pip install pyyaml
	else
		echo python yaml package is installed	
	fi
	
	echo Checking \if cutadapt is installed	
	if ! type cutadapt 2>/dev/null; then
		yes Y |sudo apt install python-cutadapt
		sudo pip install cutadapt
	else
		echo cutadapt is installed	
	fi

	echo Checking \if perldoc package is installed
	perldoc 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
	#if ! type perldoc 2>/dev/null; then
		sudo cpan perldoc
	else
		echo perldoc package is installed	
	fi
	
	echo Checking \if perl Try::Tiny::Retry module is installed	
	perldoc -l Try::Tiny::Retry 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo cpan Try::Tiny::Retry
	else
		echo perl Try::Tiny::Retry module is installed	
	fi

	echo Checking \if perl Excel::Writer::XLSX module is installed	
	perldoc -l Excel::Writer::XLSX 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo cpan Excel::Writer::XLSX
	else
		echo perl Excel::Writer::XLSX module is installed	
	fi

	echo Checking \if perl Time::Piece module is installed	
	perldoc -l Time::Piece 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo cpan Time::Piece
	else
		echo perl Time::Piece module is installed	
	fi
	
	echo Checking \if perl XML::Simple module is installed	
	perldoc -l XML::Simple 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo cpan XML::Simple
	else
		echo perl XML::Simple module is installed	
	fi
	
	echo Checking \if perl Digest::MD5 module is installed	
	perldoc -l Digest::MD5 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo cpan Digest::MD5
	else
		echo perl Digest::MD5 module is installed	
	fi

	echo Checking \if perl Bio::Perl module is installed
	perldoc -l Bio::Perl 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo cpan Bio::Perl
	else
		echo perl Bio::Perl module is installed	
	fi
	echo Checking \if git command is installed	
	if ! type git 2>/dev/null; then
		#Alternatively you can download it from here: http://git-scm.com/download/mac.
		ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
		brew doctor
		brew install git
	else
		echo git command is installed	
	fi

	echo Checking prokka on this PATH
	cd $2
	echo $2
	if [ -f "prokka/bin/prokka" ];then
		echo "Prokka exist."
		echo "Building database"
		sudo prokka/bin/prokka --setupdb
	else
		echo "Prokka does not exist, Downloading" >&2
		git clone https://github.com/tseemann/prokka.git
		echo "Building database"
		sudo prokka/bin/prokka --setupdb
	fi

	echo Checking \if tbl2asn command is installed
	if ! type tbl2asn 2>/dev/null; then
		curl -O ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/mac.tbl2asn.gz
		gunzip mac.tbl2asn.gz
		mv mac.tbl2asn tbl2asn
		chmod +x tbl2asn
		sudo cp tbl2asn /usr/local/bin
	else
		echo tbl2asn command is installed
	fi

fi

if [ $1 == "unix" ]
then
	echo Installing dependencies \for $1
	echo Checking \if \"pip\" \command is installed	
	if ! type pip 2>/dev/null; then
		yes Y | sudo apt install python-pip
	else
		echo pip is \installed
		echo upgrading pip
		sudo pip install --upgrade pip
	fi
	
	echo Checking \if python setuptools package is installed	
	python -c "import setuptools" 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		sudo apt-get install python-setuptools
	else
		echo python setuptools package is installed	
	fi
	
	echo Checking \if python tkinter package is installed	
	python -c "import Tkinter" 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes Y |sudo apt-get install python-tk
	else
		echo python tkinter package is installed	
	fi
	
	echo Checking \if python yaml package is installed	
	python -c "import yaml" 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes Y |sudo apt-get install python-yaml
	else
		echo python yaml package is installed	
	fi
	
	echo Checking \if cutadapt is installed	
	if ! type cutadapt 2>/dev/null; then
		yes Y |sudo apt install python-cutadapt
		sudo pip install cutadapt
	else
		echo cutadapt is installed	
	fi
	
	echo Checking \if perldoc package is installed	
	perldoc 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
	#if ! type perldoc 2>/dev/null; then
		yes Y |sudo apt-get install perl-doc
	else
		echo perldoc package is installed	
	fi
	
	echo Checking \if perl Try::Tiny::Retry module is installed	
	perldoc -l Try::Tiny::Retry 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo cpan Try::Tiny::Retry
	else
		echo perl Try::Tiny::Retry module is installed	
	fi
	
	echo Checking \if perl Excel::Writer::XLSX module is installed	
	perldoc -l Excel::Writer::XLSX 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo cpan Excel::Writer::XLSX
	else
		echo perl Excel::Writer::XLSX module is installed	
	fi

	echo Checking \if perl Time::Piece module is installed	
	perldoc -l Time::Piece 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo apt-get install libdatetime-perl
	else
		echo perl Time::Piece module is installed	
	fi
	
	echo Checking \if perl XML::Simple module is installed	
	perldoc -l XML::Simple 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo apt-get install libxml-simple-perl
	else
		echo perl XML::Simple module is installed	
	fi
	
	echo Checking \if perl Digest::MD5 module is installed	
	perldoc -l Digest::MD5 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo apt-get install libdigest-md5-perl
	else
		echo perl Digest::MD5 module is installed	
	fi
	
	echo Checking \if perl Bio::Perl module is installed	
	perldoc -l Bio::Perl 2>/dev/null
	if [[ $(echo $? ) == "1" ]]; then
		yes | sudo apt-get install git default-jre bioperl
	else
		echo perl Bio::Perl module is installed	
	fi
	
	echo Checking \if git command is installed	
	if ! type git 2>/dev/null; then
		yes Y | sudo apt install git
	else
		echo git command is installed	
	fi
	echo Checking prokka on this PATH
	cd $2
	echo $2
	if [ -f "prokka/bin/prokka" ];then
		echo "Prokka exist."
		echo "Building database"
		sudo prokka/bin/prokka --setupdb
	else
		echo "Prokka does not exist, Downloading" >&2
		git clone https://github.com/tseemann/prokka.git
		echo "Building database"
		sudo prokka/bin/prokka --setupdb
	fi
	#
	echo Checking \if tbl2asn command is installed
	if ! type tbl2asn 2>/dev/null; then
		curl -O ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz
		gunzip linux64.tbl2asn.gz
		mv linux64.tbl2asn tbl2asn
		chmod +x tbl2asn
		sudo cp tbl2asn /usr/local/bin
	else
		echo tbl2asn command is installed
	fi
	
fi

