
#!/bin/env bash
#

PERLBREW='http://install.perlbrew.pl'
BLASTLINUX='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-x64-linux.tar.gz'
#BLASTMAC='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-universal-macosx.tar.gz'
BLASTMAC='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-macosx.tar.gz'
BLASTFOLDER=blast

# PerlBrew needs to be installed to manage isolated perl environemnts if missing
command -v perlbrew >/dev/null 2>&1 || {
    echo 'Installing Perl brew...'
    echo "OS type: ${OSTYPE}"

    if [[ "${OSTYPE}" == 'linux'* ]]; then
        wget -O - http://install.perlbrew.pl | bash
    else
        curl -L ${PERLBREW} | bash
    fi

    source ~/perl5/perlbrew/etc/bashrc
    echo 'source ~/perl5/perlbrew/etc/bashrc' >> ~/.bash_profile
}

perlbrew init
echo "Do you want to install a local perl? [Y]/[N]"
read answer
if  [ $answer == 'Y' ]; then
    echo 'Installing perl-5.10.0 ...';
    perlbrew install -nf perl-5.10.0
else
    echo 'Local perl will be used ...';
fi

echo "Do you want to install a cpanmin locally through perlbrew? [Y]/[N]"
read answer
if  [ $answer == 'Y' ]; then
    echo 'Installing perlbrew install-cpanm...';
    perlbrew install-cpanm
else
    echo "Do you want to install a cpanmin as sudo[Y]/[N]"
    read answer
    if  [ $answer == 'Y' ]; then
        echo 'Installing cpanmin as sudo...';
        if [[ "${OSTYPE}" == 'linux'* ]]; then
            wget -O - https://cpanmin.us | perl - --sudo App::cpanminus
        else
            curl -L https://cpanmin.us | perl - --sudo App::cpanminus
        fi
    else
        echo 'Assuming cpanm is already installed...'
    fi
fi


# Installing NCBI Blast tools if missing
command -v blastall >/dev/null 2>&1 || {
    echo 'Installing Blast tools...'

    if [[ "${OSTYPE}" == 'linux'* ]]; then
        # TODO rename output folder
        wget ${BLASTLINUX}
    else
        # TODO Include versions for all OS: BSD, etc...
        curl ${BLASTMAC} > ${BLASTFOLDER}.tar.gz
        tar -zxvf ${BLASTFOLDER}.tar.gz
    fi

}
