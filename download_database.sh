#!/bin/bash

#
# Download BLAST database.
#

[ -d $BLASTDB ] || mkdir -p $BLASTDB
cd $BLASTDB

wget http://skoda.projekty.ms.mff.cuni.cz/www/conservation/2019-09-18-swissprot.tar.gz
tar -xzvf 2019-09-18-swissprot.tar.gz

wget http://skoda.projekty.ms.mff.cuni.cz/www/conservation/2019-09-18-trembl.tar.gz
tar -xzvf 2019-09-18-trembl.tar.gz

wget http://skoda.projekty.ms.mff.cuni.cz/www/conservation/2019-09-18-uniref90.tar.gz
tar -xzvf 2019-09-18-uniref90.tar.gz

rm *.tar.gz
