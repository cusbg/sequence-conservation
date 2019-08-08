#!/bin/bash

#
# Prepare BLAST database.
#

mkdir $BLASTDB
cd $BLASTDB

curl ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz | gunzip | /opt/ncbi-blast-2.9.0+/bin/makeblastdb -out uniref90 -dbtype prot -title UniRef90 -parse_seqids
curl ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz | gunzip | /opt/ncbi-blast-2.9.0+/bin/makeblastdb -out swissprot -dbtype prot -title SwissProt -parse_seqids
curl ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz | gunzip | /opt/ncbi-blast-2.9.0+/bin/makeblastdb -out trembl -dbtype prot -title TrEMBL -parse_seqids
