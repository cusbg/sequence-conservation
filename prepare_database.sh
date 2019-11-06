#!/bin/bash

#
# Prepare BLAST database.
#

[ -d $BLASTDB ] || mkdir -p $BLASTDB
cd $BLASTDB

export BLASTCMD="/opt/conservation/ncbi-blast-2.9.0+/bin/makeblastdb"

curl ftp://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz | gunzip | $BLASTCMD -out uniref90 -dbtype prot -title UniRef90 -parse_seqids

curl ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz | gunzip | $BLASTCMD -out swissprot -dbtype prot -title SwissProt -parse_seqids

curl ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz | gunzip | $BLASTCMD -out trembl -dbtype prot -title TrEMBL -parse_seqids
