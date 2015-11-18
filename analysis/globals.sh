#!/bin/bash

# to be sourced
DB=meta4
HOST=`mysql_host`

# pathes
BWA=/gscratch/esci/dacb/bwa/bin/bwa
SAMTOOLS=/gscratch/esci/dacb/samtools/bin/samtools
MYSQL=mysql
MYSQL_HOST=/gscratch/esci/dacb/mysql/bin/mysql_host
HTSEQ_COUNT=/gscratch/esci/dacb/htseq/bin/htseq-count

# job parameters
QL=walltime=999:99:99,mem=2gb,feature=8core
GROUP_LIST=hyak-lidstrom
EMAIL=dacb@uw.edu,jmatsen@uw.edu
