#!/bin/bash

# to be sourced
DB=meta4
HOST=`mysql_host`

# pathes
BWA=/gscratch/lidstrom/software/bwa/bin/bwa
SAMTOOLS=/gscratch/lidstrom/software/samtools/bin/samtools
MYSQL=mysql
MYSQL_HOST=/gscratch/esci/dacb/mysql/bin/mysql_host
HTSEQ_COUNT=/gscratch/lidstrom/software/htseq/bin/htseq-count

# job parameters
QL=walltime=999:99:99,mem=2gb,feature=8core
GROUP_LIST=hyak-lidstrom
EMAIL=dacb@uw.edu,jmatsen@uw.edu
