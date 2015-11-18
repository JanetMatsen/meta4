#!/bin/bash

#PBS -N "LakWasMet55_HOW8_2_aln"
#PBS -d /gscratch/lidstrom/meta_J_test_151111/split_read_investigation
#PBS -l walltime=999:99:99,mem=2gb,feature=8core,nodes=1:ppn=8
#PBS -W group_list=hyak-lidstrom
#PBS -M jmatsen@uw.edu #dacb@uw.edu,
#PBS -m abe
#PBS -o insert_sizes.log
#PBS -e insert_sizes.err

# based on:  http://www.genescripts.com/calculating-insert-size.php, http://genescripts.com/estimating-insert-size-from-bam-file.php
samtools view -f66 ../LakWasMet55_HOW8_2.sorted.bam|cut -f 9|sed 's/^-//' > insert_sizes.txt
# -f66 grabs only the properly paired reads. And only the 1st of each pair so you don't double count. 
# cut -f 9 grabs the 9th column, which is the insert length
