# made by copying an existing dir:
    cp -r LakWasMeta9_HOW4_2 bwa_failure_control
# copied aln script off muffin: 
jmatsen@muffin:/dacb/meta4$ cat workspace/LakWasMeta9_HOW4_2/bwa/aln
#!/bin/bash

#PBS -N "LakWasMeta9_HOW4_2_aln"
#PBS -d /gscratch/lidstrom/meta4/workspace/LakWasMeta9_HOW4_2/bwa
#PBS -l walltime=999:99:99,mem=2gb,feature=8core,nodes=1:ppn=8
#PBS -W umask=002
#PBS -W group_list=hyak-lidstrom
#PBS -M dacb@uw.edu,jmatsen@uw.edu
#PBS -m abe
#PBS -o aln.log
#PBS -e aln.err

/gscratch/lidstrom/software/bwa/bin/bwa mem -p -M -t 8 /gscratch/lidstrom/meta4/data/isolate_genomes.fasta /gscratch/lidstrom/meta4/data/LakWasMeta9_HOW4_2/Raw_Data/8904.1.115212.GATCAG.fastq.gz > LakWasMeta9_HOW4_2.sam    

