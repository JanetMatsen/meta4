#!/bin/bash

# http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/bbmap/readme.txt 

bbmap.sh ref=isolate_genomes.fasta | tee index_std_out.txt
# note that the tee step didn't work.  Is the text printed to the terminal standard error?  Or some other stream? 
